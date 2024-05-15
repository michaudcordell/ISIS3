/** This is free and unencumbered software released into the public domain.

The authors of ISIS do not claim copyright on the contents of this file.
For more details about the LICENSE terms and the AUTHORS, you will
find files of those names at the top level of this repository. **/

/* SPDX-License-Identifier: CC0-1.0 */
#include "lronaccal.h"

#include <fstream>
#include <vector>

#include <QTextStream>
#include <QDir>
#include <QRegExp>
#include <QString>

#include "Brick.h"
#include "Camera.h"
#include "IException.h"
#include "iTime.h"
#include "Message.h"
#include "ProcessByLine.h"
#include "PvlGroup.h"
#include "SpecialPixel.h"
#include "Statistics.h"
#include "Table.h"
#include "TextFile.h"
#include "UserInterface.h"

namespace Isis {
  static constexpr int MASKED_PIXEL_VALUES = 8;
  static constexpr double MAXNONLIN = 600;

  // Anonymous namespace for internal linkage of DarkFileComparison
  namespace {
    /**
     * DarkFileInfo comparison object.
     *
     * Used for sorting DarkFileInfo objects. Sort first by difference from NAC time
     *
     */
    struct DarkFileComparison {
      int nacTime;

      DarkFileComparison(int nacTime)
      {
        this->nacTime = nacTime;
      }

      // sort dark files by distance from NAC time
      bool operator() (int A, int B) {
        return (abs(nacTime - A) < abs(nacTime - B));
      }
    };
  }

  static void CopyCubeIntoVector(QString &fileString, std::vector<double> &data);
  static void ReadTextDataFile(QString &fileString, std::vector<double> &data);
  static void ReadTextDataFile(QString &fileString, std::vector<std::vector<double> > &data);
  static void GetNearestDarkFile(double imgTime, 
                                 const QString &fileString, 
                                 QString &file, 
                                 std::vector<double> &avgDarkLineCube);
  static void GetNearestDarkFilePair(double imgTime, 
                                     const QString &fileString, 
                                     QString &file0, 
                                     QString &file1, 
                                     std::vector<double> &darkTimes, 
                                     bool &nearestDark, 
                                     bool &nearestDarkPair,
                                     std::vector<double> &avgDarkLineCube0,
                                     std::vector<double> &avgDarkLineCube1);
  static void GetCalibrationDirectory(QString calibrationType, QString &calibrationDirectory);
  static void GetWeightedDarkAverages(double imgTime, 
                                      const std::vector<double> &darkTimes, 
                                      std::vector<double> &weightedDarkTimeAvgs);
  static bool AllowedSpecialPixelType(double pixelValue);

  static void RemoveMaskedOffset(Buffer &inout, 
                                 bool summed,
                                 bool maskedLeftOnly,
                                 const std::vector<int> &maskedPixelsLeft, 
                                 const std::vector<int> &maskedPixelsRight);
  static void CorrectDark(Buffer &inout, 
                          bool nearestDarkPair, 
                          const std::vector<double> &avgDarkLineCube0, 
                          const std::vector<double> &avgDarkLineCube1, 
                          const std::vector<double> &weightedDarkTimeAvgs);
  static void CorrectNonlinearity(Buffer &inout, 
                                  const std::vector<double> &linearOffsetLine, 
                                  const std::vector<std::vector<double>> &linearityCoefficients,
                                  bool summed);
  static void CorrectFlatfield(Buffer &inout, const std::vector<double> &flatfieldLine);
  static void RadiometricCalibration(Buffer &inout, 
                                     double exposure, 
                                     bool iof, 
                                     bool isLeftNac, 
                                     double solarDistance,
                                     double iofLeft, 
                                     double iofRight, 
                                     double radianceLeft, 
                                     double radianceRight);

  /**
   * @brief Calling method of the application
   *
   * Performs radiometric corrections to images acquired by the Narrow Angle
   * Camera aboard the Lunar Reconnaissance Orbiter spacecraft.
   *
   * @param ui The user interface to parse the parameters from. 
   */
  void lronaccal(UserInterface &ui) {
    Cube iCube(ui.GetCubeName("FROM"));
    lronaccal(&iCube, ui);
  }

  /**
   * This is the main constructor lronaccal method. Lronaccal is used to calibrate LRO images.
   * 
   * @internal
   *   @history 2020-01-06 Victor Silva - Added option for base calibration directory
   *   @history 2020-07-19 Victor Silva - Updated dark calibration to use dark file option
   *                                      custom, nearest dark, or nearest dark pair.
   *   @history 2021-01-09 Victor Silva - Added code to check for exp_code = zero and if so
   *                                      then use only exp_code_zero dark files for dark calibration
   *  @history 2022-04-18 Victor Silva - Refactored to make callable for GTest framework
   *  @history 2024-05-15 Cordell Michaud - Refactored to not use global variables
   */
  void lronaccal(Cube *iCube, UserInterface &ui) {
    constexpr int LINE_SIZE = 5064;
    constexpr double KM_PER_AU = 149597871;

    double exposure = 1.0; // Exposure duration
    double solarDistance = 1.01; // average distance in [AU]
    std::vector<int> maskedPixelsLeft;
    std::vector<int> maskedPixelsRight;
    double radianceLeft = 1.0;
    double radianceRight = 1.0;
    double iofLeft = 1.0;
    double iofRight = 1.0;
    bool summed = true;
    bool masked = true;
    bool dark = true;
    bool nonlinear = true;
    bool flatfield = true;
    bool radiometric = true;
    bool iof = true;
    bool isLeftNac = true;
    bool maskedLeftOnly = false;
    bool nearestDarkPair = false;
    bool nearestDark = false;
    bool customDark = false;
    std::vector<double> avgDarkLineCube0;
    std::vector<double> avgDarkLineCube1;
    std::vector<double> linearOffsetLine;
    std::vector<double> darkTimes;
    std::vector<double> weightedDarkTimeAvgs;
    std::vector<double> flatfieldLine;
    std::vector<std::vector<double>> linearityCoefficients;
    double imgTime = 0.0;

    // We will be processing by line
    ProcessByLine p;

    masked = ui.GetBoolean("MASKED");
    dark = ui.GetBoolean("DARK");
    nonlinear = ui.GetBoolean("NONLINEARITY");
    flatfield = ui.GetBoolean("FLATFIELD");
    radiometric = ui.GetBoolean("RADIOMETRIC");
    iof = (ui.GetString("RADIOMETRICTYPE") == "IOF");

    Isis::Pvl lab(ui.GetCubeName("FROM"));
    Isis::PvlGroup &inst = lab.findGroup("Instrument", Pvl::Traverse);

    // Check if it is a NAC image
    QString instId = (QString) inst["InstrumentId"];
    instId = instId.toUpper();
    if (instId != "NACL" && instId != "NACR") {
      QString msg = "This is not a NAC image.  lrocnaccal requires a NAC image.";
      throw IException(IException::User, msg, _FILEINFO_);
    }

    // And check if it has already run through calibration
    if (lab.findObject("IsisCube").hasGroup("Radiometry")) {
      QString msg = "This image has already been calibrated";
      throw IException(IException::User, msg, _FILEINFO_);
    }

    if (lab.findObject("IsisCube").hasGroup("AlphaCube")) {
      QString msg = "This application can not be run on any image that has been geometrically transformed (i.e. scaled, rotated, sheared, or reflected) or cropped.";
      throw IException(IException::User, msg, _FILEINFO_);
    }

    if (instId == "NACL") {
      isLeftNac = true;
    }
    else {
      isLeftNac = false;
    }

    if (static_cast<int>(inst["SpatialSumming"]) == 1) {
      summed = false;
    }
    else {
      summed = true;
    }

    exposure = inst["LineExposureDuration"];

    p.SetInputCube(iCube, OneBand);

    // If there is any pixel in the image with a DN > 1000
    //  then the "left" masked pixels are likely wiped out and useless
    if (iCube->statistics()->Maximum() > 1000) {
      maskedLeftOnly = true;
    }

    QString flatFile, offsetFile, coefficientFile;

    if (masked) {
      QString maskedFile = ui.GetAsString("MASKEDFILE");
      if (maskedFile.toLower() == "default" || maskedFile.length() == 0) {
        GetCalibrationDirectory("", maskedFile);
        maskedFile = maskedFile + instId + "_MaskedPixels.????.pvl";
      }
      FileName maskedFileName(maskedFile);
      if (maskedFileName.isVersioned()) {
        maskedFileName = maskedFileName.highestVersion();
      }
      if (!maskedFileName.fileExists()) {
        QString msg = maskedFile + " does not exist.";
        throw IException(IException::User, msg, _FILEINFO_);
      }
      Pvl maskedPvl(maskedFileName.expanded());
      PvlKeyword maskedPixels;
      int cutoff;
      if (summed) {
        maskedPixels = maskedPvl["Summed"];
        cutoff = LINE_SIZE / 4;
      }
      else {
        maskedPixels = maskedPvl["FullResolution"];
        cutoff = LINE_SIZE / 2;
      }

      for (int i = 0; i < maskedPixels.size(); i++) {
        if ((isLeftNac && toInt(maskedPixels[i]) < cutoff) 
            || (!isLeftNac && toInt(maskedPixels[i]) > cutoff)) {
          maskedPixelsLeft.push_back(toInt(maskedPixels[i]));
        }
        else {
          maskedPixelsRight.push_back(toInt(maskedPixels[i]));
        }
      }
    }
    
    std::vector<QString> darkFiles;

    if (dark) {
      QString darkFileType = ui.GetString("DARKFILETYPE");
      darkFileType = darkFileType.toUpper();
      if (darkFileType == "CUSTOM") {
        customDark = true;
        ui.GetAsString("DARKFILE", darkFiles);
      }
      else if (darkFileType == "PAIR" || darkFileType == "") {
        nearestDarkPair = true;
      }
      else if (darkFileType == "NEAREST") {
        nearestDark = true;
      }
      else {
        QString msg = "Error: Dark File Type selection failed.";
        throw IException(IException::User, msg, _FILEINFO_);
      }
      //Options are NEAREST, PAIR, and CUSTOM
      if (customDark) {
        if (darkFiles.size() == 1 && darkFiles[0] != "") {
          CopyCubeIntoVector(darkFiles[0], avgDarkLineCube0);
        }
        else {
          QString msg = "Custom dark file not provided. Please provide file or choose another option.";
          throw IException(IException::User, msg, _FILEINFO_);
        }
      }
      else {
        QString darkFile;
        imgTime = iTime(inst["StartTime"][0]).Et();
        GetCalibrationDirectory("nac_darks", darkFile);
        darkFile = darkFile + instId + "_AverageDarks_*T";
        
        if (summed) {
          darkFile += "_Summed";
        }
        // use exp0 dark files if cube's exp_code=0
        Isis::PvlGroup &pvl_archive_group = lab.findGroup("Archive", Pvl::Traverse);
        if (static_cast<int>(pvl_archive_group["LineExposureCode"]) == 0) {
          darkFile += "_exp0";
        }

        darkFile += ".????.cub";

        if (nearestDark) {
          darkFiles.resize(1);
          GetNearestDarkFile(imgTime, darkFile, darkFiles[0], avgDarkLineCube0);
        }
        else {
          darkFiles.resize(2);
          GetNearestDarkFilePair(imgTime, darkFile, darkFiles[0], darkFiles[1], darkTimes, 
                                 nearestDark, nearestDarkPair, avgDarkLineCube0, 
                                 avgDarkLineCube1);
          //get weigted time avgs
          if (darkTimes.size() == 2) {
            GetWeightedDarkAverages(imgTime, darkTimes, weightedDarkTimeAvgs);
          }
        }
      }
    }

    if (nonlinear) {
      offsetFile = ui.GetAsString("OFFSETFILE");

      if (offsetFile.toLower() == "default" || offsetFile.length() == 0) {
        GetCalibrationDirectory("", offsetFile);
        offsetFile = offsetFile + instId + "_LinearizationOffsets";
        if (summed) {
          offsetFile += "_Summed";
        }
        offsetFile += ".????.cub";
      }
      CopyCubeIntoVector(offsetFile, linearOffsetLine);
      coefficientFile = ui.GetAsString("NONLINEARITYFILE");
      if (coefficientFile.toLower() == "default" || coefficientFile.length() == 0) {
        GetCalibrationDirectory("", coefficientFile);
        coefficientFile = coefficientFile + instId + "_LinearizationCoefficients.????.txt";
      }
      ReadTextDataFile(coefficientFile, linearityCoefficients);
    }

    if (flatfield) {
      flatFile = ui.GetAsString("FLATFIELDFILE");

      if (flatFile.toLower() == "default" || flatFile.length() == 0) {
        GetCalibrationDirectory("", flatFile);
        flatFile = flatFile + instId + "_Flatfield";
        if (summed) {
          flatFile += "_Summed";
        }
        flatFile += ".????.cub";
      }
      CopyCubeIntoVector(flatFile, flatfieldLine);
    }

    if (radiometric) {
      QString radFile = ui.GetAsString("RADIOMETRICFILE");

      if (radFile.toLower() == "default" || radFile.length() == 0) {
        GetCalibrationDirectory("", radFile);
        radFile = radFile + "NAC_RadiometricResponsivity.????.pvl";
      }

      FileName radFileName(radFile);
      if (radFileName.isVersioned()) {
        radFileName = radFileName.highestVersion();
      }
      if (!radFileName.fileExists()) {
        QString msg = radFile + " does not exist.";
        throw IException(IException::User, msg, _FILEINFO_);
      }

      Pvl radPvl(radFileName.expanded());

      if (iof) {
        iTime startTime((QString) inst["StartTime"]);

        try {
          Camera *cam;
          cam = iCube->camera();
          cam->setTime(startTime);
          solarDistance = cam->sunToBodyDist() / KM_PER_AU;
        }
        catch(IException &e) {
          // Failed to instantiate a camera, try furnishing kernels directly
          try {
            double etStart = startTime.Et();
            // Get the distance between the Moon and the Sun at the given time in
            // Astronomical Units (AU)
            QString bspKernel1 = p.MissionData("lro", "/kernels/tspk/moon_pa_de421_1900-2050.bpc", false);
            QString bspKernel2 = p.MissionData("lro", "/kernels/tspk/de421.bsp", false);
            furnsh_c(bspKernel1.toLatin1().data());
            furnsh_c(bspKernel2.toLatin1().data());
            QString pckKernel1 = p.MissionData("base", "/kernels/pck/pck?????.tpc", true);
            QString pckKernel2 = p.MissionData("lro", "/kernels/pck/moon_080317.tf", false);
            QString pckKernel3 = p.MissionData("lro", "/kernels/pck/moon_assoc_me.tf", false);
            furnsh_c(pckKernel1.toLatin1().data());
            furnsh_c(pckKernel2.toLatin1().data());
            furnsh_c(pckKernel3.toLatin1().data());
            double sunpos[6], lt;
            spkezr_c("sun", etStart, "MOON_ME", "LT+S", "MOON", sunpos, &lt);
            solarDistance = vnorm_c(sunpos) / KM_PER_AU;
            unload_c(bspKernel1.toLatin1().data());
            unload_c(bspKernel2.toLatin1().data());
            unload_c(pckKernel1.toLatin1().data());
            unload_c(pckKernel2.toLatin1().data());
            unload_c(pckKernel3.toLatin1().data());
          }
          catch(IException &e) {
            QString msg = "Unable to find the necessary SPICE kernels for converting to IOF";
            throw IException(e, IException::User, msg, _FILEINFO_);
          }
        }
        iofLeft = radPvl["IOF_LEFT"];
        iofRight = radPvl["IOF_RIGHT"];
      }
      else {
        radianceLeft = radPvl["Radiance_LEFT"];
        radianceRight = radPvl["Radiance_RIGHT"];
      }
    }

    /**
     * @brief This method processes buffer by line to calibrate.
     *
     * @param[in] in Buffer to hold 1 line of cube data
     * @param[out] out Buffer to hold 1 line of cube data
     *
     */
    auto Calibrate = [masked, dark, nonlinear, flatfield, radiometric, summed, 
                      maskedLeftOnly, nearestDarkPair, iof, isLeftNac, exposure, 
                      solarDistance, iofLeft, iofRight, radianceLeft, radianceRight, 
                      &maskedPixelsLeft, &maskedPixelsRight, avgDarkLineCube0, 
                      &avgDarkLineCube1, &weightedDarkTimeAvgs, &linearOffsetLine, 
                      &linearityCoefficients, &flatfieldLine](Buffer &in, Buffer &out) -> void {
      for (int i = 0; i < in.size(); i++) {
        out[i] = in[i];
      }

      if (masked) {
        RemoveMaskedOffset(out, summed, maskedLeftOnly, maskedPixelsLeft, maskedPixelsRight);
      }

      if (dark) {
        CorrectDark(out, nearestDarkPair, avgDarkLineCube0, avgDarkLineCube1, weightedDarkTimeAvgs);
      }

      if (nonlinear) {
        CorrectNonlinearity(out, linearOffsetLine, linearityCoefficients, summed);
      }

      if (flatfield) {
        CorrectFlatfield(out, flatfieldLine);
      }

      if (radiometric) {
        RadiometricCalibration(out, exposure, iof, isLeftNac, solarDistance, iofLeft, iofRight, 
          radianceLeft, radianceRight);
      }
    };

    // Setup the output cube
    Cube *oCube = p.SetOutputCube(ui.GetCubeName("TO"), ui.GetOutputAttribute("TO")); 
    // Start the line-by-line calibration sequence
    p.StartProcess(Calibrate);

    PvlGroup calgrp("Radiometry");
    if (masked) {
      PvlKeyword darkColumns("DarkColumns");
      for (std::size_t i = 0; i < maskedPixelsLeft.size(); i++) {
        darkColumns += toString(maskedPixelsLeft[i]);
      }
      for (std::size_t i = 0; i < maskedPixelsRight.size(); i++) {
        darkColumns += toString(maskedPixelsRight[i]);
      }
      calgrp += darkColumns;
    }

    if (dark) {
      PvlKeyword darks("DarkFiles");
      darks.addValue(darkFiles[0]);
      if (nearestDark) {
        calgrp += PvlKeyword("DarkFileType", "NearestDarkFile");
      }
      else if (nearestDarkPair) {
        calgrp += PvlKeyword("DarkFileType", "NearestDarkFilePair");
        darks.addValue(darkFiles[1]);
      }
      else {
        calgrp += PvlKeyword("DarkFileType", "CustomDarkFile");
      }

      calgrp += darks;
    }

    if (nonlinear) {
      calgrp += PvlKeyword("NonlinearOffset", offsetFile);
      calgrp += PvlKeyword("LinearizationCoefficients", coefficientFile);
    }

    if (flatfield) {
      calgrp += PvlKeyword("FlatFile", flatFile);
    }
    if (radiometric) {
      if (iof) {
        calgrp += PvlKeyword("RadiometricType", "IOF");
        if (isLeftNac) {
          calgrp += PvlKeyword("ResponsivityValue", toString(iofLeft));
        }
        else {
          calgrp += PvlKeyword("ResponsivityValue", toString(iofRight));
        }
      }
      else {
        calgrp += PvlKeyword("RadiometricType", "AbsoluteRadiance");
        if (isLeftNac) {
          calgrp += PvlKeyword("ResponsivityValue", toString(radianceLeft));
        }
        else {
          calgrp += PvlKeyword("ResponsivityValue", toString(radianceRight));
        }
      }
      calgrp += PvlKeyword("SolarDistance", toString(solarDistance));
    }

    oCube->putGroup(calgrp);
    p.EndProcess();
  }

  /**
   * @brief This method copies a cube into a vector.
   *
   * @param[in] fileString Cube file
   * @param[out] data vector of double
   *
   */
  static void CopyCubeIntoVector(QString &fileString, std::vector<double> &data) {
    Cube cube;
    FileName filename(fileString);
    if (filename.isVersioned()) {
      filename = filename.highestVersion();
    }
    if (!filename.fileExists()) {
      QString msg = fileString + " does not exist.";
      throw IException(IException::User, msg, _FILEINFO_);
    }
    cube.open(filename.expanded());
    Brick brick(cube.sampleCount(), cube.lineCount(), cube.bandCount(), cube.pixelType());
    brick.SetBasePosition(1, 1, 1);
    cube.read(brick);
    data.clear();
    for (int i = 0; i < cube.sampleCount(); i++) {
      data.push_back(brick[i]);
    }

    fileString = filename.original();

    if (data.empty()) {
      QString msg = "Copy from + " + fileString + " into vector failed.";
      throw IException(IException::User, msg, _FILEINFO_);
    }

  }

  /**
   * Read text data file - overloaded method.
   *
   * @param[in] fileString Text file
   * @param[out] data vector of double
   */
  static void ReadTextDataFile(QString &fileString, std::vector<double> &data) {
    FileName filename(fileString);
    if (filename.isVersioned()) {
      filename = filename.highestVersion();
    }
    if (!filename.fileExists()) {
      QString msg = fileString + " does not exist.";
      throw IException(IException::User, msg, _FILEINFO_);
    }
    TextFile file(filename.expanded());
    QString lineString;
    unsigned int line = 0;
    while (file.GetLine(lineString)) {
      data.push_back(toDouble(lineString.split(QRegExp("[ ,;]")).first()));
      line++;
    }
    fileString = filename.original();
  }

  /**
   * Read text data file - overloaded method.
   *
   * @param[in] fileString Text file
   * @param[out] data multi-dimensional vector of double
   */
  static void ReadTextDataFile(QString &fileString, std::vector<std::vector<double>> &data) {
    FileName filename(fileString);
    if (filename.isVersioned()) {
      filename = filename.highestVersion();
    }
    if (!filename.fileExists()) {
      QString msg = fileString + " does not exist.";
      throw IException(IException::User, msg, _FILEINFO_);
    }
    TextFile file(filename.expanded());
    QString lineString;
    while (file.GetLine(lineString)) {
      std::vector<double> line;
      lineString = lineString.simplified().remove(QRegExp("^[ ,]*")).trimmed();

      QStringList lineTokens = lineString.split(QRegExp("[ ,]"), Qt::SkipEmptyParts);
      foreach (QString value, lineTokens) {
        line.push_back(toDouble(value));
      }

      data.push_back(line);
    }

    fileString = filename.original();
  }

  /**
   * @brief Finds the best dark file for NAC calibration and copies it into a vector.
   *
   * GetNearestDarkFile will get the dark file with the closest time (before or after) to the 
   * image time to be used for calibration abnd copy it into a vector of doubles.
   *
   * @param[in] imgTime Image time
   * @param[in] fileString String pattern defining dark files to search
   * @param[out] file Filename of dark file
   * @param[out] avgDarkLineCube Average dark line cube data for dark file
   */
  static void GetNearestDarkFile(double imgTime, 
                                 const QString &fileString, 
                                 QString &file, 
                                 std::vector<double> &avgDarkLineCube) {
    FileName filename(fileString);
    QString basename = FileName(filename.baseName()).baseName(); // We do it twice to remove the ".????.cub"
    // create a regular expression to capture time from filenames
    QString regexPattern(basename);
    regexPattern.replace("*", "([0-9\\.-]*)");
    QRegExp regex(regexPattern);
    // create a filter for the QDir to only load files matching our name
    QString filter(basename);
    filter.append(".*");
    // get a list of dark files that match our basename
    QDir dir(filename.path(), filter);
    std::vector<int> matchedDarkTimes;
    matchedDarkTimes.reserve(dir.count());
    // Loop through all files in the dir that match our basename and extract time
    for (unsigned int i = 0; i < dir.count(); i++) {
      // match against our regular expression
      int pos = regex.indexIn(dir[i]);
      if (pos == -1) {
        continue; // filename did not match basename regex (time contain non-digit)
      }
      // Get a list of regex matches. Item 0 should be the full QString, item 1 is time.
      QStringList texts = regex.capturedTexts();
      if (texts.size() < 1) {
        continue; // could not find time
      }
      // extract time from regex texts
      bool timeOK;
      int fileTime = texts[1].toInt(&timeOK);
      if (!timeOK) {
        continue; // time was not a valid numeric value
      }
      matchedDarkTimes.push_back(fileTime);
    }
    // sort the files by distance from nac time
    DarkFileComparison darkComp(static_cast<int>(imgTime));
    sort(matchedDarkTimes.begin(), matchedDarkTimes.end(), darkComp);
    int darkTime = matchedDarkTimes[0];
    int fileTimeIndex = fileString.indexOf("*T");
    file = fileString;
    file.replace(fileTimeIndex, 1, toString(darkTime));
    CopyCubeIntoVector(file, avgDarkLineCube);
  }

  /**
   * @brief Finds the best dark files for NAC calibration.
   *
   * GetNearestDarkFilePair will get the average between the two darks files
   * that the image lies between (time-wise) and copy them into vectors of doubles.
   * If this pair is not found, the nearest dark file will be used
   * for calibration.
   *
   * @param[in] imgTime Image time
   * @param[in] fileString String pattern defining dark files to search (ie. lro/calibration/nac_darks/NAC*_AverageDarks_*T_.????.cub)
   * @param[out] file0 Filename of dark file 1
   * @param[out] file1 Filename of dark file 2
   * @param[out] darkTimes Dark times
   * @param[out] nearestDark Nearest dark should be used
   * @param[out] nearestDarkPair Nearest dark pair should be used
   * @param[out] avgDarkLineCube0 Average dark line cube data for dark file 1
   * @param[out] avgDarkLineCube1 Average dark line cube data for dark file 2
   */
  static void GetNearestDarkFilePair(double imgTime, 
                                    const QString &fileString, 
                                    QString &file0, 
                                    QString &file1, 
                                    std::vector<double> &darkTimes, 
                                    bool &nearestDark, 
                                    bool &nearestDarkPair,
                                    std::vector<double> &avgDarkLineCube0,
                                    std::vector<double> &avgDarkLineCube1) {
    FileName filename(fileString);
    // We use baseName() twice to remove the ".????.cub"
    QString basename = FileName(filename.baseName()).baseName();
    // create a regular expression to capture time from filenames
    QString regexPattern(basename);
    regexPattern.replace("*", "([0-9\\.-]*)");
    QRegExp regex(regexPattern);
    // create a filter for the QDir to only load files matching our name
    QString filter(basename);
    filter.append(".*");
    // get a list of dark files that match our basename
    QDir dir( filename.path(), filter );
    std::vector<int> matchedDarkTimes;
    matchedDarkTimes.reserve(dir.count());
    if (dir.count() < 1) {
      QString msg = "Could not find any dark file of type " + filter + ".\n";
      throw IException(IException::User, msg, _FILEINFO_);
    }
    // Loop through all files in the dir that match our basename and extract time
    for (unsigned int i = 0; i < dir.count(); i++) {
      // match against our regular expression
      int pos = regex.indexIn(dir[i]);
      if (pos == -1) {
        continue; // filename did not match basename regex (time contain non-digit)
      }
      // Get a list of regex matches. Item 0 should be the full QString, item 1
      // is time.
      QStringList texts = regex.capturedTexts();
      if (texts.size() < 1) {
        continue; // could not find time
      }
      // extract time from regex texts
      bool timeOK;
      int fileTime = texts[1].toInt(&timeOK);
      if (!timeOK) {
        continue; // time was not a valid numeric value
      }
      matchedDarkTimes.push_back(fileTime);
    }
    // sort the files by distance from nac time
    DarkFileComparison darkComp(static_cast<int>(imgTime));
    sort(matchedDarkTimes.begin(), matchedDarkTimes.end(), darkComp);

    int fileTimeIndex = fileString.indexOf("*T");
    int t0 = 0;
    int t1 = 0;
    //Let's find the first time before the image
    for (size_t i = 0; i < matchedDarkTimes.size(); i++) {
      if (matchedDarkTimes[i] <= static_cast<int>(imgTime)) {
        t0 = matchedDarkTimes[i];
        break;
      }
    }
    //Let's find the second time
    for (size_t i = 0; i < matchedDarkTimes.size(); i++) {
      if (matchedDarkTimes[i] >= static_cast<int>(imgTime)) {
        t1 = matchedDarkTimes[i];
        break;
      }
    }
    if ((t0 && t1) && (t0!=t1)) {
      int timeDayDiff =  abs(t1 -t0)/86400.0;
      
    //check time range between darks is within 45 day window
    if (timeDayDiff < 0  || timeDayDiff > 45) {
        QString msg = "Could not find a pair of dark files within 45 day range that includes the image [" 
                    + basename + "]. Check to make sure your set of dark files is complete.\n";
        throw IException(IException::User, msg, _FILEINFO_);
      }
      else {
        file0 = fileString;
        file0.replace(fileTimeIndex, 1, toString(t0));
        CopyCubeIntoVector(file0, avgDarkLineCube0);
        darkTimes.push_back(t0);
        file1 = fileString;
        file1.replace(fileTimeIndex, 1, toString(t1));
        CopyCubeIntoVector(file1, avgDarkLineCube1);
        darkTimes.push_back(t1);
      }
    }
    else {
      nearestDark = true;
      nearestDarkPair = false;
      int darkTime = matchedDarkTimes[0];
      file0 = fileString;
      file0.replace(fileTimeIndex, 1, toString(darkTime));
      CopyCubeIntoVector(file0, avgDarkLineCube0);
      darkTimes.push_back(darkTime);
    }
  }

  /**
   * @brief This method returns a QString containing the path of an LRO calibration directory.
   *
   * @param[in] calibrationType
   * @param[out] calibrationDirectory Path of the calibration directory
   *
   * @internal
   *   @history 2020-01-06 Victor Silva - Added option for base calibration directory
   */
  static void GetCalibrationDirectory(QString calibrationType, QString &calibrationDirectory) {
    PvlGroup &dataDir = Preference::Preferences().findGroup("DataDirectory");
    QString missionDir = (QString) dataDir["LRO"];
    if (calibrationType != "") {
      calibrationType += "/";
    }
    calibrationDirectory = missionDir + "/calibration/" + calibrationType;
  }

  /**
   * @brief Get weighted time averages for calculating pixel dark average.
   *
   * @param[in] imgTime Image time
   * @param[in] darkTimes Dark times
   * @param[out] weightedDarkTimeAvgs Weighted time averages w0 and w1 for dark file
   *
   */
  static void GetWeightedDarkAverages(double imgTime, 
                                      const std::vector<double> &darkTimes, 
                                      std::vector<double> &weightedDarkTimeAvgs) {
    int iTime = static_cast<int>(imgTime);
    int t0 = 0;
    int t1 = 0;

    if (!darkTimes.empty()) {
      if (darkTimes.size() == 2) {
        t0 = darkTimes[0];
        t1 = darkTimes[1];
        double weight0 =
        (( t1!=iTime ) * ( (t1 > iTime ) * ( t1 - iTime) ))
        / (((( t1!=iTime ) * ( (t1 > iTime ) * ( t1 - iTime) )) +
          (( t0!=iTime ) * ( (t0 < iTime ) * ( iTime - t0) )) ) * 1.0);

        double weight1 = (( t0!=iTime ) * ( (t0 < iTime ) * ( iTime - t0) ))
        / (((( t1!=iTime ) * ( (t1 > iTime ) * ( t1 - iTime) )) +
          (( t0!=iTime ) * ( (t0 < iTime ) * ( iTime - t0) )) ) * 1.0);

        weightedDarkTimeAvgs.clear();
        weightedDarkTimeAvgs.push_back(weight0);
        weightedDarkTimeAvgs.push_back(weight1);
      }
    }
  }

  /**
   * @brief Allow special pixel types.
   *
   * @param pixelValue double
   *
   * @return bool
   */
  static bool AllowedSpecialPixelType(double pixelValue) {
    bool result = IsHisPixel(pixelValue) 
                || IsLisPixel(pixelValue) 
                || IsHrsPixel(pixelValue) 
                || IsLrsPixel(pixelValue);
    return result;
  }

  /**
   * @brief Remove masked offset in-place.
   *
   * @param[in,out] inout Buffer
   * @param[in] summed Summed
   * @param[in] maskedLeftOnly Only left is masked
   * @param[in] maskedPixelsLeft Left masked pixels
   * @param[in] maskedPixelsRight Right masked pixels
   */
  static void RemoveMaskedOffset(Buffer &inout, 
                                 bool summed,
                                 bool maskedLeftOnly,
                                 const std::vector<int> &maskedPixelsLeft, 
                                 const std::vector<int> &maskedPixelsRight) {
    int numMasked = MASKED_PIXEL_VALUES;
    if (summed) {
      numMasked /= 2;
    }

    std::vector<Statistics> statsLeft(numMasked, Statistics());
    std::vector<Statistics> statsRight(numMasked, Statistics());

    std::vector<int> leftRef(numMasked, 0);
    std::vector<int> rightRef(numMasked, 0);

    for (std::size_t i = 0; i < maskedPixelsLeft.size(); i++) {
      statsLeft[maskedPixelsLeft[i] % numMasked].AddData(&inout[maskedPixelsLeft[i]], 1);
      leftRef[maskedPixelsLeft[i] % numMasked] += maskedPixelsLeft[i];
    }

    for (std::size_t i = 0; i < maskedPixelsRight.size(); i++) {
      statsRight[maskedPixelsRight[i] % numMasked].AddData(&inout[maskedPixelsRight[i]], 1);
      rightRef[maskedPixelsRight[i] % numMasked] += maskedPixelsRight[i];
    }

    // left/rightRef is the center (average) of all the masked pixels in the set
    for (int i = 0; i < numMasked; i++) {
      leftRef[i] /= statsLeft[i].TotalPixels();
      rightRef[i] /= statsRight[i].TotalPixels();
    }

    if (maskedLeftOnly) {
      for (int i = 0; i < inout.size(); i++) {
        inout[i] -= statsLeft[i % numMasked].Average();
      }
    }
    else {
      // If we are using both sides, we interpolate between them

      for (int i = 0; i < inout.size(); i++) {
        inout[i] -= (statsLeft[i % numMasked].Average() * (rightRef[i % numMasked] - i) 
                     + statsRight[i % numMasked].Average()
                     * (i - leftRef[i % numMasked]))
                  / (rightRef[i % numMasked] - leftRef[i % numMasked]);
      }
    }
  }

  /**
   * @brief Performs dark correction of the pixel being processed in-place.
   * 
   * Dark correction will use the nearest dark pair if nearestDarkPair is true, otherwise it will 
   * use the nearest dark file.
   *
   * @param[in,out] inout Buffer
   * @param[in] nearestDarkPair Use nearest dark pair
   * @param[in] avgDarkLineCube0 Average dark line cube 1
   * @param[in] avgDarkLineCube1 Average dark line cube 2
   * @param[in] weightedDarkTimeAvgs Weighted dark time averages
   */
  static void CorrectDark(Buffer &inout, 
                          bool nearestDarkPair, 
                          const std::vector<double> &avgDarkLineCube0, 
                          const std::vector<double> &avgDarkLineCube1, 
                          const std::vector<double> &weightedDarkTimeAvgs) {
    for (int i = 0; i < inout.size(); i++) {
      if (nearestDarkPair &&
        (!IsSpecial(inout[i]) || AllowedSpecialPixelType(inout[i])) &&
        (!IsSpecial(avgDarkLineCube0[i]) || AllowedSpecialPixelType(avgDarkLineCube0[i])) &&
        (!IsSpecial(avgDarkLineCube1[i]) || AllowedSpecialPixelType(avgDarkLineCube1[i])) &&
        (!IsSpecial(inout[i]) || AllowedSpecialPixelType(inout[i])) ) {
        double w0 = weightedDarkTimeAvgs[0];
        double w1 = weightedDarkTimeAvgs[1];
        double pixelDarkAvg = (avgDarkLineCube0[i]*w0)+(avgDarkLineCube1[i]*w1);

        inout[i] -= pixelDarkAvg;

      } else if
        ((!IsSpecial(avgDarkLineCube0[i]) || AllowedSpecialPixelType(avgDarkLineCube0[i])) &&
        (!IsSpecial(inout[i]) || AllowedSpecialPixelType(inout[i])) ) {

        inout[i] -= avgDarkLineCube0[i];

      }
      else {
        inout[i] = Isis::Null;
      }
    }
  }

  /**
   * @brief Correct non-linearity of the pixel being processed in-place.
   *
   * @param[in,out] inout Buffer
   * @param[in] linearOffsetLine Linear offset line
   * @param[in] linearityCoefficients Linearity coefficients
   * @param[in] summed Summed
   */
  static void CorrectNonlinearity(Buffer &inout, 
                                  const std::vector<double> &linearOffsetLine, 
                                  const std::vector<std::vector<double>> &linearityCoefficients,
                                  bool summed) {
    for (int i = 0; i < inout.size(); i++) {
      if (!IsSpecial(inout[i])) {
        inout[i] += linearOffsetLine[i];

        if (inout[i] < MAXNONLIN) {
          if (summed) {
            inout[i] -= (1.0 / (linearityCoefficients[2* i ][0]
                                * pow(linearityCoefficients[2* i ][1], inout[i])
                                + linearityCoefficients[2* i ][2])
                         + 1.0 / (linearityCoefficients[2* i + 1][0] 
                                  * pow(linearityCoefficients[2* i + 1][1], inout[i])
                                  + linearityCoefficients[2* i + 1][2]))
                      / 2;
          }
          else {
            inout[i] -= 1.0 / (linearityCoefficients[i][0] 
                               * pow(linearityCoefficients[i][1], inout[i])
                               + linearityCoefficients[i][2]);
          }
        }
      }
      else
        inout[i] = Isis::Null;
    }
  }

  /**
   * @brief Perform flatfield correction of the pixel being processed in-place.
   * 
   * @param[in,out] inout Buffer
   * @param[in] flatfieldLine Flatfield line
   */
  static void CorrectFlatfield(Buffer &inout, const std::vector<double> &flatfieldLine) {
    for (int i = 0; i < inout.size(); i++) {
      if (!IsSpecial(inout[i]) && flatfieldLine[i] > 0) {
        inout[i] /= flatfieldLine[i];
      }
      else {
        inout[i] = Isis::Null;
      }
    }
  }

  /**
   * @brief Perform radiometric calibration of the pixel being processed in-place.
   *
   * @param[in,out] inout Buffer
   * @param[in] exposure Exposure
   * @param[in] iof IOF
   * @param[in] isLeftNac Is a left NAC
   * @param[in] solarDistance Average solar distance in AU
   * @param[in] iofLeft Left IOF
   * @param[in] iofRight Right IOF 
   * @param[in] radianceLeft Left radiance
   * @param[in] radianceRight Right radiance
   */
  static void RadiometricCalibration(Buffer &inout, 
                                     double exposure, 
                                     bool iof, 
                                     bool isLeftNac, 
                                     double solarDistance,
                                     double iofLeft, 
                                     double iofRight, 
                                     double radianceLeft, 
                                     double radianceRight) {
    for (int i = 0; i < inout.size(); i++) {
      if (!IsSpecial(inout[i])) {
        inout[i] /= exposure;
        if (iof) {
          if (isLeftNac) {
            inout[i] = inout[i] * pow(solarDistance, 2) / iofLeft;
          }
          else {
            inout[i] = inout[i] * pow(solarDistance, 2) / iofRight;
          }
        }
        else {
          if (isLeftNac) {
            inout[i] = inout[i] / radianceLeft;
          }
          else {
            inout[i] = inout[i] / radianceRight;
          }
        }
      }
      else {
        inout[i] = Isis::Null;
      }
    }
  }
}
