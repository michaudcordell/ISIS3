#include "Isis.h"

#include <set>
#include <sstream>

#include <QMap>
#include <QSet>
#include <QVector>

#include "Angle.h"
#include "Camera.h"
#include "CameraFactory.h"
#include "ControlMeasure.h"
#include "ControlNet.h"
#include "ControlPoint.h"
#include "Cube.h"
#include "CubeManager.h"
#include "FileList.h"
#include "IException.h"
#include "iString.h"
#include "Longitude.h"
#include "Latitude.h"
#include "Projection.h"
#include "ProjectionFactory.h"
#include "Pvl.h"
#include "SerialNumber.h"
#include "SurfacePoint.h"
#include "UserInterface.h"

using namespace std;
using namespace Isis;

void ExtractPointList(ControlNet &outNet, QVector<iString> nonListedPoints);
void ExtractLatLonRange(ControlNet &outNet, QVector<iString> nonLatLonPoints,
                        QVector<iString> cannotGenerateLatLonPoints,
                        QMap<iString, iString> sn2filename);
bool NotInLatLonRange(SurfacePoint surfacePt, Latitude minlat,
                      Latitude maxlat, Longitude minlon, Longitude maxlon);
void WriteCubeOutList(ControlNet cnet, QMap<iString, iString> sn2file);
void WriteResults(iString filename, QVector<iString> results);
void omit(ControlNet &cnet, int cp);
void omit(ControlPoint *point, int cm);


// Main program
void IsisMain() {
  UserInterface &ui = Application::GetUserInterface();

  if(!ui.WasEntered("FROMLIST") && ui.WasEntered("TOLIST")) {
    std::string msg = "To create a [TOLIST] the [FROMLIST] parameter must be provided.";
    throw IException(IException::User, msg, _FILEINFO_);
  }

  bool noIgnore          = ui.GetBoolean("NOIGNORE");
  bool noMeasureless     = ui.GetBoolean("NOMEASURELESS");
  bool noSingleMeasure   = ui.GetBoolean("NOSINGLEMEASURES");
  bool reference         = ui.GetBoolean("REFERENCE");
  bool fixed             = ui.GetBoolean("FIXED");
  bool noTolerancePoints = ui.GetBoolean("TOLERANCE");
  bool pointsEntered     = ui.WasEntered("POINTLIST");
  bool cubePoints        = ui.GetBoolean("CUBES");
  bool cubeMeasures      = ui.GetBoolean("CUBEMEASURES");
  bool retainReference   = ui.GetBoolean("RETAIN_REFERENCE");
  bool latLon            = ui.GetBoolean("LATLON");

  if(!(noIgnore || noMeasureless || noSingleMeasure || reference || fixed ||
       noTolerancePoints || pointsEntered || cubePoints || latLon)) {
    std::string msg = "At least one filter must be selected [";
    msg += "NOIGNORE,NOMEASURELESS,NOSINGLEMEASURE,REFERENCE,FIXED,TOLERANCE,";
    msg += "POINTLIST,CUBES,LATLON]";
    throw IException(IException::User, msg, _FILEINFO_);
  }

  if(cubeMeasures || reference) {
    noMeasureless = true;
  }

  // Gets the input parameters
  ControlNet outNet(ui.GetFileName("CNET"));
  FileList inList;
  if(ui.WasEntered("FROMLIST")) {
    //inList = ui.GetFileName("FROMLIST");
    inList.read(ui.GetFileName("FROMLIST"));
  }

  int inputPoints = outNet.GetNumPoints();
  int inputMeasures = 0;
  for (int cp = 0; cp < outNet.GetNumPoints(); cp++)
    inputMeasures += outNet.GetPoint(cp)->GetNumMeasures();

  // Set up the Serial Number to FileName mapping
  QMap<iString, iString> sn2filename;
  for(int cubeIndex = 0; cubeIndex < (int)inList.size(); cubeIndex ++) {
    iString sn = SerialNumber::Compose(inList[cubeIndex].toString());
    sn2filename[sn] = inList[cubeIndex].toString();
  }


  Progress progress;
  progress.SetMaximumSteps(outNet.GetNumPoints());
  progress.CheckStatus();

  // Set up vector records of how points/measures are removed
  QVector<iString> ignoredPoints;
  QVector<iString> ignoredMeasures;
  QVector<iString> singleMeasurePoints;
  QVector<iString> measurelessPoints;
  QVector<iString> tolerancePoints;
  QVector<iString> nonReferenceMeasures;
  QVector<iString> nonFixedPoints;
  QVector<iString> nonCubePoints;
  QVector<iString> noCubeMeasures;
  QVector<iString> noMeasurePoints;
  QVector<iString> nonListedPoints;
  QVector<iString> nonLatLonPoints;
  QVector<iString> cannotGenerateLatLonPoints;

  // Set up comparison data
  QVector<iString> serialNumbers;
  if(cubePoints) {
    FileList cubeList(ui.GetFileName("CUBELIST"));
    for(int cubeIndex = 0; cubeIndex < (int)cubeList.size(); cubeIndex ++) {
      iString sn = SerialNumber::Compose(cubeList[cubeIndex].toString());
      serialNumbers.push_back(sn);
    }
  }

  double tolerance = 0.0;
  if(noTolerancePoints) {
    tolerance = ui.GetDouble("PIXELTOLERANCE");
  }

  // Set up extracted network values
  if(ui.WasEntered("NETWORKID"))
    outNet.SetNetworkId(ui.GetString("NETWORKID"));

  outNet.SetUserName(Isis::Application::UserName());
  outNet.SetDescription(ui.GetString("DESCRIPTION"));

  for(int cp = outNet.GetNumPoints() - 1; cp >= 0; cp --) {
    progress.CheckStatus();

    ControlPoint *controlpt = outNet.GetPoint(cp);

    // Do preliminary exclusion checks
    if(noIgnore && controlpt->IsIgnored()) {
      ignoredPoints.append(controlpt->GetId());
      omit(outNet, cp);
      continue;
    }
    if(fixed && !(controlpt->GetType() == ControlPoint::Fixed)) {
      nonFixedPoints.append(controlpt->GetId());
      omit(outNet, cp);
      continue;
    }

    if(noSingleMeasure) {
      bool invalidPoint = false;
      invalidPoint |= noIgnore && (controlpt->GetNumValidMeasures() < 2);
      invalidPoint |= controlpt->GetNumMeasures() < 2 && (controlpt->GetType() != ControlPoint::Fixed);

      if(invalidPoint) {
        singleMeasurePoints.append(controlpt->GetId());
        omit(outNet, cp);
        continue;
      }
    }

    // Change the current point into a new point by manipulation of its control measures
    ControlPoint *newPoint = outNet.GetPoint(cp);
    bool replaceLock = false;
    if (newPoint->IsEditLocked()) {
      newPoint->SetEditLock(false);
      replaceLock = true;
    }

    for(int cm = newPoint->GetNumMeasures() - 1; cm >= 0; cm --) {
      const ControlMeasure *newMeasure = newPoint->GetMeasure(cm);

      if(noIgnore && newMeasure->IsIgnored()) {
        ignoredMeasures.append(newPoint->GetId() + "," + newMeasure->GetCubeSerialNumber());
        //New error with deleting Reference Measures
        if(newPoint->GetRefMeasure() != newMeasure)
          omit(newPoint, cm);
      }
      else if(reference && newPoint->GetRefMeasure() != newMeasure) {
        nonReferenceMeasures.append(newPoint->GetId() + "," + newMeasure->GetCubeSerialNumber());
        omit(newPoint, cm);
      }
      else if(cubeMeasures) {
        bool hasSerialNumber = false;
        std::string serialNum = newMeasure->GetCubeSerialNumber();
        for(int sn = 0; sn < serialNumbers.size(); sn ++) {
          if(serialNumbers[sn] == serialNum) {
            hasSerialNumber = true;
            break;
          }
        }

        if(!hasSerialNumber) {
          string msg = newPoint->GetId() + "," + newMeasure->GetCubeSerialNumber();
          //Delete Reference Measures not in the list, if retainReference is turned off
          if(newPoint->GetRefMeasure() != newMeasure ||
             (newPoint->GetRefMeasure() == newMeasure && !retainReference))
            omit(newPoint, cm);
          else {
            if(newPoint->GetRefMeasure() == newMeasure && retainReference) {
              msg += ", Reference not in the list but Retained";
            }
          }

          noCubeMeasures.append(msg);
        }
      }
    }

    if (replaceLock)
      newPoint->SetEditLock(true);

    //outNet.UpdatePoint(newPoint); // Fixed by redesign

    // Check for line/sample errors above provided tolerance
    if(noTolerancePoints) {
      bool hasLowTolerance = true;

      for(int cm = 0; cm < newPoint->GetNumMeasures() && hasLowTolerance; cm ++) {
        const ControlMeasure *newMeasure = newPoint->GetMeasure(cm);
        if(newMeasure->GetSampleResidual() >= tolerance ||
            newMeasure->GetLineResidual() >= tolerance) {
          hasLowTolerance = false;
        }
      }

      if(hasLowTolerance) {
        tolerancePoints.append(newPoint->GetId());
        omit(outNet, cp);
        continue;
      }
    }

    // Do not add outPoint if it has too few measures
    if(noSingleMeasure) {
      bool invalidPoint = false;
      invalidPoint |= noIgnore && (newPoint->GetNumValidMeasures() < 2);
      invalidPoint |= newPoint->GetNumMeasures() < 2 && newPoint->GetType() != ControlPoint::Fixed;

      if(invalidPoint) {
        singleMeasurePoints.append(controlpt->GetId());
        omit(outNet, cp);
        continue;
      }
    }

    // Do not add outPoint if it does not have a cube in CUBELIST as asked
    if(cubePoints) {
      bool hasSerialNumber = false;

      for(int cm = 0; cm < newPoint->GetNumMeasures() && !hasSerialNumber; cm ++) {
        for(int sn = 0; sn < serialNumbers.size() && !hasSerialNumber; sn ++) {
          if(serialNumbers[sn] == newPoint->GetMeasure(cm)->GetCubeSerialNumber())
            hasSerialNumber = true;
        }
      }

      if(!hasSerialNumber) {
        nonCubePoints.append(newPoint->GetId());
        omit(outNet, cp);
        continue;
      }
    }

    if(noMeasureless && newPoint->GetNumMeasures() == 0) {
      noMeasurePoints.append(newPoint->GetId());
      omit(outNet, cp);
      continue;
    }
  } //! Finished with simple comparisons

  /**
   * Use another pass to check for Ids
   */
  if(pointsEntered) {
    ExtractPointList(outNet, nonListedPoints);
  }


  /**
   *  Use another pass on outNet, because this is by far the most time consuming
   *  process, and time could be saved by using the reduced size of outNet
   */
  if(latLon) {
    ExtractLatLonRange(outNet, nonLatLonPoints, cannotGenerateLatLonPoints, sn2filename);
  }


  // Write the filenames associated with outNet
  WriteCubeOutList(outNet, sn2filename);

  Progress outProgress;
  outProgress.SetText("Writing Control Network");
  outProgress.SetMaximumSteps(3);
  outProgress.CheckStatus();

  // Write the extracted Control Network
  outNet.Write(ui.GetFileName("ONET"));

  outProgress.CheckStatus();

  // Adds the remove history to the summary and results group
  PvlGroup summary("ResultSummary");
  PvlGroup results("Results");

  summary.AddKeyword(PvlKeyword("InputPoints", iString(inputPoints)));
  summary.AddKeyword(PvlKeyword("InputMeasures", iString(inputMeasures)));

  int outputPoints = outNet.GetNumPoints();
  int outputMeasures = 0;
  for (int cp = 0; cp < outNet.GetNumPoints(); cp++)
    outputMeasures += outNet.GetPoint(cp)->GetNumMeasures();

  summary.AddKeyword(PvlKeyword("OutputPoints", iString(outputPoints)));
  summary.AddKeyword(PvlKeyword("OutputMeasures", iString(outputMeasures)));

  if(noIgnore) {
    summary.AddKeyword(PvlKeyword("IgnoredPoints", iString((int)ignoredPoints.size())));
    summary.AddKeyword(PvlKeyword("IgnoredMeasures", iString((int)ignoredMeasures.size())));
  }
  if(noSingleMeasure) {
    summary.AddKeyword(PvlKeyword("SingleMeasurePoints", iString((int)singleMeasurePoints.size())));
  }
  if(noMeasureless) {
    summary.AddKeyword(PvlKeyword("MeasurelessPoints", iString((int)measurelessPoints.size())));
  }
  if(noTolerancePoints) {
    summary.AddKeyword(PvlKeyword("TolerancePoints", iString((int)tolerancePoints.size())));
  }
  if(reference) {
    summary.AddKeyword(PvlKeyword("NonReferenceMeasures", iString((int)nonReferenceMeasures.size())));
  }
  if(fixed) {
    summary.AddKeyword(PvlKeyword("NonFixedPoints", iString((int)nonFixedPoints.size())));
  }
  if(cubePoints) {
    summary.AddKeyword(PvlKeyword("NonCubePoints", iString((int)nonCubePoints.size())));
  }
  if(noMeasurePoints.size() != 0) {
    summary.AddKeyword(PvlKeyword("NoCubeMeasure", iString((int)noMeasurePoints.size())));
  }
  if(cubeMeasures) {
    summary.AddKeyword(PvlKeyword("NoMeasurePoints", iString((int)noCubeMeasures.size())));
  }
  if(pointsEntered) {
    summary.AddKeyword(PvlKeyword("NonListedPoints", iString((int)nonListedPoints.size())));
  }
  if(latLon) {
    summary.AddKeyword(PvlKeyword("LatLonOutOfRange", iString((int)nonLatLonPoints.size())));
    summary.AddKeyword(PvlKeyword("NoLatLonPoints", iString((int)cannotGenerateLatLonPoints.size())));
  }

  outProgress.CheckStatus();

  // Log Control Net results
  Application::Log(summary);

  outProgress.CheckStatus();

  if(ui.WasEntered("PREFIX")) {
    Progress resultsProgress;
    resultsProgress.SetText("Writing Results");
    resultsProgress.SetMaximumSteps(11);
    resultsProgress.CheckStatus();

    std::string prefix = ui.GetString("PREFIX");

    if(noIgnore) {
      iString namecp = FileName(prefix + "IgnoredPoints.txt").expanded();
      WriteResults(namecp, ignoredPoints);
      iString namecm = FileName(prefix + "IgnoredMeasures.txt").expanded();
      WriteResults(namecm, ignoredMeasures);
    }

    resultsProgress.CheckStatus();

    if(noSingleMeasure) {
      iString name = FileName(prefix + "SingleMeasurePoints.txt").expanded();
      WriteResults(name, singleMeasurePoints);
    }

    resultsProgress.CheckStatus();

    if(noMeasureless) {
      iString name = FileName(prefix + "MeasurelessPoints.txt").expanded();
      WriteResults(name, measurelessPoints);
    }

    resultsProgress.CheckStatus();

    if(noTolerancePoints) {
      iString name = FileName(prefix + "TolerancePoints.txt").expanded();
      WriteResults(name, tolerancePoints);
    }

    resultsProgress.CheckStatus();

    if(reference) {
      iString name = FileName(prefix + "NonReferenceMeasures.txt").expanded();
      WriteResults(name, nonReferenceMeasures);
    }

    resultsProgress.CheckStatus();

    if(fixed) {
      iString name = FileName(prefix + "NonFixedPoints.txt").expanded();
      WriteResults(name, nonFixedPoints);
    }

    resultsProgress.CheckStatus();

    if(cubePoints) {
      iString name = FileName(prefix + "NonCubePoints.txt").expanded();
      WriteResults(name, nonCubePoints);
    }

    resultsProgress.CheckStatus();

    if(noMeasurePoints.size() != 0) {
      iString name = FileName(prefix + "NoMeasurePoints.txt").expanded();
      WriteResults(name, noMeasurePoints);
    }

    resultsProgress.CheckStatus();

    if(cubeMeasures) {
      iString name = FileName(prefix + "NonCubeMeasures.txt").expanded();
      WriteResults(name, noCubeMeasures);
    }

    resultsProgress.CheckStatus();

    if(pointsEntered) {
      iString name = FileName(prefix + "NonListedPoints.txt").expanded();
      WriteResults(name, nonListedPoints);
    }

    resultsProgress.CheckStatus();

    if(latLon) {
      iString namenon = FileName(prefix + "LatLonOutOfRange.txt").expanded();
      WriteResults(namenon, nonLatLonPoints);
      iString namegen = FileName(prefix + "NoLatLonPoints.txt").expanded();
      WriteResults(namegen, cannotGenerateLatLonPoints);
    }

    results.AddComment("Each keyword represents a filter parameter used." \
                       " Check the documentation for specific keyword descriptions.");
    Application::Log(results);

    resultsProgress.CheckStatus();
  }

}


/**
 * Removes control points not listed in POINTLIST
 *
 * @param outNet The output control net being removed from
 * @param nonListedPoints The keyword recording all of the control points
 *                        removed due to not being listed
 */
void ExtractPointList(ControlNet &outNet, QVector<iString> nonListedPoints) {
  UserInterface &ui = Application::GetUserInterface();

  FileList listedPoints(ui.GetFileName("POINTLIST"));

  for(int cp = outNet.GetNumPoints() - 1; cp >= 0; cp --) {
    ControlPoint *controlpt = outNet.GetPoint(cp);
    bool isInList = false;
    for(int pointId = 0; pointId < (int)listedPoints.size()  &&  !isInList; pointId ++) {
      isInList = controlpt->GetId().compare(listedPoints[pointId].toString()) == 0;
    }

    if(!isInList) {
      nonListedPoints.append(controlpt->GetId());
      omit(outNet, cp);
    }
  }
}


/**
 * Removes control points not in the lat/lon range provided in the unput
 * parameters.
 *
 * @param outNet The output control net being removed from
 * @param noLanLonPoint The keyword recording all of the control points removed
 *                      due to the provided lat/lon range
 * @param noLanLonPoint The keyword recording all of the control points removed
 *                      due to the inability to calculate the lat/lon for that
 *                      point
 */
void ExtractLatLonRange(ControlNet &outNet, QVector<iString> nonLatLonPoints,
                        QVector<iString> cannotGenerateLatLonPoints,  QMap<iString, iString> sn2filename) {
  if(outNet.GetNumPoints() == 0) {
    return;
  }

  UserInterface &ui = Application::GetUserInterface();

  // Get the lat/lon and fix the range for the internal 0/360
  Latitude minlat(ui.GetDouble("MINLAT"), Angle::Degrees);
  Latitude maxlat(ui.GetDouble("MAXLAT"), Angle::Degrees);
  Longitude minlon = Longitude(ui.GetDouble("MINLON"), Angle::Degrees);
  Longitude maxlon = Longitude(ui.GetDouble("MAXLON"), Angle::Degrees);

  Progress progress;
  progress.SetText("Calculating lat/lon");
  progress.SetMaximumSteps(outNet.GetNumPoints());
  progress.CheckStatus();

  CubeManager manager;
  manager.SetNumOpenCubes(50);   //Should keep memory usage to around 1GB

  bool hasFromList = ui.WasEntered("FROMLIST");
  for(int cp = outNet.GetNumPoints() - 1; cp >= 0; cp --) {
    progress.CheckStatus();
    const ControlPoint *controlPt = outNet.GetPoint(cp);
    SurfacePoint surfacePt = controlPt->GetBestSurfacePoint();

    // If the Contorl Network takes priority, use it
    if(surfacePt.Valid()) {
      if(NotInLatLonRange(surfacePt, minlat, maxlat, minlon, maxlon)) {
        nonLatLonPoints.push_back(controlPt->GetId());
        omit(outNet, cp);
      }
    }

    /**
     * If the lat/lon cannot be determined from the point, then we need to calculate
     * lat/lon on our own
     */
    else if(hasFromList) {

      // Find a cube in the Control Point to get the lat/lon from
      int cm = 0;
      iString sn = "";
      Latitude lat;
      Longitude lon;
      Distance radius;

      // First check the reference Measure
      //if(!sn2filename[controlPt[cm].GetCubeSerialNumber()].length() == 0) {
      if(!sn2filename[controlPt->GetReferenceSN()].length() == 0) {
        sn = controlPt->GetReferenceSN();
      }

      // Search for other Control Measures if needed
      if(sn.empty()) {
        // Find the Serial Number if it exists
        for(int cm = 0; (cm < controlPt->GetNumMeasures()) && sn.empty(); cm ++) {
          if(!sn2filename[controlPt->GetReferenceSN()].length() == 0) {
            sn = controlPt->GetReferenceSN();
          }
        }
      }

      // Connot fine a cube to get the lat/lon from
      if(sn.empty()) {
        cannotGenerateLatLonPoints.push_back(controlPt->GetId());
        omit(outNet, cp);
      }

      // Calculate the lat/lon and check for validity
      else {
        bool remove = false;

        Cube *cube = manager.OpenCube(sn2filename[sn]);
        Camera *camera = cube->getCamera();

        if(camera == NULL) {
          try {
            Projection *projection =
              ProjectionFactory::Create((*(cube->getLabel())));

            if(!projection->SetCoordinate(controlPt->GetMeasure(cm)->GetSample(),
                                          controlPt->GetMeasure(cm)->GetLine())) {
              nonLatLonPoints.push_back(controlPt->GetId());
              remove = true;
            }

            lat = Latitude(projection->Latitude(), Angle::Degrees);
            lon = Longitude(projection->Longitude(), Angle::Degrees);
            radius = Distance(projection->LocalRadius(), Distance::Meters);

            delete projection;
            projection = NULL;
          }
          catch(IException &) {
            remove = true;
          }
        }
        else {
          if(!camera->SetImage(controlPt->GetMeasure(cm)->GetSample(),
                               controlPt->GetMeasure(cm)->GetLine())) {
            nonLatLonPoints.push_back(controlPt->GetId());
            remove = true;
          }

          lat = camera->GetLatitude();
          lon = camera->GetLongitude();
          radius = camera->LocalRadius();

          camera = NULL;
        }

        cube = NULL;

        bool notInRange = false;
        bool validLatLonRadius = lat.isValid() && lon.isValid() && radius.isValid();
        if(validLatLonRadius) {
          SurfacePoint sfpt(lat, lon, radius);
          notInRange = NotInLatLonRange(sfpt, minlat, maxlat, minlon, maxlon);
        }

        if(remove || notInRange) {
          nonLatLonPoints.push_back(controlPt->GetId());
          omit(outNet, cp);
        }
        else if(validLatLonRadius) { // Add the reference lat/lon/radius to the Control Point
          outNet.GetPoint(cp)->SetAprioriSurfacePoint(SurfacePoint(lat, lon, radius));
        }
      }
    }
    else {
      cannotGenerateLatLonPoints.push_back(controlPt->GetId());
      omit(outNet, cp);
    }

  }

  manager.CleanCubes();
}


/**
 * Checks for correct lat/lon range, handling the meridian correctly
 *
 * @param lat The latitude to check
 * @param lon The longitude to check
 * @param minlat Minimum Latitude Minimum valid latitude
 * @param maxlat Maximum Latitude Maximum valid latitude
 * @param minlon Minimum Longitude Minimum valid longitude
 * @param maxlon Maximum Longitude Maximum valid longitude
 *
 * @return bool True when the range is valid
 */
bool NotInLatLonRange(SurfacePoint surfacePtToTest, Latitude minlat,
                      Latitude maxlat, Longitude minlon, Longitude maxlon) {
  Latitude lat = surfacePtToTest.GetLatitude();
  Longitude lon = surfacePtToTest.GetLongitude();

  bool outRange = false;
  try {
    outRange = !lat.inRange(minlat, maxlat) || !lon.inRange(minlon, maxlon);
  }
  catch (IException &e) {
    iString msg = "Cannot complete lat/lon range test with given filters";
    throw IException(e, IException::User, msg, _FILEINFO_);
  }

  return outRange;
}


/**
 * Finds and writes all input cubes contained within the given Control Network
 * to the output file list
 *
 * @param cnet The Control Network to list the filenames contained within
 * @param sn2file The map for converting the Control Network's serial numbers
 *                to filenames
 */
void WriteCubeOutList(ControlNet cnet, QMap<iString, iString> sn2file) {
  UserInterface &ui = Application::GetUserInterface();

  if(ui.WasEntered("TOLIST")) {

    Progress p;
    p.SetText("Writing Cube List");
    try {
      p.SetMaximumSteps(cnet.GetNumPoints());
      p.CheckStatus();
    }
    catch(IException &e) {
      std::string msg = "The provided filters have resulted in an empty Control Network.";
      throw IException(e, IException::User, msg, _FILEINFO_);
    }

    std::set<iString> outputsn;
    for(int cp = 0; cp < cnet.GetNumPoints(); cp ++) {
      for(int cm = 0; cm < cnet.GetPoint(cp)->GetNumMeasures(); cm ++) {
        outputsn.insert(cnet.GetPoint(cp)->GetMeasure(cm)->GetCubeSerialNumber());
      }
      p.CheckStatus();
    }

    std::string toList = ui.GetFileName("TOLIST");
    std::ofstream out_stream;
    out_stream.open(toList.c_str(), std::ios::out);
    out_stream.seekp(0, std::ios::beg);   //Start writing from beginning of file

    for(std::set<iString>::iterator sn = outputsn.begin(); sn != outputsn.end(); sn ++) {
      if(!sn2file[(*sn)].length() == 0) {
        out_stream << sn2file[(*sn)] << std::endl;
      }
    }

    out_stream.close();
  }
}


/**
 * Places the output
 *
 * @param filename The file to write the vector of results to
 * @param results  A vector of points and/or measures not extracted
 */
void WriteResults(iString filename, QVector<iString> results) {
  if(results.size() == 0) {
    return;
  }

  // Set up the output file for writing
  std::ofstream out_stream;
  out_stream.open(filename.c_str(), std::ios::out);
  out_stream.seekp(0, std::ios::beg);   //Start writing from beginning of file

  out_stream << results[0];
  for(int index = 1; index < results.size(); index ++) {
    out_stream << std::endl << results[index];
  }

  out_stream.close();
}


void omit(ControlNet &cnet, int cp) {
  ControlPoint *point = cnet.GetPoint(cp);
  if (point->IsEditLocked()) point->SetEditLock(false);
  cnet.DeletePoint(cp);
}


void omit(ControlPoint *point, int cm) {
  ControlMeasure *measure = point->GetMeasure(cm);
  if (measure->IsEditLocked()) measure->SetEditLock(false);
  point->Delete(cm);
}

