#include "Isis.h"

#include <QList>
#include <QString>

#include "Buffer.h"
#include "Camera.h"
#include "CubeAttribute.h"
#include "IException.h"
#include "ProcessByLine.h"
#include "Pvl.h"
#include "PvlGroup.h"
#include "PvlKeyword.h"
#include "PvlObject.h"
#include "SpecialPixel.h"
#include "UserInterface.h"
#include "iString.h"

using namespace std;
using namespace Isis;

void calibration(vector<Buffer *> &in,
                 vector<Buffer *> &out);
PvlObject fetchCoefficients(Pvl &calibration, QList<QString> &hierarchy);
void checkCoefficient(PvlObject &coefficients, iString keyName);

// Linearity correction
bool linear;

// calibration constants
double omegaNaught;
double sunDist;
double gainC;
double offsetC;
double deltaExpo;

// calculated values
double aCoef;
double newExpo;
double XMLT;
double w1;
double dist1;

// linear constants
double bHighEnd;
double kPowerOf;
double normalizingPower;

void IsisMain() {
  // Will be processing by line
  ProcessByLine p;

  UserInterface &ui = Application::GetUserInterface();

  Cube * incube = p.SetInputCube("FROM");
  Cube * ocube = p.SetOutputCube("TO");

  // Check for projection
  if (incube->isProjected()) {
    string msg = "The cube [" + ui.GetFileName("FROM") + "] has a projection" +
                 " and cannot be radiometrically calibrated";
    throw IException(IException::User, msg, _FILEINFO_);
  }

  // access important label objects, will be used later.
  PvlObject isiscube = incube->getLabel()->FindObject("IsisCube");
  PvlGroup instrument = isiscube.FindGroup("Instrument");
  PvlGroup archive = isiscube.FindGroup("Archive");
  PvlGroup bandbin = isiscube.FindGroup("BandBin");

  // Verify not radiometrically corrected
  if (isiscube.HasGroup("Radiometry")) {
    string msg = "Cube [" + ui.GetFileName("FROM") + "] has already been" +
                 " radiometrically corrected";
    throw IException(IException::User,msg,_FILEINFO_);
  }

  // Verify Voyager spacecraft and get number, 1 or 2
  string scNumber = instrument["SpacecraftName"][0];
  if (scNumber != "VOYAGER_1" && scNumber != "VOYAGER_2") {
    string msg = "The cube [" + ui.GetFileName("FROM") + "] does not appear" +
                 " to be a Voyager image";
    throw IException(IException::User, msg, _FILEINFO_);
  }
  scNumber = scNumber[8];

  // Open calibration file to find constants and files
  Pvl calibra(FileName("$voyager" + scNumber +
                           "/calibration/voycal.pvl").expanded());
  PvlObject calib;
  QList<QString> hierarchy;
  try {
    // Search voycal.pvl for appropriate object
    hierarchy.append(instrument.FindKeyword("SpacecraftName")[0]);
    hierarchy.append(instrument.FindKeyword("InstrumentId")[0]);
    hierarchy.append(
        iString("ShutterMode" + instrument.FindKeyword("CameraState2")[0]));
    hierarchy.append(archive.FindKeyword("MissionPhaseName")[0]);
    hierarchy.append(
        iString("ScanRate" + instrument.FindKeyword("CameraState1")[0]));
    hierarchy.append(
        iString(bandbin["FilterName"][0] + "_" + bandbin["FilterNumber"][0]));

    calib = fetchCoefficients(calibra, hierarchy);

    checkCoefficient(calib, "OmegaNaught");
    checkCoefficient(calib, "SunDistance");
    checkCoefficient(calib, "GainCorrection");
    checkCoefficient(calib, "OffsetCorrection");
    checkCoefficient(calib, "DeltaExposureTime");
  }
  catch (IException &e) {
    string msg = "Could not find match in [voycal.pvl] calibration file,";
    msg += " the error was: ";
    msg += e.what();
    throw IException(IException::Unknown, msg, _FILEINFO_);
  }

  // Get appropriate calibration files
  CubeAttributeInput in1;
  p.SetInputCube(FileName("$voyager" + scNumber + "/calibration/" +
                          (string)calib["OffsetCorrectionFile"]).expanded(), in1);
  CubeAttributeInput in2;
  p.SetInputCube(FileName("$voyager" + scNumber + "/calibration/" +
                          (string)calib["GainCorrectionFile"]).expanded(), in2);

  // Constants from voycal.pvl for correction
  omegaNaught = calib["OmegaNaught"];
  sunDist = calib["SunDistance"];
  gainC = calib["GainCorrection"];
  offsetC = calib["OffsetCorrection"];
  deltaExpo = calib["DeltaExposureTime"];

  // If we are doing a linear correction as well, go here
  linear = ui.GetBoolean("LINEAR");
  if (linear) {
    Pvl linearity(FileName("$voyager" + scNumber +
                     "/calibration/voylin.pvl").expanded());

    PvlObject lin;
    try {
      // Search voylin.pvl for appropriate object
      lin = fetchCoefficients(linearity, hierarchy);

      checkCoefficient(lin, "NormalizingPower");
      checkCoefficient(lin, "B_HighEndNon-LinearityCorrection");
      checkCoefficient(lin, "K_PowerOfNon-Linearity");
    }
    catch (IException &e) {
      string msg = "Could not find match in [voylin.pvl] calibration file,";
      msg += " the error was: ";
      msg += e.what();
      throw IException(IException::Unknown, msg, _FILEINFO_);
    }
    // Calculated constant for linear correction
    aCoef = ((double)lin["NormalizingPower"] - (double)lin["B_HighEndNon-LinearityCorrection"]) /
            (double)lin["NormalizingPower"];

    // Constants from file for linear correction
    bHighEnd = (double)lin["B_HighEndNon-LinearityCorrection"];
    kPowerOf = (double)lin["K_PowerOfNon-Linearity"];
    normalizingPower = (double)lin["NormalizingPower"];
  }

  // Other calculated constants for all corrections
  // Modified exposure duration
  newExpo = ((double)instrument["ExposureDuration"] * 1000) + (double)calib["DeltaExposureTime"];

  // Camera for sun position/distance, so we have to have spice data
  Camera * cam = incube->getCamera();
  double sunPos[3];
  cam->sunPosition(sunPos);
  // Magnitude of sun vector ||s||
  dist1 = sqrt(sunPos[0]*sunPos[0] + sunPos[1]*sunPos[1] + sunPos[2]*sunPos[2]);

  w1 = (omegaNaught/1000) * ((sunDist * sunDist) / (dist1 * dist1));

  XMLT = 1.0/(newExpo * w1);

  PvlGroup calgrp("Radiometry");
  // Regular calibration equation and constants
  calgrp.AddComment("Calibration equation in voycal:");
  calgrp.AddComment("OUT(i,j) = (IN(i,j)*GAIN)+DCF(i,j)*XMLT*FFF(i,j)");
  calgrp.AddComment("XMLT = 1.0/(EXPO*W1)");
  calgrp.AddComment("EXPO = EXPODUR + DELTAEXPO");
  calgrp.AddComment("W1 = (W0/1000)*(SUNDIST^2/CALCDIST^2)");
  calgrp.AddComment("DCF = OffsetCorrectionFile, FFF = GainCorrectionFile");
  calgrp.AddComment("IN = InputCube, GAIN = GainCorrection");
  calgrp.AddKeyword(calib[0]);
  calgrp.AddKeyword(calib[1]);
  calgrp.AddKeyword(calib[2]);
  calgrp.AddKeyword(calib[3]);
  calgrp.AddKeyword(calib[4]);
  calgrp.AddKeyword(calib[5]);
  calgrp.AddKeyword(calib[6]);
  calgrp.AddKeyword(PvlKeyword("CalcSunDistance",dist1));
  calgrp.AddKeyword(instrument["ExposureDuration"]);
  calgrp.AddKeyword(PvlKeyword("XMLT",XMLT));
  calgrp.AddKeyword(PvlKeyword("Omega_W1",w1));
  calgrp.AddKeyword(PvlKeyword("CalcExpoDuration",newExpo));

  // Linear correction equation and constants
  if (linear) {
    PvlKeyword linearity = PvlKeyword("LinearityCorrection","True");
    linearity.AddComment("Linearity correction equation:");
    linearity.AddComment("OUT(i,j) = LIN( (IN(i,j)*GAIN)+DCF(i,j) )*XMLT*FFF(i,j)");
    linearity.AddComment("LIN(X) = ACOEF*X+BCOEF*(X/XNORM)^KPOWER");
    linearity.AddComment("BCOEF = B_HighEndeNon-LinearityCorrection");
    linearity.AddComment("XNORM = NormalizingPower");
    linearity.AddComment("KPOWER = K_PowerOfNon-Linearity");
    calgrp.AddKeyword(linearity);
    calgrp.AddKeyword(PvlKeyword("ACoefficient",aCoef));
    calgrp.AddKeyword(PvlKeyword("B_HighEndNon-LinearityCorrection",bHighEnd));
    calgrp.AddKeyword(PvlKeyword("K_PowerOfNon-Linearity",kPowerOf));
    calgrp.AddKeyword(PvlKeyword("NormalizingPower",normalizingPower));
  }
  else {
    calgrp.AddKeyword(PvlKeyword("LinearityCorrection","False"));
  }

  // Add Radiometry group
  ocube->putGroup(calgrp);

  p.StartProcess(calibration);
  p.EndProcess();
}

// calibration method
void calibration(vector<Buffer *> &in, vector<Buffer *> &out) {

  Buffer &inp = *in[0];      // Input Cube
  Buffer &dcf = *in[1];      // Dark Current File
  Buffer &fff = *in[2];      // Flat Field File
  Buffer &outp = *out[0];    // Output Cube

  for (int i = 0; i < inp.size(); i++) {
    if(IsSpecial(inp[i])) {
      outp[i] = inp[i];
    }
    else if(IsSpecial(fff[i]) || IsSpecial(dcf[i])) {
      outp[i] = Isis::Null;
    }
    else {
      // Initial calculation
      double dnraw = (inp[i] * gainC) + dcf[i];
      // Linear option
      if (linear) {
        dnraw = aCoef * dnraw + bHighEnd * (pow((dnraw / normalizingPower), kPowerOf));
      }
      outp[i] = XMLT * fff[i] * dnraw;
    }
  }
}


/**
 * Fetches all the coefficients for a calibration in a single top-level
 * PvlObject.  The calibration PVL will be explored from the top-down looking
 * for keywords.  Objects named in the hierarchy list will be explored down
 * the chain in the order listed.  So, for example, if hierarchy is the list
 * ["A", "B", "C"], then this function will assume that "A" is a child of the
 * top-level calibration PvlObject, "B" is a child of "A", and "C" is a child
 * of "B".
 *
 * If any keyword found in a child object further down the chain conflicts
 * with a keyword found higher up, the child object's keyword will overwrite
 * the value of that keyword in the coefficients object.  In this way,
 * calibration files can define top-level "default" coefficients that apply to
 * the entire mission, and can provide camera- or filter-specific coefficients
 * to overwrite them as needed.
 *
 * @param calibration The mission-specific calibration PVL file
 * @param hierarchy A list of hierarchical PvlObject names to search
 *
 * @return A single object containing all coefficients
 */
PvlObject fetchCoefficients(Pvl &calibration, QList<QString> &hierarchy) {
  // All coefficients go into one top-level object without any children
  PvlObject coefficients;

  // Add all the keywords from the calibration PVL top-level object
  for (int k = 0; k < calibration.Keywords(); k++)
    coefficients.AddKeyword(calibration[k]);

  // Iterate over every object in the hierarchy looking for coefficient
  // keywords.  The first string is the name of the first object, the second
  // the name of the next object down, etc. We start from the top-level PVL
  // object and work our way down.
  bool validHierarchy = true;
  PvlObject *parent = &calibration;
  for (int o = 0; o < hierarchy.size() && validHierarchy; o++) {
    iString objectName = hierarchy[o].toStdString();

    if (parent->HasObject(objectName)) {
      // The object named in the hierarchy exists in the calibration file, so
      // grab it
      PvlObject &object = parent->FindObject(objectName);

      // Find all the keywords at the object level
      for (int k = 0; k < object.Keywords(); k++) {
        PvlKeyword &keyword = object[k];
        if (coefficients.HasKeyword(keyword.Name())) {
          // The coefficients object already has a value for this coefficient
          // keyword.  Because this one is lower down the chain, it is more
          // specifically defined, thus we should use it instead.
          PvlKeyword &coefficient = coefficients.FindKeyword(keyword.Name());
          for (int i = 0; i < coefficient.Size(); i++) {
            coefficient[i] = keyword[i];
          }
        }
        else {
          // This is a new keyword, so add it to the coefficients object
          coefficients.AddKeyword(object[k]);
        }
      }

      // Update the parent object for the next object down the chain
      parent = &object;
    }
    else {
      // We've reached a dead end, our hierarchy no longer makes sense, so
      // bail out
      validHierarchy = false;
    }
  }

  return coefficients;
}


void checkCoefficient(PvlObject &coefficients, iString keyName) {
  if (!coefficients.HasKeyword(keyName))
    throw IException(
        IException::Programmer,
        "Coefficient [" + keyName + "] was not found in the calibration PVL "
        "for the input data.  Consider adding a default.",
        _FILEINFO_);
}

