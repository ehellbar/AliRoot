/*************************************************************************
* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
*                                                                        *
* Author: The ALICE Off-line Project.                                    *
* Contributors are mentioned in the code where appropriate.              *
*                                                                        *
* Permission to use, copy, modify and distribute this software and its   *
* documentation strictly for non-commercial purposes is hereby granted   *
* without fee, provided that the above copyright notice appears in all   *
* copies and that both the copyright notice and this permission notice   *
* appear in the supporting documentation. The authors make no claims     *
* about the suitability of this software for any purpose. It is          *
* provided "as is" without express or implied warranty.                  *
**************************************************************************/


/* $Id$ */

/// \class AliTPCSpaceCharge3DDriftLine
/// \brief This class provides an interface for the usage of the distortion and correction maps calculated by AliTPCSpaceCharge3DDriftLine
///
/// \author Rifki Sadikin <rifki.sadikin@cern.ch>, Indonesian Institute of Sciences
/// \date Jan 31, 2018

#include "TGeoGlobalMagField.h"
#include "TFile.h"
#include "TTreeStream.h"
#include "TMath.h"
#include "TH2F.h"

#include "AliTPCSpaceCharge3DDriftLine.h"
#include "AliMagF.h"
#include "AliTPCParam.h"
#include "AliTPCParamSR.h"
#include "AliTPCcalibDB.h"
#include "AliLog.h"

/// \cond CLASSIMP
ClassImp(AliTPCSpaceCharge3DDriftLine)
/// \endcond

/// Construction for AliTPCSpaceCharge3DDriftLine class
/// Default values
AliTPCSpaceCharge3DDriftLine::AliTPCSpaceCharge3DDriftLine()
  : AliTPCCorrection(), AliTPCSpaceCharge3DCalc()
{
};

/// Construction for AliTPCSpaceCharge3DDriftLine class
/// Default values
AliTPCSpaceCharge3DDriftLine::AliTPCSpaceCharge3DDriftLine(const char *name, const char *title)
  : AliTPCCorrection(name, title), AliTPCSpaceCharge3DCalc()
{
};

/// Construction for AliTPCSpaceCharge3DDriftLine class
/// Member values from params
///
/// \param nRRow Int_t number of grid in r direction
/// \param nZColumn Int_t number of grid in z direction
/// \param nPhiSlice Int_t number of grid in \f$ \phi \f$ direction
///
AliTPCSpaceCharge3DDriftLine::AliTPCSpaceCharge3DDriftLine(const char *name, const char *title, Int_t nRRow, Int_t nZColumn, Int_t nPhiSlice)
  : AliTPCCorrection(name, title), AliTPCSpaceCharge3DCalc(nRRow, nZColumn, nPhiSlice)
{
}

AliTPCSpaceCharge3DDriftLine::AliTPCSpaceCharge3DDriftLine(const char *name, const char *title, Int_t nRRow, Int_t nZColumn, Int_t nPhiSlice, Int_t interpolationOrder, Int_t irregularGridSize, Int_t rbfKernelType)
  : AliTPCCorrection(name, title), AliTPCSpaceCharge3DCalc(nRRow, nZColumn, nPhiSlice, interpolationOrder, irregularGridSize, rbfKernelType)
{
}

/// Init copy from AliTPCSpaceCharge3D
void AliTPCSpaceCharge3DDriftLine::Init()
{
  AliMagF *magF = (AliMagF *) TGeoGlobalMagField::Instance()->GetField();
  if (!magF) AliError("Magnetic field - not initialized");
  Double_t bzField = magF->SolenoidField() / 10.; //field in T
  AliTPCParam *param = AliTPCcalibDB::Instance()->GetParameters();
  if (!param) AliError("Parameters - not initialized");
  Double_t vDrift = param->GetDriftV() / 1000000.; // [cm/us]   // From dataBase: to be updated: per second (ideally)
  Double_t ezField = 400; // [V/cm]   // to be updated: never (hopefully)
  Double_t wt = -10.0 * (bzField * 10) * vDrift / ezField;
  // Correction Terms for effective omegaTau; obtained by a laser calibration run
  SetOmegaTauT1T2(wt, fT1, fT2);
}

///
/// \param z
/// \param nx
/// \param ny
/// \return
TH2F *AliTPCSpaceCharge3DDriftLine::CreateHistogramDistDRInXY(Float_t z, Int_t nx, Int_t ny)
{
  /// Simple plot functionality.
  /// Returns a 2d histogram which represents the corrections in radial direction (drDist)
  /// in respect to position z within the XY plane.
  /// The histogram nx times ny entries.

  AliTPCParam *tpcParam = new AliTPCParamSR;

  TH2F *h = CreateTH2F("dr_xy", TString::Format("%s: DRInXY Z=%2.0f", GetTitle(), z).Data(), "x [cm]", "y [cm]",
                       "drDist [cm]",
                       nx, -250., 250., ny, -250., 250.);
  Float_t x[3], dx[3], xCyl[3];
  Float_t r0, phi0;

  x[2] = z;
  Int_t roc = z > 0. ? 0 : 18; // FIXME
  for (Int_t iy = 1; iy <= ny; ++iy) {
    x[1] = h->GetYaxis()->GetBinCenter(iy);
    for (Int_t ix = 1; ix <= nx; ++ix) {
      x[0] = h->GetXaxis()->GetBinCenter(ix);
      r0 = TMath::Sqrt((x[0]) * (x[0]) + (x[1]) * (x[1]));
      phi0 = TMath::ATan2(x[1], x[0]);

      while (phi0 > TMath::Pi()) phi0 -= TMath::TwoPi();
      while (phi0 < -TMath::Pi()) phi0 += TMath::TwoPi();
      xCyl[0] = r0;
      xCyl[1] = phi0;
      GetDistortionCylAC(xCyl, roc, dx);

      //Float_t r0=TMath::Sqrt((x[0]      )*(x[0]      )+(x[1]      )*(x[1]      ));

      if (tpcParam->GetPadRowRadii(0, 0) <= r0 && r0 <= tpcParam->GetPadRowRadii(36, 95)) {
        //Float_t r1=TMath::Sqrt((x[0]+dx[0])*(x[0]+dx[0])+(x[1]+dx[1])*(x[1]+dx[1]));
        h->SetBinContent(ix, iy, dx[0]);
      } else
        h->SetBinContent(ix, iy, 0.);
    }
  }
  delete tpcParam;
  return h;
}

/// Simple plot functionality.
/// Returns a 2d histogram which represents the corrections in radial direction (drDist)
/// in respect to position z within the XY plane.
/// The histogram nx times ny entries.
///
/// \param z
/// \param nx
/// \param ny
/// \return
TH2F *AliTPCSpaceCharge3DDriftLine::CreateHistogramCorrDRInXY
  (
    Float_t z,
    Int_t nx,
    Int_t ny
  )
{

  AliTPCParam *tpcParam = new AliTPCParamSR;

  TH2F *h = CreateTH2F("dr_xy", TString::Format("%s: DRInXY Z=%2.0f", GetTitle(), z).Data(), "x [cm]", "y [cm]",
                       "drDist [cm]",
                       nx, -250., 250., ny, -250., 250.);
  Float_t x[3], dx[3], xCyl[3];
  Float_t r0, phi0;

  x[2] = z;
  Int_t roc = z > 0. ? 0 : 18; // FIXME
  for (Int_t iy = 1; iy <= ny; ++iy) {
    x[1] = h->GetYaxis()->GetBinCenter(iy);
    for (Int_t ix = 1; ix <= nx; ++ix) {
      x[0] = h->GetXaxis()->GetBinCenter(ix);
      r0 = TMath::Sqrt((x[0]) * (x[0]) + (x[1]) * (x[1]));
      phi0 = TMath::ATan2(x[1], x[0]);

      while (phi0 > TMath::Pi()) phi0 -= TMath::TwoPi();
      while (phi0 < -TMath::Pi()) phi0 += TMath::TwoPi();
      xCyl[0] = r0;
      xCyl[1] = phi0;
      GetCorrectionCylAC(xCyl, roc, dx);

      //Float_t r0=TMath::Sqrt((x[0]      )*(x[0]      )+(x[1]      )*(x[1]      ));

      if (tpcParam->GetPadRowRadii(0, 0) <= r0 && r0 <= tpcParam->GetPadRowRadii(36, 95)) {
        //Float_t r1=TMath::Sqrt((x[0]+dx[0])*(x[0]+dx[0])+(x[1]+dx[1])*(x[1]+dx[1]));
        h->SetBinContent(ix, iy, dx[0]);
      } else
        h->SetBinContent(ix, iy, 0.);
    }
  }
  delete tpcParam;
  return h;
}

/// Simple plot functionality.
/// Returns a 2d histogram which represents the corrections in r Phi direction (drPhi)
/// in respect to position z within the XY plane.
/// The histogram nx times ny entries.
///
/// \param z
/// \param nx
/// \param ny
/// \return
TH2F *AliTPCSpaceCharge3DDriftLine::CreateHistogramDistDRPhiInXY(Float_t z, Int_t nx, Int_t ny)
{

  AliTPCParam *tpcParam = new AliTPCParamSR;

  TH2F *h = CreateTH2F("drPhi_xy", TString::Format("%s: DRPhiInXY Z=%2.0f", GetTitle(), z).Data(), "x [cm]", "y [cm]",
                       "drPhi [cm]",
                       nx, -250., 250., ny, -250., 250.);
  Float_t x[3], dx[3];
  x[2] = z;
  Int_t roc = z > 0. ? 0 : 18; // FIXME
  for (Int_t iy = 1; iy <= ny; ++iy) {
    x[1] = h->GetYaxis()->GetBinCenter(iy);
    for (Int_t ix = 1; ix <= nx; ++ix) {
      x[0] = h->GetXaxis()->GetBinCenter(ix);
      GetDistortion(x, roc, dx);
      Float_t r0 = TMath::Sqrt((x[0]) * (x[0]) + (x[1]) * (x[1]));
      if (tpcParam->GetPadRowRadii(0, 0) <= r0 && r0 <= tpcParam->GetPadRowRadii(36, 95)) {
        Float_t phi0 = TMath::ATan2(x[1], x[0]);
        Float_t phi1 = TMath::ATan2(x[1] + dx[1], x[0] + dx[0]);

        Float_t dPhi = phi1 - phi0;
        if (dPhi < TMath::Pi()) dPhi += TMath::TwoPi();
        if (dPhi > TMath::Pi()) dPhi -= TMath::TwoPi();

        h->SetBinContent(ix, iy, r0 * dPhi);
      } else
        h->SetBinContent(ix, iy, 0.);
    }
  }

  delete tpcParam;
  return h;
}

/// Simple plot functionality.
/// Returns a 2d histogram which represents the corrections in r phi direction (drPhi)
/// in respect to position z within the XY plane.
/// The histogram nx times ny entries.
/// \param z
/// \param nx
/// \param ny
/// \return
TH2F *AliTPCSpaceCharge3DDriftLine::CreateHistogramCorrDRPhiInXY(Float_t z, Int_t nx, Int_t ny)
{

  AliTPCParam *tpcParam = new AliTPCParamSR;

  TH2F *h = CreateTH2F("drPhi_xy", TString::Format("%s: DRPhiInXY Z=%2.0f", GetTitle(), z).Data(), "x [cm]", "y [cm]",
                       "drPhi [cm]",
                       nx, -250., 250., ny, -250., 250.);
  Float_t x[3], dx[3];
  x[2] = z;
  Int_t roc = z > 0. ? 0 : 18; // FIXME
  for (Int_t iy = 1; iy <= ny; ++iy) {
    x[1] = h->GetYaxis()->GetBinCenter(iy);
    for (Int_t ix = 1; ix <= nx; ++ix) {
      x[0] = h->GetXaxis()->GetBinCenter(ix);
      GetCorrection(x, roc, dx);
      Float_t r0 = TMath::Sqrt((x[0]) * (x[0]) + (x[1]) * (x[1]));
      if (tpcParam->GetPadRowRadii(0, 0) <= r0 && r0 <= tpcParam->GetPadRowRadii(36, 95)) {
        Float_t phi0 = TMath::ATan2(x[1], x[0]);
        Float_t phi1 = TMath::ATan2(x[1] + dx[1], x[0] + dx[0]);

        Float_t dPhi = phi1 - phi0;
        if (dPhi < TMath::Pi()) dPhi += TMath::TwoPi();
        if (dPhi > TMath::Pi()) dPhi -= TMath::TwoPi();

        h->SetBinContent(ix, iy, r0 * dPhi);
      } else
        h->SetBinContent(ix, iy, 0.);
    }
  }

  delete tpcParam;
  return h;
}

/// Simple plot functionality.
/// Returns a 2d histogram which represents the corrections in longitudinal direction (dzDist)
/// in respect to position z within the XY plane.
/// The histogram nx times ny entries.
/// \param z
/// \param nx
/// \param ny
/// \return
TH2F *AliTPCSpaceCharge3DDriftLine::CreateHistogramDistDZInXY(Float_t z, Int_t nx, Int_t ny)
{

  AliTPCParam *tpcParam = new AliTPCParamSR;

  TH2F *h = CreateTH2F("dz_xy", TString::Format("%s: DZInXY Z=%2.0f", GetTitle(), z).Data(), "x [cm]", "y [cm]",
                       "dzDist [cm]",
                       nx, -250., 250., ny, -250., 250.);
  Float_t x[3], dx[3];
  x[2] = z;
  Int_t roc = z > 0. ? 0 : 18; // FIXME
  for (Int_t iy = 1; iy <= ny; ++iy) {
    x[1] = h->GetYaxis()->GetBinCenter(iy);
    for (Int_t ix = 1; ix <= nx; ++ix) {
      x[0] = h->GetXaxis()->GetBinCenter(ix);
      GetDistortion(x, roc, dx);
      Float_t r0 = TMath::Sqrt((x[0]) * (x[0]) + (x[1]) * (x[1]));
      if (tpcParam->GetPadRowRadii(0, 0) <= r0 && r0 <= tpcParam->GetPadRowRadii(36, 95)) {
        h->SetBinContent(ix, iy, dx[2]);
      } else
        h->SetBinContent(ix, iy, 0.);
    }
  }
  delete tpcParam;
  return h;
}

///
/// \param z
/// \param nx
/// \param ny
/// \return
TH2F *AliTPCSpaceCharge3DDriftLine::CreateHistogramCorrDZInXY(Float_t z, Int_t nx, Int_t ny)
{
  /// Simple plot functionality.
  /// Returns a 2d histogram which represents the corrections in longitudinal direction (dzDist)
  /// in respect to position z within the XY plane.
  /// The histogram nx times ny entries.

  AliTPCParam *tpcParam = new AliTPCParamSR;

  TH2F *h = CreateTH2F("dz_xy", TString::Format("%s: DZInXY Z=%2.0f", GetTitle(), z).Data(), "x [cm]", "y [cm]",
                       "dzDist [cm]",
                       nx, -250., 250., ny, -250., 250.);
  Float_t x[3], dx[3];
  x[2] = z;
  Int_t roc = z > 0. ? 0 : 18; // FIXME
  for (Int_t iy = 1; iy <= ny; ++iy) {
    x[1] = h->GetYaxis()->GetBinCenter(iy);
    for (Int_t ix = 1; ix <= nx; ++ix) {
      x[0] = h->GetXaxis()->GetBinCenter(ix);
      GetCorrection(x, roc, dx);
      Float_t r0 = TMath::Sqrt((x[0]) * (x[0]) + (x[1]) * (x[1]));
      if (tpcParam->GetPadRowRadii(0, 0) <= r0 && r0 <= tpcParam->GetPadRowRadii(36, 95)) {
        h->SetBinContent(ix, iy, dx[2]);
      } else
        h->SetBinContent(ix, iy, 0.);
    }
  }
  delete tpcParam;
  return h;
}

///
/// \param z
/// \param nx
/// \param ny
/// \return
TH2F *AliTPCSpaceCharge3DDriftLine::CreateHistogramSCInXY(Float_t z, Int_t nx, Int_t ny) {
  TH2F *h = CreateTH2F("spaceCharge", GetTitle(), "x [cm]", "y [cm]", "#rho_{sc} [C/m^{3}/e_{0}]",
                       nx, -250., 250., ny, -250., 250.);

  for (Int_t iy = 1; iy <= ny; ++iy) {
    Double_t yp = h->GetYaxis()->GetBinCenter(iy);
    for (Int_t ix = 1; ix <= nx; ++ix) {
      Double_t xp = h->GetXaxis()->GetBinCenter(ix);

      Float_t r = TMath::Sqrt(xp * xp + yp * yp);
      Float_t phi = TMath::ATan2(yp, xp);

      if (85. <= r && r <= 250.) {

        Float_t sc = GetSpaceChargeDensity(r, phi, z) / fgke0; // in [C/m^3/e0]
        h->SetBinContent(ix, iy, sc);
      } else {
        h->SetBinContent(ix, iy, 0.);
      }
    }
  }

  return h;
}

/// return a simple histogram containing the space charge distribution (input for the calculation)
///
/// \param phi
/// \param nz
/// \param nr
/// \return
TH2F *AliTPCSpaceCharge3DDriftLine::CreateHistogramSCInZR(Float_t phi, Int_t nz, Int_t nr) {

  TH2F *h = CreateTH2F("spaceCharge ", GetTitle(), "z [cm]", "r [cm]", "#rho_{sc} [C/m^{3}/e_{0}]",
                       nz, -250., 250., nr, 85., 250.);
  for (Int_t ir = 1; ir <= nr; ++ir) {
    Float_t r = h->GetYaxis()->GetBinCenter(ir);
    for (Int_t iz = 1; iz <= nz; ++iz) {
      Float_t z = h->GetXaxis()->GetBinCenter(iz);
      if (85. <= r && r <= 250.) {
        Float_t sc = GetSpaceChargeDensity(r, phi, z) / fgke0; // in [C/m^3/e0]
        h->SetBinContent(iz, ir, sc);
      } else {
        h->SetBinContent(iz, ir, 0.);
      }
    }
  }
  return h;
}

/// create the distortion tree on a mesh with granularity given by step
/// return the tree with distortions at given position
/// Map is created on the mesh with given step size
/// type - 0: Call GetDistortion()
///        1: Call GetDistortionIntegralDz()
///
/// \param step
/// \return
TTree *AliTPCSpaceCharge3DDriftLine::CreateDistortionTree(Double_t step)
{
  TTreeSRedirector *pcStream = new TTreeSRedirector(Form("distortion%s.root", GetName()));
  Float_t xyz[3];     // current point
  Float_t dist[3];    // distortion
  Float_t localDist[3];    // distortion
  Float_t corr[3];    // correction
  Float_t xyzDist[3]; // distorted point
  Float_t xyzCorr[3]; // corrected point

  //AliMagF* mag= (AliMagF*)TGeoGlobalMagField::Instance()->GetField();
  //if (!mag) AliError("Magnetic field - not initialized");

  Int_t roc;
  AliTPCParam *tpcParam = new AliTPCParamSR;
  Double_t r, phi, rDist, phiDist, drDist, drPhiDist, rCorr, phiCorr, drCorr, drPhiCorr;
  for (Double_t x = -250; x < 250; x += step) {
    for (Double_t y = -250; y < 250; y += step) {

      r = TMath::Sqrt(x * x + y * y);

      if (tpcParam->GetPadRowRadii(0, 0) > r || r > tpcParam->GetPadRowRadii(36, 95)) continue;

      phi = TMath::ATan2(y, x);

      for (Double_t z = -250; z < 250; z += step) {
        roc = (z > 0) ? 0 : 18;
        xyz[0] = x;
        xyz[1] = y;
        xyz[2] = z;

        GetDistortion(xyz, roc, dist);

        for (Int_t i = 0; i < 3; ++i) {
          xyzDist[i] = xyz[i] + dist[i];
        }

        xyz[0] = r;
        xyz[1] = phi;
        xyz[2] = z;

        GetLocalDistortionCylAC(xyz, roc, localDist);

        GetCorrection(xyzDist, roc, corr);

        for (Int_t i = 0; i < 3; ++i) {
          xyzCorr[i] = xyzDist[i] + corr[i];
        }

        // === r, rPhi + residuals for the distorted point =========================
        rDist = TMath::Sqrt(xyzDist[0] * xyzDist[0] + xyzDist[1] * xyzDist[1]);
        phiDist = TMath::ATan2(xyzDist[1], xyzDist[0]);
        //rDist = xyzDist[0];
        //phiDist = xyzDist[1];

        while ((phiDist - phi) > TMath::Pi()) phiDist -= TMath::TwoPi();
        while ((phiDist - phi) < -TMath::Pi()) phiDist += TMath::TwoPi();

        drDist = rDist - r;
        drPhiDist = (phiDist - phi) * r;

        // === r, rPhi + residuals for the corrected point =========================
        rCorr = TMath::Sqrt(xyzCorr[0] * xyzCorr[0] + xyzCorr[1] * xyzCorr[1]);
        phiCorr = TMath::ATan2(xyzCorr[1], xyzCorr[0]);
        //rCorr = xyzCorr[0];
        //phiCorr = xyzCorr[1];

        while ((phiCorr - phiDist) > TMath::Pi()) phiCorr -= TMath::TwoPi();
        while ((phiCorr - phiDist) < -TMath::Pi()) phiCorr += TMath::TwoPi();

        drCorr = rCorr - rDist;
        drPhiCorr = (phiCorr - phiDist) * r;
        (*pcStream) << "distortion" <<
                    "x=" << x <<           // original position
                    "y=" << y <<
                    "z=" << z <<
                    "r=" << r <<
                    "phi=" << phi <<
                    "xDist=" << xyzDist[0] <<      // distorted position
                    "yDist=" << xyzDist[1] <<
                    "zDist=" << xyzDist[2] <<
                    "rDist=" << rDist <<
                    "phiDist=" << phiDist <<
                    "dxDist=" << dist[0] <<     // distortion
                    "dyDist=" << dist[1] <<
                    "dzDist=" << dist[2] <<
                    "drDist=" << drDist <<
                    "drPhiDist=" << drPhiDist <<
                    "drLocalDist=" << localDist[0] <<
                    "drPhiLocalDist=" << localDist[1] <<
                    "dzLocalDist=" << localDist[2] <<
                    "xCorr=" << xyzCorr[0] <<      // corrected position
                    "yCorr=" << xyzCorr[1] <<
                    "zCorr=" << xyzCorr[2] <<
                    "rCorr=" << rCorr <<
                    "phiCorr=" << phiCorr <<
                    //
                    "dxCorr=" << corr[0] <<     // correction
                    "dyCorr=" << corr[1] <<
                    "dzCorr=" << corr[2] <<
                    "drCorr=" << drCorr <<
                    "drPhiCorr=" << drPhiCorr <<
                    "\n";
      }
    }
  }
  delete pcStream;
  TFile f(Form("distortion%s.root", GetName()));
  TTree *tree = (TTree *) f.Get("distortion");

  return tree;
}

/// create the distortion tree on a mesh with granularity given nRRow, nZColumn and nPhiSlice
///
/// \param nRRowTest
/// \param nZColTest
/// \param nPhiSliceTest
/// \return
TTree *AliTPCSpaceCharge3DDriftLine::CreateDistortionTree(const Int_t nRRowTest, const Int_t nZColTest,
                                                     const Int_t nPhiSliceTest)
{
  TTreeSRedirector *pcStream = new TTreeSRedirector(Form("distortion%s.root", GetName()));

  Int_t nRRow;
  Int_t nZColumn;
  Int_t nPhiSlice;
  Int_t orderInterpolation;
  Int_t rbfKernel;
  Int_t irregularGridSize;

  const Int_t nPoint = nRRowTest * nZColTest * nPhiSliceTest;

  Float_t ofcRadius = AliTPCPoissonSolver::fgkOFCRadius;
  Float_t ifcRadius = AliTPCPoissonSolver::fgkIFCRadius;
  Float_t tpcZ0 = AliTPCPoissonSolver::fgkTPCZ0;

  Float_t dRadius = (AliTPCPoissonSolver::fgkOFCRadius - AliTPCPoissonSolver::fgkIFCRadius) / (nRRowTest - 1);
  Float_t dZ = AliTPCPoissonSolver::fgkTPCZ0 / (nZColTest - 1);
  Float_t dPhi = TMath::TwoPi() / nPhiSliceTest;

  Double_t phi0;
  Double_t z0;
  Double_t r0;

  Float_t drDist;
  Float_t drPhiDist;
  Float_t dzDist;

  Float_t drCorr;
  Float_t drPhiCorr;
  Float_t dzCorr;

  Double_t charge;
  Double_t potential;
  Double_t chargeFormula;
  Double_t potentialFormula;

  nRRow = GetNRRows();
  nZColumn = GetNZColumns();
  nPhiSlice = GetNPhiSlices();
  orderInterpolation = GetInterpolationOrder();
  irregularGridSize = GetIrregularGridSize();
  rbfKernel = GetRBFKernelType();

  Int_t rocNum;

  Float_t point0[] = {0.0, 0.0, 0.0};
  Float_t point1[] = {0.0, 0.0, 0.0};
  Float_t dist[] = {0.0, 0.0, 0.0};
  Float_t corr[] = {0.0, 0.0, 0.0};
  Float_t localDist[3];    // distortion

  Float_t r1;
  Float_t phi1;
  Float_t z1;

  for (Int_t m = 0; m < nPhiSliceTest; m++) {
    for (Int_t i = 0; i < nRRowTest; i++) {
      for (Int_t j = 0; j < nZColTest; j++) {

        phi0 = m * dPhi;
        r0 = ifcRadius + (dRadius * i);
        z0 = dZ * j;
        charge = GetSpaceChargeDensity(r0, phi0, z0);
        potential = GetPotential(r0, phi0, z0);
        if (!(fFormulaChargeRho == NULL))
          chargeFormula = fFormulaChargeRho->Eval(r0, phi0, z0);
        if (!(fFormulaPotentialV == NULL))
          potentialFormula = fFormulaPotentialV->Eval(r0, phi0, z0);
        point0[0] = r0;
        point0[1] = phi0;
        point0[2] = z0;

        rocNum = (z0 > 0.0) ? 0 : 18;

        GetDistortionCyl(point0, rocNum, dist);
        GetLocalDistortionCylAC(point0, rocNum, localDist);
        drDist = dist[0];
        drPhiDist = dist[1];
        dzDist = dist[2];

        r1 = r0 + dist[0];
        phi1 = phi0 + dist[1] / r0;
        z1 = z0 + dist[2];
        point1[0] = r1;
        point1[1] = phi1;
        point1[2] = z1;
        GetCorrectionCyl(point1, rocNum, corr);

        drCorr = corr[0];
        drPhiCorr = corr[1];
        dzCorr = corr[2];

        (*pcStream) << "distortion" <<
                    "z=" << z0 <<
                    "r=" << r0 <<
                    "phi=" << phi0 <<
                    "dzDist=" << dzDist <<
                    "drDist=" << drDist <<
                    "drPhiDist=" << drPhiDist <<
                    "drLocalDist=" << localDist[0] <<
                    "drPhiLocalDist=" << localDist[1] <<
                    "dzLocalDist=" << localDist[2] <<
                    "dzCorr=" << dzCorr <<
                    "drCorr=" << drCorr <<
                    "drPhiCorr=" << drPhiCorr <<
                    "charge=" << charge <<
                    "potential=" << potential <<
                    "chargeFormula=" << chargeFormula <<
                    "potentialFormula=" << potentialFormula <<
                    "\n";
      }
    }
  }

  delete pcStream;
  TFile f(Form("distortion%s.root", GetName()));
  TTree *tree = (TTree *) f.Get("distortion");

  return tree;
}
