/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id$ */

/// \class AliTPCSpaceCharge3DDriftLine
/// \brief This class provides an interface for the usage of the distortion and correction maps calculated by AliTPCSpaceCharge3DCalc
///
/// \author Rifki Sadikin <rifki.sadikin@cern.ch>, Indonesian Institute of Sciences
/// \date Jan 31, 2018

#ifndef ALITPCSPACECHARGE3DDRIFTLINE_H
#define ALITPCSPACECHARGE3DDRIFTLINE_H

#include "AliTPCCorrection.h"
#include "AliTPCSpaceCharge3DCalc.h"

class TH2F;
class TTree;

class AliTPCSpaceCharge3DDriftLine : public AliTPCCorrection, public AliTPCSpaceCharge3DCalc {
public:
  AliTPCSpaceCharge3DDriftLine();
  AliTPCSpaceCharge3DDriftLine(const char *name, const char *title);
  AliTPCSpaceCharge3DDriftLine(const char *name, const char *title, Int_t nRRow, Int_t nZColumn, Int_t nPhiSlice);
  AliTPCSpaceCharge3DDriftLine(const char *name, const char *title, Int_t nRRow, Int_t nZColumn, Int_t nPhiSlice, Int_t interpolationOrder, Int_t irregularGridSize, Int_t rbfKernelType);
  virtual ~AliTPCSpaceCharge3DDriftLine() {};

  void GetDistortion(const Float_t x[], Short_t roc, Float_t dx[]) { AliTPCSpaceCharge3DCalc::GetDistortion(x, roc, dx); }

  void GetCorrection(const Float_t x[], Short_t roc, Float_t dx[]) { AliTPCSpaceCharge3DCalc::GetCorrection(x, roc, dx); }

  void SetOmegaTauT1T2(Float_t omegaTau, Float_t t1, Float_t t2)
  {
    fT1 = t1;
    fT2 = t2;
    const Float_t wt0 = t2 * omegaTau;
    const Float_t c0 = 1. / (1. + wt0 * wt0);
    const Float_t wt1 = t1 * omegaTau;
    const Float_t c1 = wt1 / (1. + wt1 * wt1);
    SetC0C1(c0, c1);
  };

  TH2F *CreateHistogramDistDRInXY(Float_t z, Int_t nx, Int_t ny);

  TH2F *CreateHistogramDistDRPhiInXY(Float_t z, Int_t nx, Int_t ny);

  TH2F *CreateHistogramDistDZInXY(Float_t z, Int_t nx, Int_t ny);

  TH2F *CreateHistogramCorrDRInXY(Float_t z, Int_t nx, Int_t ny);

  TH2F *CreateHistogramCorrDRPhiInXY(Float_t z, Int_t nx, Int_t ny);

  TH2F *CreateHistogramCorrDZInXY(Float_t z, Int_t nx, Int_t ny);

  TH2F *CreateHistogramSCInXY(Float_t z, Int_t nx, Int_t ny);

  TH2F *CreateHistogramSCInZR(Float_t phi, Int_t nz, Int_t nr);

  TTree *CreateDistortionTree(Double_t step);

  TTree *CreateDistortionTree(const Int_t nRRowTest, const Int_t nZColTest, const Int_t nPhiSliceTest);

  void Init();

/// \cond CLASSIMP
ClassDef(AliTPCSpaceCharge3DDriftLine, 1);
/// \endcond
};


#endif //ALIROOT_ALITPCSPACECHARGE3DDRIFTLINE_H
