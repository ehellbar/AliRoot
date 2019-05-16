//**************************************************************************\
//* This file is property of and copyright by the ALICE Project            *\
//* ALICE Experiment at CERN, All rights reserved.                         *\
//*                                                                        *\
//* Primary Authors: Matthias Richter <Matthias.Richter@ift.uib.no>        *\
//*                  for The ALICE HLT Project.                            *\
//*                                                                        *\
//* Permission to use, copy, modify and distribute this software and its   *\
//* documentation strictly for non-commercial purposes is hereby granted   *\
//* without fee, provided that the above copyright notice appears in all   *\
//* copies and that both the copyright notice and this permission notice   *\
//* appear in the supporting documentation. The authors make no claims     *\
//* about the suitability of this software for any purpose. It is          *\
//* provided "as is" without express or implied warranty.                  *\
//**************************************************************************

/// \file GPUChainITS.cxx
/// \author David Rohr

#include "GPUChainITS.h"
#include "GPUReconstructionIncludesITS.h"
#include "DataFormatsITS/TrackITS.h"
#include <algorithm>

using namespace GPUCA_NAMESPACE::gpu;
using namespace o2::its;

GPUChainITS::~GPUChainITS()
{
  mITSTrackerTraits.reset();
  mITSVertexerTraits.reset();
}

GPUChainITS::GPUChainITS(GPUReconstruction* rec) : GPUChain(rec) {}

void GPUChainITS::RegisterPermanentMemoryAndProcessors() { mRec->RegisterGPUProcessor(&processors()->itsFitter, GetRecoStepsGPU() & RecoStep::ITSTracking); }

void GPUChainITS::RegisterGPUProcessors()
{
  if (GetRecoStepsGPU() & RecoStep::ITSTracking) {
    mRec->RegisterGPUDeviceProcessor(&processorsShadow()->itsFitter, &processors()->itsFitter);
  }
}

void GPUChainITS::MemorySize(size_t& gpuMem, size_t& pageLockedHostMem)
{
  gpuMem = 1024 * 1024 * 1024;
  pageLockedHostMem = 1024 * 1024 * 1024;
}

int GPUChainITS::Init()
{
  mRec->GetITSTraits(mITSTrackerTraits, mITSVertexerTraits);
  mITSTrackerTraits->SetRecoChain(this, &GPUChainITS::PrepareAndRunITSTrackFit);
  return 0;
}

int GPUChainITS::PrepareEvent() { return 0; }

int GPUChainITS::Finalize() { return 0; }

int GPUChainITS::RunStandalone() { return 0; }

int GPUChainITS::PrepareAndRunITSTrackFit(std::vector<Road>& roads, std::array<const Cluster*, 7> clusters, std::array<const Cell*, 5> cells, const std::array<std::vector<TrackingFrameInfo>, 7>& tf, std::vector<TrackITS>& tracks)
{
  mRec->PrepareEvent();
  return RunITSTrackFit(roads, clusters, cells, tf, tracks);
}

int GPUChainITS::RunITSTrackFit(std::vector<Road>& roads, std::array<const Cluster*, 7> clusters, std::array<const Cell*, 5> cells, const std::array<std::vector<TrackingFrameInfo>, 7>& tf, std::vector<TrackITS>& tracks)
{
  ActivateThreadContext();
  mRec->SetThreadCounts(RecoStep::ITSTracking);
  bool doGPU = GetRecoStepsGPU() & RecoStep::ITSTracking;
  GPUITSFitter& Fitter = processors()->itsFitter;
  GPUITSFitter& FitterShadow = doGPU ? processorsShadow()->itsFitter : Fitter;

  Fitter.clearMemory();
  Fitter.SetNumberOfRoads(roads.size());
  for (int i = 0; i < 7; i++) {
    Fitter.SetNumberTF(i, tf[i].size());
  }
  Fitter.SetMaxData();
  std::copy(clusters.begin(), clusters.end(), Fitter.clusters());
  std::copy(cells.begin(), cells.end(), Fitter.cells());
  SetupGPUProcessor(&Fitter, true);
  std::copy(roads.begin(), roads.end(), Fitter.roads());
  for (int i = 0; i < 7; i++) {
    std::copy(tf[i].begin(), tf[i].end(), Fitter.trackingFrame()[i]);
  }

  WriteToConstantMemory((char*)&processors()->itsFitter - (char*)processors(), &FitterShadow, sizeof(FitterShadow), 0);
  TransferMemoryResourcesToGPU(&Fitter, 0);
  runKernel<GPUITSFitterKernel>({ BlockCount(), ThreadCount(), 0 }, nullptr, krnlRunRangeNone, krnlEventNone);
  TransferMemoryResourcesToHost(&Fitter, 0);

  SynchronizeGPU();

  for (unsigned int i = 0; i < Fitter.NumberOfTracks(); i++) {
    auto& trkin = Fitter.tracks()[i];

    tracks.emplace_back(TrackITS{ { trkin.X(),
                                    trkin.mAlpha,
                                    { trkin.Par()[0], trkin.Par()[1], trkin.Par()[2], trkin.Par()[3], trkin.Par()[4] },
                                    { trkin.Cov()[0], trkin.Cov()[1], trkin.Cov()[2], trkin.Cov()[3], trkin.Cov()[4], trkin.Cov()[5], trkin.Cov()[6], trkin.Cov()[7], trkin.Cov()[8], trkin.Cov()[9], trkin.Cov()[10], trkin.Cov()[11], trkin.Cov()[12], trkin.Cov()[13], trkin.Cov()[14] } },
                                  (short int)((trkin.NDF() + 5) / 2),
                                  0.139f,
                                  trkin.Chi2(),
                                  0,
                                  { trkin.mOuterParam.X,
                                    trkin.mOuterParam.alpha,
                                    { trkin.mOuterParam.P[0], trkin.mOuterParam.P[1], trkin.mOuterParam.P[2], trkin.mOuterParam.P[3], trkin.mOuterParam.P[4] },
                                    { trkin.mOuterParam.C[0], trkin.mOuterParam.C[1], trkin.mOuterParam.C[2], trkin.mOuterParam.C[3], trkin.mOuterParam.C[4], trkin.mOuterParam.C[5], trkin.mOuterParam.C[6], trkin.mOuterParam.C[7], trkin.mOuterParam.C[8], trkin.mOuterParam.C[9],
                                      trkin.mOuterParam.C[10], trkin.mOuterParam.C[11], trkin.mOuterParam.C[12], trkin.mOuterParam.C[13], trkin.mOuterParam.C[14] } },
                                  { { trkin.mClusters[0], trkin.mClusters[1], trkin.mClusters[2], trkin.mClusters[3], trkin.mClusters[4], trkin.mClusters[5], trkin.mClusters[6] } } });
  }

  ReleaseThreadContext();
  return 0;
}
