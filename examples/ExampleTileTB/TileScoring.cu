// SPDX-FileCopyrightText: 2022 CERN
// SPDX-License-Identifier: Apache-2.0

#include "TileScoring.h"
#include "AdeptIntegration.h"

#include <CopCore/Global.h>
#include <CopCore/PhysicalConstants.h>

#include "Track.cuh" // not nice - we expose the track model here, interface of DepositEnergy to be changed
#include "TileBasicID.cuh"

#include <iostream>
#include <iomanip>
#include <stdio.h>

TileScoring *TileScoring::InitializeOnGPU()
{
  fAuxData_dev = AdeptIntegration::VolAuxArray::GetInstance().fAuxData_dev;
  // Allocate memory to score charged track length and energy deposit per volume.
  COPCORE_CUDA_CHECK(cudaMalloc(&fChargedTrackLength_dev, sizeof(double) * fNumSensitive));
  COPCORE_CUDA_CHECK(cudaMemset(fChargedTrackLength_dev, 0, sizeof(double) * fNumSensitive));
  COPCORE_CUDA_CHECK(cudaMalloc(&fEnergyDeposit_dev, sizeof(double) * fNumSensitive));
  COPCORE_CUDA_CHECK(cudaMemset(fEnergyDeposit_dev, 0, sizeof(double) * fNumSensitive));

  // Allocate and initialize scoring and statistics.
  COPCORE_CUDA_CHECK(cudaMalloc(&fGlobalScoring_dev, sizeof(GlobalScoring)));
  COPCORE_CUDA_CHECK(cudaMemset(fGlobalScoring_dev, 0, sizeof(GlobalScoring)));

  ScoringPerVolume scoringPerVolume_devPtrs;
  scoringPerVolume_devPtrs.chargedTrackLength = fChargedTrackLength_dev;
  scoringPerVolume_devPtrs.energyDeposit      = fEnergyDeposit_dev;
  COPCORE_CUDA_CHECK(cudaMalloc(&fScoringPerVolume_dev, sizeof(ScoringPerVolume)));
  COPCORE_CUDA_CHECK(
      cudaMemcpy(fScoringPerVolume_dev, &scoringPerVolume_devPtrs, sizeof(ScoringPerVolume), cudaMemcpyHostToDevice));
  // Now allocate space for the TileScoring placeholder on device and only copy the device pointers of components
  TileScoring *TileScoring_dev = nullptr;
  COPCORE_CUDA_CHECK(cudaMalloc(&TileScoring_dev, sizeof(TileScoring)));
  COPCORE_CUDA_CHECK(cudaMemcpy(TileScoring_dev, this, sizeof(TileScoring),  cudaMemcpyHostToDevice));
  return TileScoring_dev;
}

void TileScoring::FreeGPU()
{
  // Free resources.
  COPCORE_CUDA_CHECK(cudaFree(fChargedTrackLength_dev));
  COPCORE_CUDA_CHECK(cudaFree(fEnergyDeposit_dev));

  COPCORE_CUDA_CHECK(cudaFree(fGlobalScoring_dev));
  COPCORE_CUDA_CHECK(cudaFree(fScoringPerVolume_dev));
}

void TileScoring::ClearGPU()
{
  // Clear the device hits content
  COPCORE_CUDA_CHECK(cudaMemset(fGlobalScoring_dev, 0, sizeof(GlobalScoring)));
  COPCORE_CUDA_CHECK(cudaMemset(fChargedTrackLength_dev, 0, sizeof(double) * fNumSensitive));
  COPCORE_CUDA_CHECK(cudaMemset(fEnergyDeposit_dev, 0, sizeof(double) * fNumSensitive));
}

void TileScoring::CopyHitsToHost()
{
  // Transfer back scoring.
  COPCORE_CUDA_CHECK(cudaMemcpy(&fGlobalScoring, fGlobalScoring_dev, sizeof(GlobalScoring), cudaMemcpyDeviceToHost));

  // Transfer back the scoring per volume (charged track length and energy deposit).
  COPCORE_CUDA_CHECK(cudaMemcpy(fScoringPerVolume.chargedTrackLength, fChargedTrackLength_dev,
                                sizeof(double) * fNumSensitive, cudaMemcpyDeviceToHost));
  COPCORE_CUDA_CHECK(cudaMemcpy(fScoringPerVolume.energyDeposit, fEnergyDeposit_dev, sizeof(double) * fNumSensitive,
                                cudaMemcpyDeviceToHost));
}

__device__ void TileScoring::Score(vecgeom::NavStateIndex const &crt_state, int charge, double geomStep, double edep)
{
  assert(fGlobalScoring_dev && "Scoring not initialized on device");
  auto nLevels=crt_state.GetLevel();
 // auto volume  = crt_state.At(nLevels);
 // auto volume_2  = crt_state.At(nLevels-2);
 // auto volume_5  = crt_state.At(nLevels-5);
  int scintillator_copy_no = crt_state.At(nLevels)->GetCopyNo();
  int period_copy_no=crt_state.At(nLevels-2)->GetCopyNo();
  int module_copy_no=crt_state.At(nLevels-5)->GetCopyNo();
  //
  // Mapping of module numbers. This is done by checking the volume name in the Tile SD
  // Since we can't access this info, we get the copy no and use a similar logic
  // At some point the copy no logic between cuda and C++ needs to be aligned
  int module_no(0);
  switch(module_copy_no) {
    case 1: // Barrel 1. ID=0
      module_no=0;
      break;
    case 2: // Barrel 2. ID=1
      module_no=1;
      break;
    case 101: // EBarrel ID=2
      module_no=2; 
      break;
    case 1001: // PlugToModule ID=3
      module_no=3; 
      break;
    case 2001: // ITC ID=4
      module_no=4; 
      break;
  }   
// Some check printouts - remove comments if you want them!
//  printf("In Scoring: Number of levels %d \n",crt_state.GetLevel() );
//  printf("In Scoring: (Current level) Copy no at level %d = %d \n", nLevels, scintillator_copy_no  );
//  printf("In Scoring: (2 levels above) Copy no at level %d = %d \n", nLevels-2, period_copy_no  );
//  printf("In Scoring: (5 levels above) Copy no at level %d = %d \n", nLevels-5, module_no  );

  int volumeID = TileBasicID(scintillator_copy_no, period_copy_no, module_no);

// This is a host function only
//  printf("In Scoring: Name at level %d = %s \n", nLevels, volume->GetName()  );

//  crt_state.Print();
//   int volumeID = volume->id();
  int charged  = abs(charge);

//  int lvolID = volume->GetLogicalVolume()->id();

  // Add to charged track length, global energy deposit and deposit per volume
  atomicAdd(&fScoringPerVolume_dev->chargedTrackLength[volumeID], charged * geomStep);
  atomicAdd(&fGlobalScoring_dev->energyDeposit, edep);
  atomicAdd(&fScoringPerVolume_dev->energyDeposit[volumeID], edep);
}

__device__ void TileScoring::AccountHit()
{
  // Increment hit counter
  atomicAdd(&fGlobalScoring_dev->hits, 1);
}

__device__ void TileScoring::AccountChargedStep(int charge)
{
  // Increase counters for charged/neutral steps
  int charged = abs(charge);
  // Increment global number of steps
  atomicAdd(&fGlobalScoring_dev->chargedSteps, charged);
  atomicAdd(&fGlobalScoring_dev->neutralSteps, 1 - charged);
}

__device__ void TileScoring::AccountProduced(int num_ele, int num_pos, int num_gam)
{
  // Increment number of secondaries
  atomicAdd(&fGlobalScoring_dev->numElectrons, num_ele);
  atomicAdd(&fGlobalScoring_dev->numPositrons, num_pos);
  atomicAdd(&fGlobalScoring_dev->numGammas, num_gam);
}
