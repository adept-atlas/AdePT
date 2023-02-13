// Basic function that takes as arguments the hit definitions (scintillator, period, module) and turns it
// into a basic ID
//
// Initial implementation: returns the module number (0 Barrel, 1 Barrel, 2 EB, 3 plug module, 4 ITC)
// Can be extended to deliver a sliced set of IDs, or use a look up table for the "real" thing
//
// Implementation that can be used both on the device and the host
//
//
#ifndef TILEBASICID_CUH
#define TILEBASICID_CUH

inline
#ifdef __CUDA_ARCH__
__device__ __host__ 
#endif 
int TileBasicID(int scintillator, int period, int module) {
   return module;
}


//
// The numbers of cells that are used in this implementation
//
inline int TileBasicID_Max() {return 5; }

#endif
