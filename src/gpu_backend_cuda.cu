#include "gpu_backend.h"
#include "gk.h"

#include <cuda_runtime.h>

int gpuBackendAvailable(void) {
  int deviceCount = 0;
  cudaError_t err = cudaGetDeviceCount(&deviceCount);
  if (err != cudaSuccess) {
    return 0;
  }
  return (deviceCount > 0) ? 1 : 0;
}

int gpuNearfieldApply(ssystem *sys, double alpha, const double *sgm, double *pot) {
  (void)sys;
  (void)alpha;
  (void)sgm;
  (void)pot;

  /*
   * Stage-1 GPU integration:
   * runtime backend wiring is in place, but nearfield compute still falls back
   * to CPU until panel quadrature and neighbor-list kernels are ported.
   */
  return 0;
}
