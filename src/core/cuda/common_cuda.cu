/*
 * Copyright (C) 2010-2022 The ESPResSo project
 *
 * This file is part of ESPResSo.
 *
 * ESPResSo is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * ESPResSo is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program.  If not, see <http://www.gnu.org/licenses/>.
 */

#include "errorhandling.hpp"

#include "utils.cuh"
#include "utils.hpp"

#include <cuda.h>
#include <cuda_runtime.h>

#include <cstdio>
#include <cstdlib>
#include <sstream>
#include <string>
#include <utility>

cudaStream_t stream[1];

void cuda_check_errors_exit(const dim3 &block, const dim3 &grid,
                            const char *function, const char *file,
                            unsigned int line) {
  cudaError_t CU_err = cudaGetLastError();
  if (CU_err != cudaSuccess) {
    std::stringstream message;
    message << "CUDA error: \"" << cudaGetErrorString(CU_err) << "\" "
            << "calling " << function << " with "
            << "block: [" << block.x << "," << block.y << "," << block.z
            << "], "
            << "grid: [" << grid.x << "," << grid.y << "," << grid.z << "] "
            << "in " << file << ":" << line;
    throw cuda_fatal_error(message.str());
  }
}

void cuda_safe_mem_exit(cudaError_t CU_err, const char *file,
                        unsigned int line) {
  if (CU_err != cudaSuccess) {
    std::stringstream message;
    message << "CUDA error: \"" << cudaGetErrorString(CU_err)
            << "\" during memory operation in " << file << ":" << line;
    if (CU_err == cudaErrorInvalidValue)
      message << ". You may have tried to allocate zero memory";
    throw cuda_fatal_error(message.str());
  }
  {
    CU_err = cudaGetLastError();
    if (CU_err != cudaSuccess) {
      std::stringstream message;
      message << "CUDA error: \"" << cudaGetErrorString(CU_err) << "\" in "
              << file << ":" << line << ". Error found during memory operation"
              << ". Possibly however from a failed operation before the memory "
                 "operation";
      throw cuda_fatal_error(message.str());
    }
  }
}

cuda_fatal_error::cuda_fatal_error(std::string msg)
    : m_msg(std::move(msg)), m_terminate_handler(&errexit) {}

void cuda_fatal_error::terminate() noexcept {
  ((m_terminate_handler == nullptr) ? &std::abort : m_terminate_handler)();
}
