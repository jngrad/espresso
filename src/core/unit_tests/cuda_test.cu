/*
 * Copyright (C) 2024 The ESPResSo project
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

#define BOOST_TEST_MODULE cuda test
#define BOOST_TEST_DYN_LINK
#include <boost/test/unit_test.hpp>

#include "cuda/init.hpp"
#include "cuda/utils.cuh"
#include "cuda/utils.hpp"
#include "errorhandling.hpp"

#include "cuda/CudaHostAllocator.hpp"

#include <cuda.h>

#include <boost/test/unit_test.hpp>

#include <cstddef>
#include <optional>
#include <string>

boost::test_tools::assertion_result has_gpu(boost::unit_test::test_unit_id) {
  int n_devices = 0;
  cudaGetDeviceCount(&n_devices);
  if (n_devices > 0) {
    cudaDeviceProp prop;
    cudaGetDeviceProperties(&prop, 0);
    if (prop.major >= 3) {
      return true;
    }
  }
  return false;
}

std::optional<std::string> read_pending_cuda_errors() {
  auto const CU_err = cudaGetLastError();
  if (CU_err != cudaSuccess) {
    auto const message = std::string(cudaGetErrorString(CU_err));
    return {"There is a pending CUDA error: \"" + message + "\""};
  }
  return std::nullopt;
}

void setup() {}
void teardown() {
  auto error = read_pending_cuda_errors();
  BOOST_REQUIRE_MESSAGE(not error.has_value(), error.value_or(""));
}

namespace Testing::non_sticky_cuda_error {

/** @brief Trigger a non-sticky CUDA error for testing purposes. */
void trigger() { cudaSetDevice(-1); }

/** @brief Clear a non-sticky CUDA error raised by @ref trigger. */
void clear() {
  auto const error_code = cudaGetLastError();
  BOOST_REQUIRE_MESSAGE(error_code == cudaSuccess or
                            error_code == cudaErrorInvalidDevice,
                        "An unexpected CUDA error was pending!");
}

} // namespace Testing::non_sticky_cuda_error

static auto fixture = boost::unit_test::fixture(&setup, &teardown);

BOOST_AUTO_TEST_SUITE(suite, *boost::unit_test::precondition(has_gpu))

BOOST_AUTO_TEST_CASE(gpu_fixture, *fixture) {
  {
    auto error = read_pending_cuda_errors();
    BOOST_REQUIRE(not error.has_value());
  }
  {
    // check we can raise and clear non-sticky CUDA errors
    Testing::non_sticky_cuda_error::trigger();
    Testing::non_sticky_cuda_error::clear();
    auto error = read_pending_cuda_errors();
    BOOST_REQUIRE(not error.has_value());
  }
  {
    // check fixture can handle the default non-sticky CUDA error
    Testing::non_sticky_cuda_error::trigger();
    auto ref_what{"There is a pending CUDA error: \"invalid device ordinal\""};
    auto error = read_pending_cuda_errors();
    BOOST_REQUIRE(error.has_value());
    BOOST_REQUIRE_EQUAL(error.value(), ref_what);
    // sticky error should have been cleared
    error = read_pending_cuda_errors();
    BOOST_REQUIRE(not error.has_value());
  }
  {
    // check fixture can handle a custom non-sticky CUDA error
    cudaMallocHost(nullptr, std::size_t(0u));
    auto ref_what{"There is a pending CUDA error: \"invalid argument\""};
    auto error = read_pending_cuda_errors();
    BOOST_REQUIRE(error.has_value());
    BOOST_REQUIRE_EQUAL(error.value(), ref_what);
    // sticky error should have been cleared
    error = read_pending_cuda_errors();
    BOOST_REQUIRE(not error.has_value());
  }
}

static int fatal_error_counter = 0;
static void increment_counter() noexcept { ++fatal_error_counter; }

BOOST_AUTO_TEST_CASE(gpu_interface, *fixture) {
  fatal_error_counter = 0;
  auto local_error_counter = 0;
  {
    std::string const what = "message 1";
    try {
      throw cuda_fatal_error(what);
    } catch (cuda_fatal_error &err) {
      BOOST_CHECK_EQUAL(err.what(), what);
      BOOST_CHECK_EQUAL(err.get_terminate(), &errexit);
      err.set_terminate(nullptr);
      BOOST_CHECK_EQUAL(err.get_terminate(), nullptr);
      err.set_terminate(increment_counter);
      BOOST_CHECK_EQUAL(err.get_terminate(), &increment_counter);
      BOOST_CHECK_EQUAL(fatal_error_counter, local_error_counter);
    }
    ++local_error_counter;
    BOOST_REQUIRE_EQUAL(fatal_error_counter, local_error_counter);
  }
  {
    auto error_caught = false;
    auto const block = dim3{1, 2, 3};
    auto const grid = dim3{4, 5, 6};
    // should not throw
    cuda_check_errors_exit(block, grid, "", "", 0u);
    try {
      Testing::non_sticky_cuda_error::trigger();
      // should clear the CUDA error flag and throw a fatal error
      cuda_check_errors_exit(block, grid, "cudaSetDevice()", "filename.cu", 4u);
    } catch (cuda_fatal_error &err) {
      error_caught = true;
      err.set_terminate(increment_counter);
      std::string const what =
          "CUDA error: \"invalid device ordinal\" while calling "
          "cudaSetDevice() with block: <1,2,3>, grid: <4,5,6> in filename.cu:4";
      BOOST_CHECK_EQUAL(fatal_error_counter, local_error_counter);
      BOOST_CHECK_EQUAL(err.what(), what);
      BOOST_CHECK_EQUAL(cudaGetLastError(), cudaSuccess);
    }
    ++local_error_counter;
    BOOST_REQUIRE(error_caught);
    BOOST_REQUIRE_EQUAL(fatal_error_counter, local_error_counter);
  }
  {
    auto error_caught = false;
    // should not throw
    cuda_safe_mem_exit(cudaSuccess, "", 0u);
    try {
      Testing::non_sticky_cuda_error::trigger();
      // should throw
      cuda_safe_mem_exit(cudaSuccess, "filename.cu", 4u);
    } catch (cuda_fatal_error &err) {
      error_caught = true;
      err.set_terminate(increment_counter);
      std::string const what =
          "CUDA error: \"invalid device ordinal\" in filename.cu:4. Error "
          "found during memory operation. Possibly however from a failed "
          "operation before the memory operation";
      BOOST_CHECK_EQUAL(fatal_error_counter, local_error_counter);
      BOOST_CHECK_EQUAL(err.what(), what);
    }
    ++local_error_counter;
    BOOST_REQUIRE(error_caught);
    BOOST_REQUIRE_EQUAL(fatal_error_counter, local_error_counter);
  }
  {
    auto error_caught = false;
    try {
      cuda_safe_mem_exit(cudaErrorNotPermitted, "filename.cu", 4u);
    } catch (cuda_fatal_error &err) {
      error_caught = true;
      err.set_terminate(increment_counter);
      std::string const what = "CUDA error: \"operation not permitted\" during "
                               "memory operation in filename.cu:4";
      BOOST_CHECK_EQUAL(fatal_error_counter, local_error_counter);
      BOOST_CHECK_EQUAL(err.what(), what);
    }
    ++local_error_counter;
    BOOST_REQUIRE(error_caught);
    BOOST_REQUIRE_EQUAL(fatal_error_counter, local_error_counter);
  }
  {
    auto error_caught = false;
    try {
      cuda_safe_mem_exit(cudaErrorInvalidValue, "function_name()", 4u);
    } catch (cuda_fatal_error &err) {
      error_caught = true;
      err.set_terminate(increment_counter);
      std::string const what =
          "CUDA error: \"invalid argument\" during memory operation in "
          "function_name():4. You may have tried to allocate zero memory";
      BOOST_CHECK_EQUAL(fatal_error_counter, local_error_counter);
      BOOST_CHECK_EQUAL(err.what(), what);
    }
    ++local_error_counter;
    BOOST_REQUIRE(error_caught);
    BOOST_REQUIRE_EQUAL(fatal_error_counter, local_error_counter);
  }
  {
    BOOST_REQUIRE_EQUAL(stream[0], nullptr);
    auto error_caught = false;
    cuda_init(); // allocate
    BOOST_REQUIRE_NE(stream[0], nullptr);
    cuda_set_device(0); // reallocate, may or may not result in the same pointer
    BOOST_REQUIRE_NE(stream[0], nullptr);
    auto const old_stream = stream[0];
    try {
      cuda_set_device(-1); // fail to reallocate, pointer remains the same
    } catch (cuda_runtime_error_cuda const &err) {
      error_caught = true;
      std::string const what = "CUDA error: invalid device ordinal";
      BOOST_CHECK_EQUAL(err.what(), what);
    }
    BOOST_REQUIRE(error_caught);
    BOOST_REQUIRE_EQUAL(stream[0], old_stream);
  }
  {
    BOOST_REQUIRE_GE(cuda_get_n_gpus(), 1);
    char gpu_name_buffer[260] = {'\0'};
    cuda_get_gpu_name(0, gpu_name_buffer);
    for (int i = 255; i < 260; ++i) {
      BOOST_REQUIRE_EQUAL(gpu_name_buffer[i], '\0');
    }
  }
}

BOOST_AUTO_TEST_SUITE_END()
