/*
 * Copyright (C) 2020 The ESPResSo project
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

/* Unit tests for the FFT code. */

#define BOOST_TEST_NO_MAIN
#define BOOST_TEST_MODULE FFT test
#define BOOST_TEST_ALTERNATIVE_INIT_API
#define BOOST_TEST_DYN_LINK
#include <boost/test/unit_test.hpp>

#include <boost/mpi.hpp>

#include "communication.hpp"
#include "electrostatics_magnetostatics/fft.hpp"

#include <utils/Vector.hpp>
#include <utils/constants.hpp>
#include <utils/index.hpp>
#include <utils/math/int_pow.hpp>
#include <utils/mpi/cart_comm.hpp>

#include <algorithm>
#include <cmath>
#include <limits>
#include <numeric>

#include <iostream>

BOOST_AUTO_TEST_CASE(fft_calc) {
  extern boost::mpi::communicator comm_cart;
  int n_nodes;
  MPI_Comm_size(MPI_COMM_WORLD, &n_nodes);
  auto const dims = Utils::Mpi::dims_create<3>(n_nodes);
  comm_cart = Utils::Mpi::cart_create(comm_cart, dims, false);

  constexpr auto margin = 4;
  constexpr auto size_in = 32;
  constexpr auto size_out = size_in + 2 * margin;
  constexpr auto volume_in = Utils::int_pow<3>(size_in);
  constexpr auto volume_out = Utils::int_pow<3>(size_out);
  auto const global_mesh_dim = Utils::Vector3i::broadcast(size_in);
  auto const global_mesh_off = Utils::Vector3d::broadcast(0.5);
  auto const ca_mesh_dim = Utils::Vector3i::broadcast(size_out);
  auto const ca_mesh_margin = Utils::VectorXi<6>::broadcast(margin);
  int ks_pnum;
  fft_data_struct fft;
  auto const retval =
      fft_init(ca_mesh_dim, ca_mesh_margin.data(), global_mesh_dim.data(),
               global_mesh_off.data(), ks_pnum, fft, dims, comm_cart);

  BOOST_CHECK(fft.init_tag);
  BOOST_CHECK_EQUAL(fft.max_mesh_size, retval);
  BOOST_CHECK_EQUAL(fft.max_mesh_size, fft.data_buf.size());
  BOOST_CHECK_GE(fft.max_mesh_size, volume_out);
  BOOST_CHECK_EQUAL(fft.plan[1].new_size, volume_in);

  /* cosines */
  {
    auto data = std::vector<double>(fft.data_buf.size());
    constexpr double wave_vector = 2 * Utils::pi() / size_in;
    Utils::Vector3i pos;
    for (int i = 0; i < size_out; i++) {
      for (int j = 0; j < size_out; j++) {
        for (int k = 0; k < size_out; k++) {
          auto const index = Utils::get_linear_index(i, j, k, ca_mesh_dim);
          data[index] = std::sin(wave_vector * i * 2) +
                        std::sin(wave_vector * j * 2) +
                        std::sin(wave_vector * k * 2);
        }
      }
    }
    fft_perform_forw(data.data(), fft, comm_cart);
    constexpr double ref_val = Utils::int_pow<2>(128);
    std::array<Utils::Vector3i, 6> coord_freq = {{{4, 0, 0},
                                                  {20, 1, 0},
                                                  {8, 3, 0},
                                                  {0, 8, 1},
                                                  {16, 22, 2},
                                                  {0, 16, 38}}};
    for (auto const &coord : coord_freq) {
      auto const pos = Utils::get_linear_index(coord, ca_mesh_dim);
      BOOST_CHECK_EQUAL(data[pos], ref_val);
      data[pos] = 0.;
    }
    auto const highest_non_zero = std::abs(
        *std::max_element(data.begin(), data.end(), [](double a, double b) {
          return (std::abs(a) < std::abs(b));
        }));
    BOOST_CHECK_LT(highest_non_zero, 1e-10);
  }

  /* Dirac function */
  {
    auto data = std::vector<double>(fft.data_buf.size());
    constexpr double height = 2;
    auto const center = Utils::get_linear_index(ca_mesh_dim / 2, ca_mesh_dim);
    data[center] = height;
    fft_perform_forw(data.data(), fft, comm_cart);
    auto const count_positive_heights =
        std::accumulate(data.begin(), data.end(), 0, [](int acc, double x) {
          return acc + (std::abs(x - height) <= 1e-10);
        });
    auto const count_negative_heights =
        std::accumulate(data.begin(), data.end(), 0, [](int acc, double x) {
          return acc + (std::abs(x + height) <= 1e-10);
        });
    BOOST_CHECK_EQUAL(count_positive_heights, volume_in / 2);
    BOOST_CHECK_EQUAL(count_negative_heights, volume_in / 2);
  }

  /* Compare with fftw3 */
  {
    auto in = std::vector<fftw_complex>(volume_out);
    auto out = std::vector<fftw_complex>(volume_out);
    auto const plan = fftw_plan_dft_3d(size_out, size_out, size_out, in.data(),
                                       out.data(), FFTW_FORWARD, FFTW_MEASURE);
    constexpr double wave_vector = 2 * Utils::pi() / size_out;
    Utils::Vector3i pos;
    for (int i = 0; i < size_out; i++) {
      for (int j = 0; j < size_out; j++) {
        for (int k = 0; k < size_out; k++) {
          auto const index = Utils::get_linear_index(i, j, k, ca_mesh_dim);
          in[index][0] = std::sin(wave_vector * i * 2) +
                         std::sin(wave_vector * j * 2) +
                         std::sin(wave_vector * k * 2);
          in[index][1] = 0;
        }
      }
    }
    fftw_execute(plan);
    for (int i = 0; i < size_out; i++) {
      for (int j = 0; j < size_out; j++) {
        for (int k = 0; k < size_out; k++) {
          auto const index = Utils::get_linear_index(i, j, k, ca_mesh_dim);
          if (std::abs(out[index][0]) >= 1e-9)
            std::cout << "i = " << i << " j = " << j << " k = " << k
                      << " val = " << out[index][0] << std::endl;
        }
      }
    }
  }
}

int main(int argc, char **argv) {
  boost::mpi::environment mpi_env(argc, argv);
  return boost::unit_test::unit_test_main(init_unit_test, argc, argv);
}
