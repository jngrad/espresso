/*
 * Copyright (C) 2010-2026 The ESPResSo project
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

#pragma once

/**
 * @file
 * Shared waLBerla scaffolding for both single-component
 * (LBWalberlaImplSingleComponent) and color-gradient
 * (LBWalberlaImplColorGradient) lattice-Boltzmann leaves.
 * Used as a CRTP base; the leaf's typed field IDs are reached via
 * static_cast<Derived &>(*this).
 */

#include <blockforest/Initialization.h>
#include <blockforest/StructuredBlockForest.h>
#include <domain_decomposition/BlockDataID.h>
#include <domain_decomposition/IBlock.h>
#include <field/AddToStorage.h>
#include <field/vtk/FlagFieldCellFilter.h>
#include <field/vtk/VTKWriter.h>
#include <stencil/D3Q19.h>
#include <stencil/D3Q27.h>
#include <waLBerlaDefinitions.h>
#if defined(__CUDACC__) and defined(WALBERLA_BUILD_WITH_CUDA)
#include <gpu/AddGPUFieldToStorage.h>
#include <gpu/HostFieldAllocator.h>
#endif

#include "../BoundaryHandling.hpp"
#include "../BoundaryPackInfo.hpp"
#include "../utils/boundary.hpp"
#include "../utils/types_conversion.hpp"
#include "../electrokinetics/generated_kernels/EK_FieldAccessors_double_precision.h"
#include "../electrokinetics/generated_kernels/EK_FieldAccessors_single_precision.h"
#include "lb_fields.hpp"
#include "lb_kernels.hpp"
#if defined(__CUDACC__) and defined(WALBERLA_BUILD_WITH_CUDA)
// the GPU accessors are called with qualified names in the class templates
// below, so they must be declared before the class definitions
#include "../electrokinetics/generated_kernels/EK_FieldAccessors_double_precision_CUDA.cuh"
#include "../electrokinetics/generated_kernels/EK_FieldAccessors_single_precision_CUDA.cuh"
#include "lb_fields.cuh"
#include "lb_kernels.cuh"
#endif

#include <walberla_bridge/Architecture.hpp>
#include <walberla_bridge/BlockAndCell.hpp>
#include <walberla_bridge/LatticeWalberla.hpp>
#include <walberla_bridge/lattice_boltzmann/LBWalberlaBase.hpp>
#include <walberla_bridge/utils/ResourceManager.hpp>
#include <walberla_bridge/walberla_init.hpp>

#include <utils/Vector.hpp>
#include <utils/interpolation/bspline_3d.hpp>

#include <algorithm>
#include <array>
#include <bitset>
#include <cassert>
#include <cstddef>
#include <cstdint>
#include <functional>
#include <iostream>
#include <iterator>
#include <limits>
#include <memory>
#include <optional>
#include <stdexcept>
#include <string>
#include <type_traits>
#include <utility>
#include <variant>
#include <vector>

namespace walberla {

/**
 * @brief Exception for accessing a lattice node outside the local domain
 *  and ghost layers during B-spline interpolation.
 */
class interpolation_illegal_access : public std::runtime_error {
public:
  interpolation_illegal_access(std::string const &field,
                               Utils::Vector3d const &pos,
                               std::array<int, 3> const &node, double weight)
      : std::runtime_error("Access to LB " + field + " field failed") {
    std::cerr << "pos [" << pos << "], node [" << Utils::Vector3i(node)
              << "], weight " << weight << "\n";
  }
};

inline void interpolate_bspline_at_pos(Utils::Vector3d const &pos,
                                       auto const &&f) {
  Utils::Interpolation::bspline_3d<2>(
      pos, f, Utils::Vector3d::broadcast(1.), // grid spacing
      Utils::Vector3d::broadcast(.5));        // offset
}

/**
 * @brief CRTP mixin providing the shared waLBerla scaffolding for
 *        @ref LBWalberlaImpl variants.
 *
 * The leaf (@p Derived) supplies typed BlockDataIDs for PDF fields and any
 * additional, layout-dependent fields. The mixin owns the lattice handle,
 * UBB boundary scaffolding, off-lattice read-side interpolation, and the
 * non-CG ghost-comm scaffolding.
 */
template <class Derived, typename FloatType, lbmpy::Arch Architecture,
         typename BoundaryKernel =
             typename detail::KernelTrait<FloatType, Architecture>::DynamicUBB>
class LBWalberlaCommon : public virtual LBWalberlaBase {
#if not defined(WALBERLA_BUILD_WITH_CUDA)
  static_assert(Architecture != lbmpy::Arch::GPU,
                "waLBerla was compiled without CUDA support");
#endif
protected:
  Derived &derived() noexcept { return static_cast<Derived &>(*this); }
  Derived const &derived() const noexcept {
    return static_cast<Derived const &>(*this);
  }

  // ---- Types & Constants ----

  using Kernels = detail::KernelTrait<FloatType, Architecture>;
  using BoundaryModel =
      BoundaryHandling<FloatType, Vector3<FloatType>, BoundaryKernel>;

public:
  /** @brief Stencil for collision and streaming operations. */
  using Stencil = stencil::D3Q19;
  /** @brief Stencil for ghost communication (includes domain corners). */
  using StencilFull = stencil::D3Q27;
  /** @brief Lattice model (e.g. blockforest). */
  using BlockStorage = LatticeWalberla::Lattice_T;

protected:
  // "underlying" field types (`GPUField` has no f-size info at compile time)
  using _PdfField = FieldTrait<FloatType, Stencil>::PdfField;
  using _VectorField = FieldTrait<FloatType, Stencil>::VectorField;
  using _ScalarField = FieldTrait<FloatType, Stencil>::ScalarField;

public:
  using ScalarField = FieldTrait<FloatType, Stencil, Architecture>::ScalarField;
  using PdfField = FieldTrait<FloatType, Stencil, Architecture>::PdfField;
  using VectorField = FieldTrait<FloatType, Stencil, Architecture>::VectorField;
  using FlagField = BoundaryModel::FlagField;
#if defined(__CUDACC__) and defined(WALBERLA_BUILD_WITH_CUDA)
  using GPUField = gpu::GPUField<FloatType>;
  using PdfFieldCpu =
      FieldTrait<FloatType, Stencil, lbmpy::Arch::CPU>::PdfField;
  using ScalarFieldCpu =
      FieldTrait<FloatType, Stencil, lbmpy::Arch::CPU>::ScalarField;
  using VectorFieldCpu =
      FieldTrait<FloatType, Stencil, lbmpy::Arch::CPU>::VectorField;
#endif

  struct GhostComm {
    /** @brief Ghost communication operations. */
    enum GhostCommFlags : unsigned {
      PDF, ///< PDFs communication
      VEL, ///< velocities communication
      LAF, ///< last applied forces communication
      UBB, ///< boundaries communication
      PHI, ///< phasefield communication (two_component)
      RHO, ///< per-component density communication (two_component)
      SIZE
    };
  };

protected:
  /**
   * @brief Full communicator.
   * We use the D3Q27 directions to update cells along the diagonals during
   * a full ghost communication. This is needed to properly update the corners
   * of the ghost layer when setting cell velocities or populations.
   */
  using RegularFullCommunicator =
      FieldTrait<FloatType, Stencil,
                 Architecture>::template RegularCommScheme<stencil::D3Q27>;
  using BoundaryFullCommunicator =
      FieldTrait<FloatType, Stencil,
                 Architecture>::template BoundaryCommScheme<stencil::D3Q27>;
  /**
   * @brief Regular communicator.
   * We use the same directions as the stencil during integration.
   */
  using PDFStreamingCommunicator =
      FieldTrait<FloatType, Stencil,
                 Architecture>::template RegularCommScheme<Stencil>;
  template <class Field>
  using PackInfo =
      FieldTrait<FloatType, Stencil, Architecture>::template PackInfo<Field>;

protected:
  // ---- Member Variables (shared infrastructure) ----

  FloatType m_density;
  double m_zc_to_md; // zero-centered conversion factor to MD units
  double m_zc_to_lb; // zero-centered conversion factor to LB units
  unsigned int m_seed{0u};

  // lattice
  std::shared_ptr<LatticeWalberla> m_lattice;

  // Shared block data IDs (single-instance, layout-independent)
  BlockDataID m_flag_field_id;
  BlockDataID m_velocity_field_id;
  BlockDataID m_vel_tmp_field_id;

  /** Flag for boundary cells. */
  FlagUID const Boundary_flag{"boundary"};
  bool m_has_boundaries{false};

  // boundaries: one BoundaryModel per PDF field (size 1 for SC, size 2 for
  // CG); all instances share the same flag field and store
  // the same wall velocity value per boundary cell
  std::vector<std::shared_ptr<BoundaryModel>> m_boundary;

  // communicators (single-component scaffolding; CG comms live in the leaf)
  std::shared_ptr<BoundaryFullCommunicator> m_boundary_communicator;
  std::shared_ptr<RegularFullCommunicator> m_full_communicator;
  std::shared_ptr<RegularFullCommunicator> m_pdf_communicator;
  std::shared_ptr<RegularFullCommunicator> m_vel_communicator;
  std::shared_ptr<RegularFullCommunicator> m_laf_communicator;
  std::shared_ptr<PDFStreamingCommunicator> m_pdf_streaming_communicator;
  std::bitset<GhostComm::SIZE> m_pending_ghost_comm;
  ResourceObserver m_mpi_cart_comm_observer;

protected:
  LBWalberlaCommon(std::shared_ptr<LatticeWalberla> lattice, FloatType density)
      : m_density(density), m_zc_to_md(static_cast<double>(density)),
        m_zc_to_lb(1. / static_cast<double>(density)),
        m_lattice(std::move(lattice)),
        m_mpi_cart_comm_observer(get_mpi_cart_comm_observer()) {}

public:
  template <typename T> FloatType FloatType_c(T t) const {
    return numeric_cast<FloatType>(t);
  }

  [[nodiscard]] std::size_t stencil_size() const noexcept override {
    return static_cast<std::size_t>(Stencil::Size);
  }

  [[nodiscard]] bool is_double_precision() const noexcept override {
    return std::is_same_v<FloatType, double>;
  }

  [[nodiscard]] bool is_gpu() const noexcept override {
    return Architecture == lbmpy::Arch::GPU;
  }

  [[nodiscard]] LatticeWalberla const &get_lattice() const noexcept override {
    return *m_lattice;
  }

  [[nodiscard]] std::size_t get_velocity_field_id() const noexcept override {
    return m_velocity_field_id;
  }

  [[nodiscard]] std::size_t get_force_field_id() const noexcept override {
    // component a only on two-component (color-gradient) LB; see the
    // doc comment on LBWalberlaBase::get_force_field_id
    return derived().m_force_to_be_applied_id[0];
  }

  [[nodiscard]] double get_density() const noexcept override {
    return static_cast<double>(m_density);
  }

  // ---- Global external force ----

  void set_external_force(Utils::Vector3d const &ext_force) override {
    derived().m_reset_force->set_ext_force(zero_centered_to_lb(ext_force));
  }

  [[nodiscard]] Utils::Vector3d get_external_force() const noexcept override {
    return zero_centered_to_md(derived().m_reset_force->get_ext_force());
  }

  // ---- RNG seed / state (thermalized collision model) ----

  [[nodiscard]] unsigned int get_seed() const noexcept override {
    return m_seed;
  }

  [[nodiscard]] std::optional<uint64_t> get_rng_state() const override {
    auto const cm =
        std::get_if<typename Kernels::StreamCollisionModelThermalized>(
            &*derived().m_collision_model);
    if (!cm or derived().m_kT == 0.) {
      return std::nullopt;
    }
    return {static_cast<uint64_t>(cm->getTime_step())};
  }

  void set_rng_state(uint64_t counter) override {
    auto const cm =
        std::get_if<typename Kernels::StreamCollisionModelThermalized>(
            &*derived().m_collision_model);
    if (!cm or derived().m_kT == 0.) {
      throw std::runtime_error("This LB instance is unthermalized");
    }
    assert(counter <=
           static_cast<uint32_t>(std::numeric_limits<uint_t>::max()));
    cm->setTime_step(static_cast<uint32_t>(counter));
  }

  // ---- Lattice position checker ----

  std::function<bool(Utils::Vector3d const &)>
  make_lattice_position_checker(bool consider_points_in_halo) const override {
    auto const &lat = *m_lattice;
    if (consider_points_in_halo) {
      return [&](Utils::Vector3d const &p) { return lat.pos_in_local_halo(p); };
    }
    return [&](Utils::Vector3d const &p) { return lat.pos_in_local_domain(p); };
  }

  // ---- Boundary handling (UBB scaffolding) ----

  std::optional<Utils::Vector3d>
  get_node_velocity_at_boundary(Utils::Vector3i const &node,
                                bool consider_ghosts = false) const override {
    assert(not(consider_ghosts and m_pending_ghost_comm.test(GhostComm::UBB)));
    auto const bc = get_block_and_cell(get_lattice(), node, consider_ghosts);
    if (!bc or !m_boundary[0]->node_is_boundary(node))
      return std::nullopt;

    return {to_vector3d(m_boundary[0]->get_node_value_at_boundary(node))};
  }

  bool set_node_velocity_at_boundary(Utils::Vector3i const &node,
                                     Utils::Vector3d const &velocity) override {
    derived().on_boundary_add();
    m_pending_ghost_comm.set(GhostComm::UBB);
    auto bc = get_block_and_cell(get_lattice(), node, false);
    if (!bc) {
      bc = get_block_and_cell(get_lattice(), node, true);
    }
    if (bc) {
      for (auto &boundary : m_boundary) {
        boundary->set_node_value_at_boundary(
            node, to_vector3<FloatType>(velocity), *bc);
      }
    }
    return bc.has_value();
  }

  std::vector<std::optional<Utils::Vector3d>> get_slice_velocity_at_boundary(
      Utils::Vector3i const &lower_corner,
      Utils::Vector3i const &upper_corner) const override {
    std::vector<std::optional<Utils::Vector3d>> out;
    for_each_block_in_slice(
        get_lattice(), lower_corner, upper_corner,
        [&](auto const &, auto const &bci, auto const &ci,
            auto const &block_offset) {
          if (out.empty())
            out.resize(ci.numCells());

          auto kernel = [&out, this](unsigned const, unsigned const local_index,
                                     Utils::Vector3i const &node) {
            if (m_boundary[0]->node_is_boundary(node)) {
              out[local_index] =
                  to_vector3d(m_boundary[0]->get_node_value_at_boundary(node));
            } else {
              out[local_index] = std::nullopt;
            }
          };

          copy_block_buffer(bci, ci, block_offset, lower_corner, kernel);
        });
    return out;
  }

  void set_slice_velocity_at_boundary(
      Utils::Vector3i const &lower_corner, Utils::Vector3i const &upper_corner,
      std::vector<std::optional<Utils::Vector3d>> const &velocity) override {
    derived().on_boundary_add();
    m_pending_ghost_comm.set(GhostComm::UBB);
    auto const &lattice = get_lattice();
    for_each_block_in_slice(
        lattice, lower_corner, upper_corner,
        [&](auto &block, auto const &bci, auto const &ci,
            auto const &block_offset) {
          assert(velocity.size() == ci.numCells());

          auto kernel = [&, this](unsigned const, unsigned const local_index,
                                  Utils::Vector3i const &node) {
            auto const bc = get_block_and_cell(lattice, node, false);
            assert(bc->block->getAABB() == block.getAABB());
            auto const &opt = velocity[local_index];
            if (opt) {
              for (auto &boundary : m_boundary) {
                boundary->set_node_value_at_boundary(
                    node, to_vector3<FloatType>(*opt), *bc);
              }
            } else {
              for (auto &boundary : m_boundary) {
                boundary->remove_node_from_boundary(node, *bc);
              }
            }
          };

          copy_block_buffer(bci, ci, block_offset, lower_corner, kernel);
        });
  }

  std::optional<std::vector<Utils::Vector3d>>
  get_node_boundary_force(Utils::Vector3i const &node) const override {
    auto const bc = get_block_and_cell(get_lattice(), node, true);
    if (!bc or !m_boundary[0]->node_is_boundary(node))
      return std::nullopt;

    return this->get_node_last_applied_force(node, true);
  }

  bool remove_node_from_boundary(Utils::Vector3i const &node) override {
    auto bc = get_block_and_cell(get_lattice(), node, false);
    if (!bc) {
      bc = get_block_and_cell(get_lattice(), node, true);
    }
    if (bc) {
      for (auto &boundary : m_boundary) {
        boundary->remove_node_from_boundary(node, *bc);
      }
    }
    return bc.has_value();
  }

  std::optional<bool>
  get_node_is_boundary(Utils::Vector3i const &node,
                       bool consider_ghosts = false) const override {
    assert(not(consider_ghosts and m_pending_ghost_comm.test(GhostComm::UBB)));
    auto const bc = get_block_and_cell(get_lattice(), node, consider_ghosts);
    if (!bc)
      return std::nullopt;

    return {m_boundary[0]->node_is_boundary(node)};
  }

  std::vector<bool>
  get_slice_is_boundary(Utils::Vector3i const &lower_corner,
                        Utils::Vector3i const &upper_corner) const override {
    std::vector<bool> out;
    for_each_block_in_slice(
        get_lattice(), lower_corner, upper_corner,
        [&](auto const &, auto const &bci, auto const &ci,
            auto const &block_offset) {
          if (out.empty())
            out.resize(ci.numCells());

          auto kernel = [&out, this](unsigned const, unsigned const local_index,
                                     Utils::Vector3i const &node) {
            out[local_index] = m_boundary[0]->node_is_boundary(node);
          };

          copy_block_buffer(bci, ci, block_offset, lower_corner, kernel);
        });
    return out;
  }

  void reallocate_ubb_field() override {
    for (auto &boundary : m_boundary) {
      boundary->boundary_update();
    }
  }

  void clear_boundaries() override {
    derived().reset_boundary_handling(get_lattice().get_blocks());
    m_pending_ghost_comm.set(GhostComm::UBB);
    this->ghost_communication();
    m_has_boundaries = false;
    derived().setup_streaming_communicator();
  }

  /**
   * @brief Set boundary conditions from a rasterized shape.
   * @param raster_flat  Flattened 3D mask (non-zero = boundary node).
   * @param data_flat    Flattened 3D array of slip velocities (3 components
   *                     per node, same ordering as @p raster_flat).
   */
  void
  update_boundary_from_shape(std::vector<int> const &raster_flat,
                             std::vector<double> const &data_flat) override {
    derived().on_boundary_add();
    m_pending_ghost_comm.set(GhostComm::UBB);
    auto const &grid_size = get_lattice().get_grid_dimensions();
    auto data = fill_3D_vector_array(data_flat, grid_size);
    for (auto &boundary : m_boundary) {
      set_boundary_from_grid(*boundary, get_lattice(), raster_flat, data);
    }
    this->ghost_communication();
    reallocate_ubb_field();
  }

  /** @brief Convert a flat (row-major) index to a 3D grid coordinate. */
  [[nodiscard]] Utils::Vector3i flat_index_to_node(int index) const {
    Utils::Vector3i node({0, 0, 0});
    auto const &grid_size = get_lattice().get_grid_dimensions();
    node[2] = index % grid_size[2];
    int tmp = index / grid_size[2];
    node[1] = tmp % grid_size[1];
    node[0] = tmp / grid_size[1];
    return node;
  }

  /**
   * @brief Get the fluid neighbor of a boundary node along stencil direction
   *        @p dir, with periodic wrapping.
   */
  [[nodiscard]] Utils::Vector3i get_neighbor_node(Utils::Vector3i const &node,
                                                  int dir) const {
    Utils::Vector3i neighbor({0, 0, 0});
    auto const &grid_size = get_lattice().get_grid_dimensions();
    auto constexpr neighbor_offset = BoundaryKernel::neighborOffset;
    for (int i = 0; i < neighbor.size(); i++) {
      neighbor[i] =
          (node[i] - neighbor_offset[i][dir] + grid_size[i]) % grid_size[i];
    }
    return neighbor;
  }

  /**
   * @brief Total force exerted by the fluid on a subset of boundary nodes.
   * Only boundary nodes where @p raster_flat is non-zero are included.
   */
  [[nodiscard]] Utils::Vector3d get_boundary_force_from_shape(
      std::vector<int> const &raster_flat) const override {
    Utils::Vector3d force({0, 0, 0});
    auto const &grid_size = get_lattice().get_grid_dimensions();
    for (auto &block : *get_lattice().get_blocks()) {
      auto const offset = get_lattice().get_block_corner(block, true);
      for (auto &boundary : m_boundary) {
        auto const &force_field = boundary->get_force_vector(&block);
        auto const &index_field = boundary->get_index_vector(&block);
        for (int i = 0; i < raster_flat.size(); i++) {
          if (raster_flat[i] != 0) {
            auto node = flat_index_to_node(i);
            if (get_lattice().node_in_local_halo(node)) {
              // shift node to local frame
              node = (node - offset + grid_size) % grid_size;
              for (int j = 0; j < index_field.size(); j++) {
                auto neighbor_node =
                    get_neighbor_node(node, index_field[j].dir);
                if (index_field[j].x == neighbor_node[0] &&
                    index_field[j].y == neighbor_node[1] &&
                    index_field[j].z == neighbor_node[2]) {
                  force[0] += force_field[j].F_0;
                  force[1] += force_field[j].F_1;
                  force[2] += force_field[j].F_2;
                }
              }
            }
          }
        }
      }
    }
    return zero_centered_to_md(force);
  }

  [[nodiscard]] Utils::Vector3d get_boundary_force() const override {
    Vector3<double> force(0.);
    for (auto &block : *get_lattice().get_blocks()) {
      for (auto &boundary : m_boundary) {
        force += boundary->get_total_force(&block);
      }
    }
    return zero_centered_to_md(to_vector3d(force));
  }

  // ---- Read-side interpolation (velocity + density at positions) ----

  auto make_velocity_interpolation_kernel() const {
    auto const &lattice = *m_lattice;
    auto const &blocks = *lattice.get_blocks();
    assert(lattice.get_ghost_layers() == 1u);
    return [&](Utils::Vector3d const &pos) {
      Utils::Vector3d acc{0., 0., 0.};
      interpolate_bspline_at_pos(pos, [&, field_id = m_velocity_field_id](
                                          std::array<int, 3> const node,
                                          double weight) {
        // Nodes with zero weight might not be accessible, because they can
        // be outside ghost layers
        if (weight != 0.) {
          auto block = get_block_extended(lattice, node, 1u);
          if (!block)
            throw interpolation_illegal_access("velocity", pos, node, weight);
          Vector3<FloatType> vel;
          if (m_has_boundaries and m_boundary[0]->node_is_boundary(node)) {
            vel = m_boundary[0]->get_node_value_at_boundary(node);
          } else {
            auto cell = to_cell(node);
            blocks.transformGlobalToBlockLocalCell(cell, *block);
            auto field =
                block->template uncheckedFastGetData<VectorField>(field_id);
            vel = lbm::accessor::Vector::get(field, cell);
          }
          acc += to_vector3d(vel) * weight;
        }
      });
      return acc;
    };
  }

  auto make_density_interpolation_kernel() const {
    auto const &lattice = *m_lattice;
    auto const &blocks = *lattice.get_blocks();
    assert(lattice.get_ghost_layers() == 1u);
    return [&, this](Utils::Vector3d const &pos) {
      double acc = 0.;
      interpolate_bspline_at_pos(pos, [&, density = m_density,
                                       field_id = derived().m_pdf_field_id[0]](
                                          std::array<int, 3> const node,
                                          double weight) {
        // Nodes with zero weight might not be accessible, because they can
        // be outside ghost layers
        if (weight != 0.) {
          auto block = get_block_extended(lattice, node, 1u);
          if (!block)
            throw interpolation_illegal_access("density", pos, node, weight);
          auto cell = to_cell(node);
          blocks.transformGlobalToBlockLocalCell(cell, *block);
          auto field = block->template uncheckedFastGetData<PdfField>(field_id);
          auto const rho = lbm::accessor::Density::get(field, density, cell);
          acc += rho * weight;
        }
      });
      return acc;
    };
  }

  std::optional<Utils::Vector3d>
  get_velocity_at_pos(Utils::Vector3d const &pos,
                      bool consider_points_in_halo = false) const override {
    assert(not m_pending_ghost_comm.test(GhostComm::VEL));
    assert(not m_pending_ghost_comm.test(GhostComm::UBB));
    if (!consider_points_in_halo and !m_lattice->pos_in_local_domain(pos))
      return std::nullopt;
    if (consider_points_in_halo and !m_lattice->pos_in_local_halo(pos))
      return std::nullopt;
    auto const kernel = make_velocity_interpolation_kernel();
    return {kernel(pos)};
  }

  std::optional<double>
  get_density_at_pos(Utils::Vector3d const &pos,
                     bool consider_points_in_halo = false) const override {
    assert(not m_pending_ghost_comm.test(GhostComm::PDF));
    if (!consider_points_in_halo and !m_lattice->pos_in_local_domain(pos))
      return std::nullopt;
    if (consider_points_in_halo and !m_lattice->pos_in_local_halo(pos))
      return std::nullopt;
    auto const kernel = make_density_interpolation_kernel();
    return {kernel(pos)};
  }

  std::vector<Utils::Vector3d>
  get_velocities_at_pos(std::vector<Utils::Vector3d> const &pos) override {
    if (pos.empty()) {
      return {};
    }
    std::vector<Utils::Vector3d> vel{};
    vel.reserve(pos.size());
    if constexpr (Architecture == lbmpy::Arch::CPU) {
      auto const kernel = make_velocity_interpolation_kernel();
      std::ranges::transform(pos, std::back_inserter(vel), kernel);
    }
#if defined(__CUDACC__) and defined(WALBERLA_BUILD_WITH_CUDA)
    if constexpr (Architecture == lbmpy::Arch::GPU) {
      auto const &lattice = get_lattice();
      auto const &block = *(lattice.get_blocks()->begin());
      auto const origin = block.getAABB().min();
      std::vector<FloatType> host_pos;
      host_pos.reserve(3ul * pos.size());
      assert(lattice.get_blocks()->getNumberOfBlocks() == 1u);
      for (auto const &vec : pos) {
#pragma unroll
        for (std::size_t i : {0ul, 1ul, 2ul}) {
          host_pos.emplace_back(static_cast<FloatType>(vec[i] - origin[i]));
        }
      }
      auto const gl = lattice.get_ghost_layers();
      auto field =
          block.template uncheckedFastGetData<VectorField>(m_velocity_field_id);
      // the velocity field has indeterminate values inside boundary regions;
      // we overwrite them with boundary slip velocities before interpolation
      auto const [dev_idx, dev_vel] = m_boundary[0]->get_flattened_map_device();
      if (not dev_idx->empty()) {
        lbm::accessor::Vector::set_from_list(field, *dev_idx, *dev_vel, gl);
      }
      auto const res =
          lbm::accessor::Interpolation::get_vel(field, host_pos, gl);
      for (auto it = res.begin(); it != res.end(); it += 3) {
        vel.emplace_back(Utils::Vector3d{static_cast<double>(*(it + 0)),
                                         static_cast<double>(*(it + 1)),
                                         static_cast<double>(*(it + 2))});
      }
    }
#endif
    return vel;
  }

  std::vector<double>
  get_densities_at_pos(std::vector<Utils::Vector3d> const &pos) override {
    if (pos.empty()) {
      return {};
    }
    std::vector<double> rho{};
    rho.reserve(pos.size());
    if constexpr (Architecture == lbmpy::Arch::CPU) {
      auto const kernel = make_density_interpolation_kernel();
      std::ranges::transform(pos, std::back_inserter(rho), kernel);
    }
#if defined(__CUDACC__) and defined(WALBERLA_BUILD_WITH_CUDA)
    if constexpr (Architecture == lbmpy::Arch::GPU) {
      auto const &lattice = get_lattice();
      auto const &block = *(lattice.get_blocks()->begin());
      auto const origin = block.getAABB().min();
      std::vector<FloatType> host_pos;
      host_pos.reserve(3ul * pos.size());
      assert(lattice.get_blocks()->getNumberOfBlocks() == 1u);
      for (auto const &vec : pos) {
#pragma unroll
        for (std::size_t i : {0ul, 1ul, 2ul}) {
          host_pos.emplace_back(static_cast<FloatType>(vec[i] - origin[i]));
        }
      }
      auto const gl = lattice.get_ghost_layers();
      auto field = block.template uncheckedFastGetData<PdfField>(
          derived().m_pdf_field_id[0]);
      auto res =
          lbm::accessor::Interpolation::get_rho(field, host_pos, m_density, gl);
      if constexpr (std::is_same_v<FloatType, double>) {
        std::swap(rho, res);
      } else {
        for (auto const &v : res) {
          rho.emplace_back(static_cast<double>(v));
        }
      }
    }
#endif
    return rho;
  }

protected:
  /**
   * @brief Scale data by a conversion factor (in-place).
   * Used for zero-centered density representation: LB internally stores
   * density fluctuations around zero, while the user interface uses
   * absolute densities. The conversion factors @ref m_zc_to_md and
   * @ref m_zc_to_lb translate between these representations.
   */
  template <typename T>
  void zero_centered_transform_impl(T &data, auto const factor) const {
    if constexpr (std::is_arithmetic_v<T>) {
      static_assert(std::is_floating_point_v<T>);
      data *= static_cast<T>(factor);
    } else {
      auto const coef = static_cast<typename T::value_type>(factor);
      std::transform(std::begin(data), std::end(data), std::begin(data),
                     [coef](auto value) { return value * coef; });
    }
  }

  void zero_centered_to_lb_in_place(auto &data) const {
    zero_centered_transform_impl(data, m_zc_to_lb);
  }

  void zero_centered_to_md_in_place(auto &data) const {
    zero_centered_transform_impl(data, m_zc_to_md);
  }

  auto zero_centered_to_lb(auto const &data) const {
    auto transformed_data = data;
    zero_centered_to_lb_in_place(transformed_data);
    return transformed_data;
  }

  auto zero_centered_to_md(auto const &data) const {
    auto transformed_data = data;
    zero_centered_to_md_in_place(transformed_data);
    return transformed_data;
  }

protected:
  // ---- Field allocation helper ----

  /**
   * @brief Convenience function to add a field with a custom allocator.
   *
   * When vectorization is off, let waLBerla decide which memory allocator
   * to use. When vectorization is on, the aligned memory allocator is
   * required, otherwise <tt>cpu_vectorize_info["assume_aligned"]</tt> will
   * trigger assertions. That is because for single-precision kernels the
   * waLBerla heuristic in <tt>src/field/allocation/FieldAllocator.h</tt>
   * will fall back to @c StdFieldAlloc, yet @c AllocateAligned is needed
   * for intrinsics to work.
   */
  template <typename Field> auto add_to_storage(std::string const tag) {
    auto const &blocks = m_lattice->get_blocks();
    auto const n_ghost_layers = m_lattice->get_ghost_layers();
    if constexpr (Architecture == lbmpy::Arch::CPU) {
#ifdef ESPRESSO_BUILD_WITH_AVX_KERNELS
      constexpr auto alignment = field::SIMDAlignment();
      using value_type = Field::value_type;
      using Allocator = field::AllocateAligned<value_type, alignment>;
      auto const allocator = std::make_shared<Allocator>();
      auto const empty_set = Set<SUID>::emptySet();
      return field::addToStorage<Field>(
          blocks, tag, field::internal::defaultSize, FloatType{0}, field::fzyx,
          n_ghost_layers, false, {}, empty_set, empty_set, allocator);
#else  // ESPRESSO_BUILD_WITH_AVX_KERNELS
      return field::addToStorage<Field>(blocks, tag, FloatType{0}, field::fzyx,
                                        n_ghost_layers);
#endif // ESPRESSO_BUILD_WITH_AVX_KERNELS
    }
#if defined(__CUDACC__) and defined(WALBERLA_BUILD_WITH_CUDA)
    else {
      auto field_id = gpu::addGPUFieldToStorage<GPUField>(
          blocks, tag, Field::F_SIZE, field::fzyx, n_ghost_layers);
      if constexpr (std::is_same_v<Field, _VectorField>) {
        for (auto &block : *blocks) {
          auto field = block.template getData<GPUField>(field_id);
          lbm::accessor::Vector::initialize(field, Vector3<FloatType>{0});
        }
      } else if constexpr (std::is_same_v<Field, _PdfField>) {
        for (auto &block : *blocks) {
          auto field = block.template getData<GPUField>(field_id);
          lbm::accessor::Population::initialize(
              field, std::array<FloatType, Stencil::Size>{});
        }
      } else if constexpr (std::is_same_v<Field, _ScalarField>) {
        for (auto &block : *blocks) {
          auto field = block.template getData<GPUField>(field_id);
          ek::accessor::Scalar::initialize(field, FloatType{0});
        }
      }
      return field_id;
    }
#endif
  }
};

} // namespace walberla
