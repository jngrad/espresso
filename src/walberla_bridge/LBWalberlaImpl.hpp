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
#ifndef LB_WALBERLA_H
#define LB_WALBERLA_H

/**
 * @file
 * @ref walberla::LBWalberlaImpl implements the interface of the LB
 * waLBerla bridge. It is a templated class that is specialized by lattice
 * models created by lbmpy (see <tt>maintainer/walberla_kernels</tt>).
 */

#include "blockforest/Initialization.h"
#include "blockforest/StructuredBlockForest.h"
#include "blockforest/communication/UniformBufferedScheme.h"
#include "boundary/BoundaryHandling.h"
#include "field/GhostLayerField.h"
#include "field/adaptors/GhostLayerFieldAdaptor.h"
#include "field/vtk/FlagFieldCellFilter.h"
#include "field/vtk/VTKWriter.h"
#include "lbm/lattice_model/CollisionModel.h"
#include "lbm/lattice_model/D3Q19.h"
#include "lbm/vtk/all.h"

#include "domain_decomposition/SharedSweep.h"

#include "field/AddToStorage.h"
#include "field/FlagField.h"
#include "field/adaptors/AdaptorCreators.h"
#include "field/communication/PackInfo.h"
#include "lbm/boundary/NoSlip.h"
#include "lbm/boundary/UBB.h"
#include "lbm/communication/PdfFieldPackInfo.h"
#include "lbm/field/Adaptors.h"
#include "lbm/field/AddToStorage.h"
#include "lbm/field/PdfField.h"
#include "lbm/lattice_model/D3Q19.h"
#include "lbm/sweeps/CellwiseSweep.h"

#include "stencil/D3Q19.h"
#include "stencil/D3Q27.h"

#include "LBWalberlaBase.hpp"
#include "ResetForce.hpp"
#include "generated_kernels/InitialPDFsSetter.h"
#include "generated_kernels/StreamSweep.h"
#include "generated_kernels/UBB.h"
#include "walberla_utils.hpp"

#ifdef __AVX2__
#include "generated_kernels/CollideSweepAVX.h"
#include "generated_kernels/CollideSweepThermalizedAVX.h"
#define ThermalizedCollisionModel                                              \
  walberla::pystencils::CollideSweepThermalizedAVX
#define UnthermalizedCollisionModel walberla::pystencils::CollideSweepAVX
#else
#include "generated_kernels/CollideSweep.h"
#include "generated_kernels/CollideSweepThermalized.h"
#define ThermalizedCollisionModel walberla::pystencils::CollideSweepThermalized
#define UnthermalizedCollisionModel walberla::pystencils::CollideSweep
#endif

#include <utils/Vector.hpp>
#include <utils/interpolation/bspline_3d.hpp>
#include <utils/math/make_lin_space.hpp>

#include <boost/optional.hpp>
#include <boost/tuple/tuple.hpp>
#include <boost/variant.hpp>

#include <array>
#include <cmath>
#include <cstddef>
#include <limits>
#include <map>
#include <memory>
#include <sstream>
#include <stdexcept>
#include <string>
#include <type_traits>
#include <utility>
#include <vector>

namespace walberla {
// forward declare
class LBWalberlaImpl;
} // namespace walberla

#include "generated_kernels/macroscopic_values_accessors.h"

namespace walberla {

// Flags marking fluid and boundaries
const FlagUID Fluid_flag("fluid");
const FlagUID UBB_flag("velocity bounce back");

/** Class that runs and controls the LB on WaLBerla
 */
class LBWalberlaImpl : public LBWalberlaBase {
protected:
  using CollisionModel =
      boost::variant<UnthermalizedCollisionModel, ThermalizedCollisionModel>;

private:
  auto generate_collide_sweep() const {
    real_t const omega = shear_mode_relaxation_rate();
    real_t const omega_odd = odd_mode_relaxation_rate(omega);
    std::shared_ptr<CollisionModel> ptr;
    if (m_kT == 0.) {
      auto obj = UnthermalizedCollisionModel(m_last_applied_force_field_id,
                                             m_pdf_field_id, omega, omega,
                                             omega_odd, omega);
      ptr = std::make_shared<CollisionModel>(std::move(obj));
    } else {
      auto obj = ThermalizedCollisionModel(
          m_last_applied_force_field_id, m_pdf_field_id, uint32_t(0u),
          uint32_t(0u), uint32_t(0u), m_kT, omega, omega, omega_odd, omega,
          m_seed, 0u);
      obj.block_offset_generator =
          [this](IBlock *const block, uint32_t &block_offset_0,
                 uint32_t &block_offset_1, uint32_t &block_offset_2) {
            block_offset_0 = m_blocks->getBlockCellBB(*block).xMin();
            block_offset_1 = m_blocks->getBlockCellBB(*block).yMin();
            block_offset_2 = m_blocks->getBlockCellBB(*block).zMin();
          };
      ptr = std::make_shared<CollisionModel>(std::move(obj));
    }
    return ptr;
  }

  class : public boost::static_visitor<> {
  public:
    template <typename CM> void operator()(CM &cm, IBlock *b) const { cm(b); }

  } run_collide_sweep;

  class : public boost::static_visitor<> {
  public:
    void operator()(UnthermalizedCollisionModel &cm, uint32_t) const {
      throw std::runtime_error("The LB does not use a random number generator");
    }

    void operator()(ThermalizedCollisionModel &cm, uint32_t time_step) const {
      cm.time_step_ = time_step;
    }

  } set_time_step_impl;

  class : public boost::static_visitor<uint32_t> {
  public:
    uint32_t operator()(UnthermalizedCollisionModel const &cm) const {
      throw std::runtime_error("The LB does not use a random number generator");
    }

    uint32_t operator()(ThermalizedCollisionModel const &cm) const {
      return cm.time_step_;
    }

  } get_time_step_impl;

  class : public boost::static_visitor<> {
  public:
    void operator()(UnthermalizedCollisionModel &cm) const {}

    void operator()(ThermalizedCollisionModel &cm) const { cm.time_step_++; }

  } increment_time_step_if_thermalized_impl;

  double shear_mode_relaxation_rate() const {
    return 2 / (6 * m_viscosity + 1);
  }

  double odd_mode_relaxation_rate(double shear_relaxation,
                                  double magic_number = 3. / 16.) const {
    return (4 - 2 * shear_relaxation) /
           (4 * magic_number * shear_relaxation + 2 - shear_relaxation);
  }

public:
  // Type definitions
  using LatticeModel_T = LBWalberlaImpl;
  typedef stencil::D3Q19 Stencil;
  using VectorField = GhostLayerField<real_t, 3u>;
  using FlagField = field::FlagField<uint8_t>;
  using PdfField = GhostLayerField<real_t, Stencil::Size>;

  static constexpr real_t w[19] = {
      1. / 3.,  1. / 18., 1. / 18., 1. / 18., 1. / 18., 1. / 18., 1. / 18.,
      1. / 36., 1. / 36., 1. / 36., 1. / 36., 1. / 36., 1. / 36., 1. / 36.,
      1. / 36., 1. / 36., 1. / 36., 1. / 36., 1. / 36.};
  static constexpr real_t wInv[19] = {3.,  18., 18., 18., 18., 18., 18.,
                                      36., 36., 36., 36., 36., 36., 36.,
                                      36., 36., 36., 36., 36.};
  static constexpr bool compressible = true;

protected:
  /** Velocity boundary condition */
  using UBB = lbm::espresso::UBB<LatticeModel_T, uint8_t, true, true>;
  using Boundaries = BoundaryHandling<FlagField, Stencil, UBB>;

private:
  // Backend classes can access private members:
  template <class LM, class Enable> friend class lbm::EquilibriumDistribution;
  template <class LM, class Enable> friend struct lbm::Equilibrium;
  template <class LM, class Enable>
  friend struct lbm::internal::AdaptVelocityToForce;
  template <class LM, class Enable> friend struct lbm::Density;
  template <class LM> friend struct lbm::DensityAndVelocity;
  template <class LM, class Enable>
  friend struct lbm::DensityAndMomentumDensity;
  template <class LM, class Enable> friend struct lbm::MomentumDensity;
  template <class LM, class It, class Enable>
  friend struct lbm::DensityAndVelocityRange;
  template <class LM, typename flag_t, bool A, bool B>
  friend class lbm::espresso::UBB;

  real_t getDensity(const BlockAndCell &bc) const {
    auto pdf_field = bc.block->template getData<PdfField>(m_pdf_field_id);
    return lbm::Density<LatticeModel_T>::get(*this, *pdf_field, bc.cell.x(),
                                             bc.cell.y(), bc.cell.z());
  }

  real_t getDensityAndVelocity(const BlockAndCell &bc,
                               Vector3<real_t> &velocity) const {
    auto const pdf_field = bc.block->template getData<PdfField>(m_pdf_field_id);
    auto const force_field =
        bc.block->template getData<VectorField>(m_last_applied_force_field_id);
    return getDensityAndVelocity(pdf_field, force_field, bc.cell.x(),
                                 bc.cell.y(), bc.cell.z(), velocity);
  }

  real_t getDensityAndVelocity(const PdfField *pdf_field,
                               const VectorField *force_field,
                               const cell_idx_t x, const cell_idx_t y,
                               const cell_idx_t z,
                               Vector3<real_t> &velocity) const {
    const real_t rho = lbm::DensityAndMomentumDensity<LatticeModel_T>::get(
        velocity, *force_field, *pdf_field, x, y, z);
    if constexpr (compressible) {
      const real_t invRho = real_t(1) / rho;
      velocity *= invRho;
    }
    return rho;
  }

  void setDensityAndVelocity(const BlockAndCell &bc,
                             Vector3<real_t> const &velocity, real_t rho) {
    auto pdf_field = bc.block->template getData<PdfField>(m_pdf_field_id);
    auto force_field =
        bc.block->template getData<VectorField>(m_last_applied_force_field_id);
    lbm::DensityAndVelocity<LatticeModel_T>::set(*pdf_field, bc.cell.x(),
                                                 bc.cell.y(), bc.cell.z(),
                                                 *force_field, velocity, rho);
  }

  Matrix3<real_t> getPressureTensor(const BlockAndCell &bc) const {
    Matrix3<real_t> pressureTensor;
    auto pdf_field = bc.block->template getData<PdfField>(m_pdf_field_id);
    lbm::PressureTensor<LatticeModel_T>::get(pressureTensor, *this, *pdf_field,
                                             bc.cell.x(), bc.cell.y(),
                                             bc.cell.z());
    return pressureTensor;
  }

protected:
  /** VTK writers that are executed automatically */
  std::map<std::string, std::pair<std::shared_ptr<vtk::VTKOutput>, bool>>
      m_vtk_auto;
  /** VTK writers that are executed manually */
  std::map<std::string, std::shared_ptr<vtk::VTKOutput>> m_vtk_manual;

  // Member variables
  Utils::Vector3i m_grid_dimensions;
  unsigned int m_n_ghost_layers;
  double m_viscosity;
  double m_density;
  double m_kT;
  unsigned int m_seed;

  // Block data access handles
  BlockDataID m_pdf_field_id;
  BlockDataID m_pdf_tmp_field_id;
  BlockDataID m_flag_field_id;

  BlockDataID m_last_applied_force_field_id;
  BlockDataID m_force_to_be_applied_id;

  BlockDataID m_velocity_field_id;

  BlockDataID m_boundary_handling_id;

  using FullCommunicator = blockforest::communication::UniformBufferedScheme<
      typename stencil::D3Q27>;
  std::shared_ptr<FullCommunicator> m_full_communication;
  using PDFStreamingCommunicator =
      blockforest::communication::UniformBufferedScheme<
          typename stencil::D3Q19>;
  std::shared_ptr<PDFStreamingCommunicator> m_pdf_streaming_communication;

  /** Block forest */
  std::shared_ptr<blockforest::StructuredBlockForest> m_blocks;

  // ResetForce sweep + external force handling
  std::shared_ptr<ResetForce<PdfField, VectorField>> m_reset_force;

  // Stream sweep
  std::shared_ptr<pystencils::StreamSweep> m_stream;

  // Collision sweep
  std::shared_ptr<CollisionModel> m_collision_model;

  // Boundary handling
  class LBBoundaryHandling {
  public:
    LBBoundaryHandling(const BlockDataID &flag_field_id,
                       const BlockDataID &pdf_field_id,
                       const BlockDataID &force_field_id,
                       const LatticeModel_T *lm)
        : m_flag_field_id(flag_field_id), m_pdf_field_id(pdf_field_id),
          m_force_field_id(force_field_id), m_lm(lm) {}

    Boundaries *operator()(IBlock *const block) {

      auto flag_field = block->template getData<FlagField>(m_flag_field_id);
      auto pdf_field = block->template getData<PdfField>(m_pdf_field_id);
      auto force_field = block->template getData<VectorField>(m_force_field_id);

      const auto fluid = flag_field->flagExists(Fluid_flag)
                             ? flag_field->getFlag(Fluid_flag)
                             : flag_field->registerFlag(Fluid_flag);

      return new Boundaries("boundary handling", flag_field, fluid,
                            UBB("velocity bounce back", UBB_flag, m_lm,
                                pdf_field, force_field, nullptr));
    }

  private:
    const BlockDataID m_flag_field_id;
    const BlockDataID m_pdf_field_id;
    const BlockDataID m_force_field_id;
    const LatticeModel_T *m_lm;
  };

  // Boundary sweep
  std::shared_ptr<UBB> m_boundary;

  class LeesEdwardsUpdate {
  public:
    LeesEdwardsUpdate(const std::shared_ptr<StructuredBlockForest> &blocks,
                      BlockDataID pdf_field_id, BlockDataID pdf_tmp_field_id,
                      LeesEdwardsCallbacks callbacks)
        : blocks_(blocks), m_pdf_field_id(pdf_field_id),
          m_pdf_tmp_field_id(pdf_tmp_field_id), m_callbacks(callbacks) {}

    void operator()(IBlock *block) {
      // TODO should dimension_x contain the ghost layers or not. At the moment
      // value is 64 with GL it is 66. In the lbmpy Leed Edwards this is mixed.
      // Probably not good
      auto const offset = m_callbacks->get_pos_offset();

      // Top cells
      if (blocks_->atDomainYMaxBorder(*block)) {
        uint_t dimension_x = blocks_->getNumberOfXCells(*block);
        real_t weight = fmod(offset + real_c(dimension_x), 1.0);

        auto pdf_field = block->template getData<PdfField>(m_pdf_field_id);
        auto pdf_tmp_field =
            block->template getData<PdfField>(m_pdf_tmp_field_id);

        CellInterval ci;
        pdf_field->getGhostRegion(stencil::N, ci, 1, true);

        for (auto cell = ci.begin(); cell != ci.end(); ++cell) {
          cell_idx_t x = cell->x();

          uint_t ind1 = uint_c(floor(x - offset)) % dimension_x;
          uint_t ind2 = uint_c(ceil(x - offset)) % dimension_x;

          for (uint_t q = 0; q < Stencil::Q; ++q) {
            pdf_tmp_field->get(*cell, 0) =
                (1 - weight) *
                    pdf_field->get(cell_idx_c(ind1), cell->y(), cell->z(), q) +
                weight *
                    pdf_field->get(cell_idx_c(ind2), cell->y(), cell->z(), q);
            // printf("%f -> %f\n", pdf_field->get(*cell, 0), pdf_tmp_field->get(*cell, 0));
          }
        }
      }
      // Bottom cells
      if (blocks_->atDomainYMinBorder(*block)) {
        uint_t dimension_x = blocks_->getNumberOfXCells(*block);
        real_t weight = fmod(offset + real_c(dimension_x), 1.0);

        auto pdf_field = block->template getData<PdfField>(m_pdf_field_id);
        auto pdf_tmp_field =
            block->template getData<PdfField>(m_pdf_tmp_field_id);

        CellInterval ci;
        pdf_field->getGhostRegion(stencil::S, ci, 1, true);

        for (auto cell = ci.begin(); cell != ci.end(); ++cell) {
          cell_idx_t x = cell->x();

          uint_t ind1 = uint_c(floor(x + offset)) % dimension_x;
          uint_t ind2 = uint_c(ceil(x + offset)) % dimension_x;

          for (uint_t q = 0; q < Stencil::Q; ++q) {
            pdf_tmp_field->get(*cell, 0) =
                (1 - weight) *
                    pdf_field->get(cell_idx_c(ind1), cell->y(), cell->z(), q) +
                weight *
                    pdf_field->get(cell_idx_c(ind2), cell->y(), cell->z(), q);
          }
        }
      }
    }

  private:
    const std::shared_ptr<StructuredBlockForest> &blocks_;
    BlockDataID m_pdf_field_id;
    BlockDataID m_pdf_tmp_field_id;
    boost::optional<LeesEdwardsCallbacks> m_callbacks;
  };

  /** Sweep that swaps @c pdf_field and @c pdf_tmp_field.
   */
  class LeesEdwardsSwap {
  public:
    LeesEdwardsSwap(BlockDataID pdf_field_id, BlockDataID pdf_tmp_field_id)
        : m_pdf_field_id(pdf_field_id), m_pdf_tmp_field_id(pdf_tmp_field_id) {}

    void operator()(IBlock *block) {
      auto pdf_field = block->template getData<PdfField>(m_pdf_field_id);
      auto m_pdf_tmp_field =
          block->template getData<PdfField>(m_pdf_tmp_field_id);
      m_pdf_tmp_field->swapDataPointers(pdf_field);
    }

  private:
    const BlockDataID m_pdf_field_id;
    const BlockDataID m_pdf_tmp_field_id;
  };

  // Lees-Edwards sweep
  std::shared_ptr<LeesEdwardsUpdate> m_lees_edwards_update_sweep;
  std::shared_ptr<LeesEdwardsSwap> m_lees_edwards_swap_sweep;
  boost::optional<LeesEdwardsCallbacks> m_lees_edwards_callbacks;

  template <typename LatticeModel_T, typename OutputType = float>
  class DensityVTKWriter : public vtk::BlockCellDataWriter<OutputType> {
  public:
    using PdfField_T = typename LatticeModel_T::PdfField;

    DensityVTKWriter(const ConstBlockDataID &pdf, const std::string &id)
        : vtk::BlockCellDataWriter<OutputType>(id), bdid_(pdf), pdf_(nullptr) {}

  protected:
    void configure() override {
      WALBERLA_ASSERT_NOT_NULLPTR(this->block_);
      pdf_ = this->block_->template getData<PdfField_T>(bdid_);
    }

    OutputType evaluate(const cell_idx_t x, const cell_idx_t y,
                        const cell_idx_t z, const cell_idx_t /*f*/) override {
      WALBERLA_ASSERT_NOT_NULLPTR(pdf_);
      return numeric_cast<OutputType>(
          lbm::Density<LatticeModel_T>::get(*pdf_, x, y, z));
    }

    const ConstBlockDataID bdid_;
    const PdfField_T *pdf_;
  };

  std::size_t stencil_size() const override {
    return static_cast<std::size_t>(Stencil::Size);
  }

public:
  LBWalberlaImpl(double viscosity, double density,
                 const Utils::Vector3i &grid_dimensions,
                 const Utils::Vector3i &node_grid, unsigned int n_ghost_layers,
                 double kT, unsigned int seed)
      : m_grid_dimensions(grid_dimensions), m_n_ghost_layers(n_ghost_layers),
        m_viscosity(viscosity), m_density(density), m_kT(kT), m_seed(seed) {

    if (n_ghost_layers == 0)
      throw std::runtime_error("At least one ghost layer must be used");
    for (int i : {0, 1, 2}) {
      if (m_grid_dimensions[i] % node_grid[i] != 0) {
        throw std::runtime_error(
            "LB grid dimensions and mpi node grid are not compatible.");
      }
    }

    m_blocks = blockforest::createUniformBlockGrid(
        uint_c(node_grid[0]), // blocks in x direction
        uint_c(node_grid[1]), // blocks in y direction
        uint_c(node_grid[2]), // blocks in z direction
        uint_c(m_grid_dimensions[0] /
               node_grid[0]), // number of cells per block in x direction
        uint_c(m_grid_dimensions[1] /
               node_grid[1]), // number of cells per block in y direction
        uint_c(m_grid_dimensions[2] /
               node_grid[2]), // number of cells per block in z direction
        1,                    // Lattice constant
        uint_c(node_grid[0]), uint_c(node_grid[1]),
        uint_c(node_grid[2]), // cpus per direction
        true, true, true);

    // Init and register fields
    m_pdf_field_id = field::addToStorage<PdfField>(
        m_blocks, "pdfs", real_t{0}, field::fzyx, m_n_ghost_layers);
    m_pdf_tmp_field_id = field::addToStorage<PdfField>(
        m_blocks, "pdfs_tmp", real_t{0}, field::fzyx, m_n_ghost_layers);
    m_last_applied_force_field_id = field::addToStorage<VectorField>(
        m_blocks, "force field", real_t{0}, field::fzyx, m_n_ghost_layers);
    m_force_to_be_applied_id = field::addToStorage<VectorField>(
        m_blocks, "force field", real_t{0}, field::fzyx, m_n_ghost_layers);
    m_velocity_field_id = field::addToStorage<VectorField>(
        m_blocks, "velocity field", real_t{0}, field::fzyx, m_n_ghost_layers);

    // Init and register flag field (fluid/boundary)
    m_flag_field_id = field::addFlagFieldToStorage<FlagField>(
        m_blocks, "flag field", m_n_ghost_layers);

    // Init and register pdf field
    auto pdf_setter =
        pystencils::InitialPDFsSetter(m_force_to_be_applied_id, m_pdf_field_id,
                                      m_velocity_field_id, real_c(m_density));
    for (auto b = m_blocks->begin(); b != m_blocks->end(); ++b) {
      pdf_setter(&(*b));
    }

    // Register boundary handling
    m_boundary_handling_id = m_blocks->addBlockData<Boundaries>(
        LBBoundaryHandling(m_flag_field_id, m_pdf_field_id,
                           m_last_applied_force_field_id, this),
        "boundary handling");
    clear_boundaries();

    // sets up the communication and registers pdf field and force field to it
    m_pdf_streaming_communication =
        std::make_shared<PDFStreamingCommunicator>(m_blocks);
    m_pdf_streaming_communication->addPackInfo(
        std::make_shared<field::communication::PackInfo<PdfField>>(
            m_pdf_field_id, m_n_ghost_layers));
    m_pdf_streaming_communication->addPackInfo(
        std::make_shared<field::communication::PackInfo<VectorField>>(
            m_last_applied_force_field_id, m_n_ghost_layers));

    m_full_communication = std::make_shared<FullCommunicator>(m_blocks);
    m_full_communication->addPackInfo(
        std::make_shared<field::communication::PackInfo<PdfField>>(
            m_pdf_field_id, m_n_ghost_layers));
    m_full_communication->addPackInfo(
        std::make_shared<field::communication::PackInfo<VectorField>>(
            m_last_applied_force_field_id, m_n_ghost_layers));
    m_full_communication->addPackInfo(
        std::make_shared<field::communication::PackInfo<VectorField>>(
            m_velocity_field_id, m_n_ghost_layers));

    // Instance the sweep responsible for force double buffering and
    // external forces
    m_reset_force = std::make_shared<ResetForce<PdfField, VectorField>>(
        m_pdf_field_id, m_last_applied_force_field_id,
        m_force_to_be_applied_id);

    // Prepare LB sweeps
    // Note: For now, combined collide-stream sweeps cannot be used,
    // because the collide-push variant is not supported by lbmpy.
    // The following functors are individual in-place collide and stream steps
    m_stream = std::make_shared<pystencils::StreamSweep>(
        m_last_applied_force_field_id, m_pdf_field_id, m_velocity_field_id);
    m_collision_model = generate_collide_sweep();

    // Synchronize ghost layers
    (*m_full_communication)();
  }

  void add_lees_edwards(LeesEdwardsCallbacks &&lees_edwards_callbacks) {
    m_lees_edwards_callbacks = std::move(lees_edwards_callbacks);
    m_lees_edwards_update_sweep = std::make_shared<LeesEdwardsUpdate>(
        m_blocks, m_pdf_field_id, m_pdf_tmp_field_id,
        *m_lees_edwards_callbacks);
    m_lees_edwards_swap_sweep =
        std::make_shared<LeesEdwardsSwap>(m_pdf_field_id, m_pdf_tmp_field_id);
  }

  void integrate() override {
    // Reset force fields
    for (auto b = m_blocks->begin(); b != m_blocks->end(); ++b)
      (*m_reset_force)(&*b);
    // LB collide
    for (auto b = m_blocks->begin(); b != m_blocks->end(); ++b)
      boost::apply_visitor(run_collide_sweep, *m_collision_model,
                           boost::variant<IBlock *>(&*b));
    (*m_pdf_streaming_communication)();

    // Lees-Edwards shift
    if (m_lees_edwards_callbacks) {
      for (auto b = m_blocks->begin(); b != m_blocks->end(); ++b)
        (*m_lees_edwards_update_sweep)(&*b);
      for (auto b = m_blocks->begin(); b != m_blocks->end(); ++b)
        (*m_lees_edwards_swap_sweep)(&*b);
    }
    // Handle boundaries
    for (auto b = m_blocks->begin(); b != m_blocks->end(); ++b)
      Boundaries::getBlockSweep(m_boundary_handling_id)(&*b);
    // LB stream
    for (auto b = m_blocks->begin(); b != m_blocks->end(); ++b)
      (*m_stream)(&*b);
    (*m_full_communication)();

    // Handle VTK writers
    for (auto it = m_vtk_auto.begin(); it != m_vtk_auto.end(); ++it) {
      if (it->second.second)
        vtk::writeFiles(it->second.first)();
    }

    boost::apply_visitor(increment_time_step_if_thermalized_impl,
                         *m_collision_model);
  }

  void ghost_communication() override { (*m_full_communication)(); }

  void set_viscosity(double viscosity) override { m_viscosity = viscosity; }

  double get_viscosity() const override { return m_viscosity; }

  // Velocity
  boost::optional<Utils::Vector3d>
  get_node_velocity(const Utils::Vector3i &node,
                    bool consider_ghosts = false) const override {
    boost::optional<bool> is_boundary =
        get_node_is_boundary(node, consider_ghosts);
    if (is_boundary)    // is info available locally
      if (*is_boundary) // is the node a boundary
        return get_node_velocity_at_boundary(node);
    auto const bc =
        get_block_and_cell(node, consider_ghosts, m_blocks, n_ghost_layers());
    if (!bc)
      return {};
    auto const &vel_field =
        (*bc).block->template getData<VectorField>(m_velocity_field_id);
    return Utils::Vector3d{double_c(vel_field->get((*bc).cell, uint_t(0u))),
                           double_c(vel_field->get((*bc).cell, uint_t(1u))),
                           double_c(vel_field->get((*bc).cell, uint_t(2u)))};
  }
  bool set_node_velocity(const Utils::Vector3i &node,
                         const Utils::Vector3d &v) override {
    auto bc = get_block_and_cell(node, false, m_blocks, n_ghost_layers());
    if (!bc)
      return false;
    // We have to set both, the pdf and the stored velocity field
    auto const density = getDensity(*bc);
    auto const vel = Vector3<real_t>{real_c(v[0]), real_c(v[1]), real_c(v[2])};
    setDensityAndVelocity(*bc, vel, density);
    auto vel_field =
        (*bc).block->template getData<VectorField>(m_velocity_field_id);
    for (uint_t f = 0u; f < 3u; ++f) {
      vel_field->get((*bc).cell, f) = real_c(v[f]);
    }

    return true;
  }
  boost::optional<Utils::Vector3d>
  get_velocity_at_pos(const Utils::Vector3d &pos,
                      bool consider_points_in_halo = false) const override {
    if (!consider_points_in_halo and !pos_in_local_domain(pos))
      return {};
    if (consider_points_in_halo and !pos_in_local_halo(pos))
      return {};
    Utils::Vector3d v{0.0, 0.0, 0.0};
    interpolate_bspline_at_pos(
        pos, [this, &v, pos](const std::array<int, 3> node, double weight) {
          // Nodes with zero weight might not be accessible, because they can be
          // outside ghost layers
          if (weight != 0) {
            auto res = get_node_velocity(
                Utils::Vector3i{{node[0], node[1], node[2]}}, true);
            if (!res) {
              printf("Pos: %g %g %g, Node %d %d %d, weight %g\n", pos[0],
                     pos[1], pos[2], node[0], node[1], node[2], weight);
              throw std::runtime_error("Access to LB velocity field failed.");
            }
            v += *res * weight;
          }
        });
    return {v};
  }

  boost::optional<double> get_interpolated_density_at_pos(
      const Utils::Vector3d &pos,
      bool consider_points_in_halo = false) const override {
    if (!consider_points_in_halo and !pos_in_local_domain(pos))
      return {};
    if (consider_points_in_halo and !pos_in_local_halo(pos))
      return {};
    double dens = 0.0;
    interpolate_bspline_at_pos(
        pos, [this, &dens, pos](const std::array<int, 3> node, double weight) {
          // Nodes with zero weight might not be accessible, because they can be
          // outside ghost layers
          if (weight != 0) {
            auto res =
                get_node_density(Utils::Vector3i{{node[0], node[1], node[2]}});
            if (!res) {
              printf("Pos: %g %g %g, Node %d %d %d, weight %g\n", pos[0],
                     pos[1], pos[2], node[0], node[1], node[2], weight);
              throw std::runtime_error("Access to LB density field failed.");
            }
            dens += *res * weight;
          }
        });
    return {dens};
  }

  // Local force
  bool add_force_at_pos(const Utils::Vector3d &pos,
                        const Utils::Vector3d &force) override {
    if (!pos_in_local_halo(pos))
      return false;
    auto force_at_node = [this, force](const std::array<int, 3> node,
                                       double weight) {
      auto const bc = get_block_and_cell(to_vector3i(node), true, m_blocks,
                                         n_ghost_layers());
      if (bc) {
        auto force_field = (*bc).block->template getData<VectorField>(
            m_force_to_be_applied_id);
        for (int i : {0, 1, 2})
          force_field->get((*bc).cell, i) += real_c(force[i] * weight);
      }
    };
    interpolate_bspline_at_pos(pos, force_at_node);
    return true;
  }

  boost::optional<Utils::Vector3d>
  get_node_force_to_be_applied(const Utils::Vector3i &node) const override {
    auto const bc = get_block_and_cell(node, true, m_blocks, n_ghost_layers());
    if (!bc)
      return {};

    auto const &force_field =
        (*bc).block->template getData<VectorField>(m_force_to_be_applied_id);
    return Utils::Vector3d{double_c(force_field->get((*bc).cell, uint_t(0u))),
                           double_c(force_field->get((*bc).cell, uint_t(1u))),
                           double_c(force_field->get((*bc).cell, uint_t(2u)))};
  }

  bool set_node_last_applied_force(Utils::Vector3i const &node,
                                   Utils::Vector3d const &force) override {
    auto bc = get_block_and_cell(node, false, m_blocks, n_ghost_layers());
    if (!bc)
      return false;

    auto force_field = (*bc).block->template getData<VectorField>(
        m_last_applied_force_field_id);
    for (uint_t f = 0u; f < 3u; ++f) {
      force_field->get((*bc).cell, f) = real_c(force[f]);
    }

    return true;
  }

  boost::optional<Utils::Vector3d>
  get_node_last_applied_force(const Utils::Vector3i &node,
                              bool consider_ghosts = false) const override {
    auto const bc =
        get_block_and_cell(node, consider_ghosts, m_blocks, n_ghost_layers());
    if (!bc)
      return {};

    auto const force_field = (*bc).block->template getData<VectorField>(
        m_last_applied_force_field_id);
    return Utils::Vector3d{double_c(force_field->get((*bc).cell, uint_t(0u))),
                           double_c(force_field->get((*bc).cell, uint_t(1u))),
                           double_c(force_field->get((*bc).cell, uint_t(2u)))};
  }

  // Population
  bool set_node_pop(const Utils::Vector3i &node,
                    std::vector<double> const &population) override {
    auto bc = get_block_and_cell(node, false, m_blocks, n_ghost_layers());
    if (!bc)
      return false;

    auto pdf_field = (*bc).block->template getData<PdfField>(m_pdf_field_id);
    assert(population.size() == Stencil::Size);
    for (uint_t f = 0u; f < Stencil::Size; ++f) {
      pdf_field->get((*bc).cell, f) = real_c(population[f]);
    }

    return true;
  }

  boost::optional<std::vector<double>>
  get_node_pop(const Utils::Vector3i &node) const override {
    auto bc = get_block_and_cell(node, false, m_blocks, n_ghost_layers());
    if (!bc)
      return {boost::none};

    auto pdf_field = bc->block->template getData<PdfField>(m_pdf_field_id);
    std::vector<double> population(Stencil::Size);
    for (uint_t f = 0u; f < Stencil::Size; ++f) {
      population[f] = double_c(pdf_field->get((*bc).cell, f));
    }

    return {population};
  }

  // Density
  bool set_node_density(const Utils::Vector3i &node, double density) override {
    auto bc = get_block_and_cell(node, false, m_blocks, n_ghost_layers());
    if (!bc)
      return false;

    Vector3<real_t> vel;
    getDensityAndVelocity(*bc, vel);
    setDensityAndVelocity(*bc, vel, density);

    return true;
  }

  boost::optional<double>
  get_node_density(const Utils::Vector3i &node) const override {
    auto bc = get_block_and_cell(node, false, m_blocks, n_ghost_layers());
    if (!bc)
      return {boost::none};

    auto const density = getDensity(*bc);
    return {double_c(density)};
  }

  // Boundary related
  boost::optional<Utils::Vector3d>
  get_node_velocity_at_boundary(const Utils::Vector3i &node) const override {
    auto bc = get_block_and_cell(node, true, m_blocks, n_ghost_layers());
    if (!bc)
      return {boost::none};
    const Boundaries *boundary_handling =
        (*bc).block->template getData<Boundaries>(m_boundary_handling_id);
    boundary::BoundaryUID uid = boundary_handling->getBoundaryUID(UBB_flag);

    if (!boundary_handling->isBoundary((*bc).cell))
      return {boost::none};

    return {to_vector3d(
        boundary_handling->template getBoundaryCondition<UBB>(uid).getValue(
            (*bc).cell[0], (*bc).cell[1], (*bc).cell[2]))};
  }

  bool set_node_velocity_at_boundary(const Utils::Vector3i &node,
                                     const Utils::Vector3d &v) override {
    auto bc = get_block_and_cell(node, true, m_blocks, n_ghost_layers());
    if (!bc)
      return false;

    const typename UBB::Velocity velocity(real_c(v[0]), real_c(v[1]),
                                          real_c(v[2]));

    auto *boundary_handling =
        (*bc).block->template getData<Boundaries>(m_boundary_handling_id);
    boundary_handling->forceBoundary(UBB_flag, bc->cell[0], bc->cell[1],
                                     bc->cell[2], velocity);
    return true;
  }

  boost::optional<Utils::Vector3d>
  get_node_boundary_force(const Utils::Vector3i &node) const override {
    auto bc = get_block_and_cell(node, true, m_blocks,
                                 n_ghost_layers()); // including ghosts
    if (!bc)
      return {boost::none};
    // Get boundary handling
    auto const &bh =
        (*bc).block->template getData<Boundaries>(m_boundary_handling_id);
    auto const &ff = (*bc).block->template getData<FlagField>(m_flag_field_id);
    try {
      if (!ff->isFlagSet((*bc).cell, ff->getFlag(UBB_flag)))
        return {boost::none};
    } catch (std::exception &e) {
      return {boost::none};
    }

    auto const uid = bh->getBoundaryUID(UBB_flag);
    auto const &ubb = bh->template getBoundaryCondition<UBB>(uid);
    return {to_vector3d(
        ubb.getForce((*bc).cell.x(), (*bc).cell.y(), (*bc).cell.z()))};
  }

  bool remove_node_from_boundary(const Utils::Vector3i &node) override {
    auto bc = get_block_and_cell(node, true, m_blocks, n_ghost_layers());
    if (!bc)
      return false;
    auto *boundary_handling =
        (*bc).block->template getData<Boundaries>(m_boundary_handling_id);
    boundary_handling->removeBoundary((*bc).cell[0], (*bc).cell[1],
                                      (*bc).cell[2]);
    return true;
  }

  boost::optional<bool>
  get_node_is_boundary(const Utils::Vector3i &node,
                       bool consider_ghosts = false) const override {
    auto bc =
        get_block_and_cell(node, consider_ghosts, m_blocks, n_ghost_layers());
    if (!bc)
      return {boost::none};

    auto *boundary_handling =
        (*bc).block->template getData<Boundaries>(m_boundary_handling_id);
    return {boundary_handling->isBoundary((*bc).cell)};
  }

  void clear_boundaries() override {
    const CellInterval &domain_bb_in_global_cell_coordinates =
        m_blocks->getCellBBFromAABB(
            m_blocks->begin()->getAABB().getExtended(real_c(n_ghost_layers())));
    for (auto block = m_blocks->begin(); block != m_blocks->end(); ++block) {

      auto *boundary_handling =
          block->template getData<Boundaries>(m_boundary_handling_id);

      CellInterval domain_bb(domain_bb_in_global_cell_coordinates);
      m_blocks->transformGlobalToBlockLocalCellInterval(domain_bb, *block);

      boundary_handling->fillWithDomain(domain_bb);
    }
  }

  // Pressure tensor
  boost::optional<Utils::Vector6d>
  get_node_pressure_tensor(const Utils::Vector3i &node) const override {
    auto bc = get_block_and_cell(node, false, m_blocks, n_ghost_layers());
    if (!bc)
      return {boost::none};
    return to_vector6d(getPressureTensor(*bc));
  }

  // Global momentum
  Utils::Vector3d get_momentum() const override {
    Vector3<real_t> mom;
    for (auto block_it = m_blocks->begin(); block_it != m_blocks->end();
         ++block_it) {
      auto pdf_field = block_it->template getData<PdfField>(m_pdf_field_id);
      auto force_field = block_it->template getData<VectorField>(
          m_last_applied_force_field_id);
      Vector3<real_t> local_v;
      WALBERLA_FOR_ALL_CELLS_XYZ(pdf_field, {
        real_t local_dens =
            getDensityAndVelocity(pdf_field, force_field, x, y, z, local_v);
        mom += local_dens * local_v;
      });
    }
    return to_vector3d(mom);
  }
  // Global external force
  void set_external_force(const Utils::Vector3d &ext_force) override {
    m_reset_force->set_ext_force(ext_force);
  }
  Utils::Vector3d get_external_force() const override {
    return m_reset_force->get_ext_force();
  }

  double get_kT() const override { return m_kT; }

  uint64_t get_rng_state() const override {
    return boost::apply_visitor(get_time_step_impl, *m_collision_model);
  }
  void set_rng_state(uint64_t counter) override {
    boost::apply_visitor(set_time_step_impl, *m_collision_model,
                         boost::variant<uint64_t>(counter));
  }

  // Grid, domain, halo
  int n_ghost_layers() const override {
    return static_cast<int>(m_n_ghost_layers);
  }
  Utils::Vector3i get_grid_dimensions() const override {
    return m_grid_dimensions;
  }
  std::pair<Utils::Vector3d, Utils::Vector3d>
  get_local_domain() const override {
    // We only have one block per mpi rank
    assert(++(m_blocks->begin()) == m_blocks->end());

    auto const ab = m_blocks->begin()->getAABB();
    return {to_vector3d(ab.min()), to_vector3d(ab.max())};
  }

  bool node_in_local_domain(const Utils::Vector3i &node) const override {
    // Note: Lattice constant =1, cell centers offset by .5
    return get_block_and_cell(node, false, m_blocks, n_ghost_layers()) !=
           boost::none;
  }
  bool node_in_local_halo(const Utils::Vector3i &node) const override {
    return get_block_and_cell(node, true, m_blocks, n_ghost_layers()) !=
           boost::none;
  }
  bool pos_in_local_domain(const Utils::Vector3d &pos) const override {
    return get_block(pos, false, m_blocks, n_ghost_layers()) != nullptr;
  }
  bool pos_in_local_halo(const Utils::Vector3d &pos) const override {
    return get_block(pos, true, m_blocks, n_ghost_layers()) != nullptr;
  }

  std::vector<std::pair<Utils::Vector3i, Utils::Vector3d>>
  node_indices_positions(bool include_ghosts = false) const override {
    int ghost_offset = 0;
    if (include_ghosts)
      ghost_offset = static_cast<int>(m_n_ghost_layers);
    std::vector<std::pair<Utils::Vector3i, Utils::Vector3d>> res;
    for (auto block = m_blocks->begin(); block != m_blocks->end(); ++block) {
      auto left = block->getAABB().min();
      // Lattice constant is 1, node centers are offset by .5
      Utils::Vector3d pos_offset =
          to_vector3d(left) + Utils::Vector3d::broadcast(.5);

      // Lattice constant is 1, so cast left corner position to ints
      Utils::Vector3i index_offset =
          Utils::Vector3i{int(left[0]), int(left[1]), int(left[2])};

      // Get field data which knows about the indices
      // In the loop, x,y,z are in block-local coordinates
      auto pdf_field = block->template getData<PdfField>(m_pdf_field_id);
      for (int x = -ghost_offset; x < int(pdf_field->xSize()) + ghost_offset;
           x++) {
        for (int y = -ghost_offset; y < int(pdf_field->ySize()) + ghost_offset;
             y++) {
          for (int z = -ghost_offset;
               z < int(pdf_field->zSize()) + ghost_offset; z++) {
            res.push_back({index_offset + Utils::Vector3i{x, y, z},
                           pos_offset + Utils::Vector3d{double(x), double(y),
                                                        double(z)}});
          }
        }
      }
    }
    return res;
  }

  void create_vtk(unsigned delta_N, unsigned initial_count,
                  unsigned flag_observables, std::string const &identifier,
                  std::string const &base_folder,
                  std::string const &prefix) override {
    // VTKOuput object must be unique
    std::stringstream unique_identifier;
    unique_identifier << base_folder << "/" << identifier;
    std::string const vtk_uid = unique_identifier.str();
    if (m_vtk_auto.find(vtk_uid) != m_vtk_auto.end() or
        m_vtk_manual.find(vtk_uid) != m_vtk_manual.end()) {
      throw std::runtime_error("VTKOutput object " + vtk_uid +
                               " already exists");
    }

    // instantiate VTKOutput object
    unsigned const write_freq = (delta_N) ? static_cast<unsigned>(delta_N) : 1u;
    auto pdf_field_vtk = vtk::createVTKOutput_BlockData(
        m_blocks, identifier, uint_c(write_freq), uint_c(0), false, base_folder,
        prefix, true, true, true, true, uint_c(initial_count));
    field::FlagFieldCellFilter<FlagField> fluid_filter(m_flag_field_id);
    fluid_filter.addFlag(Fluid_flag);
    pdf_field_vtk->addCellInclusionFilter(fluid_filter);

    // add writers
    if (static_cast<unsigned>(OutputVTK::density) & flag_observables) {
      pdf_field_vtk->addCellDataWriter(
          make_shared<DensityVTKWriter<LatticeModel_T, float>>(
              m_pdf_field_id, "DensityFromPDF"));
    }
    if (static_cast<unsigned>(OutputVTK::velocity_vector) & flag_observables) {
      pdf_field_vtk->addCellDataWriter(
          make_shared<field::VTKWriter<VectorField, float>>(
              m_velocity_field_id, "VelocityFromVelocityField"));
    }
    /* TODO WALBERLA: pressure tensor
    if (static_cast<unsigned>(OutputVTK::pressure_tensor) & flag_observables) {
      pdf_field_vtk->addCellDataWriter(
          make_shared<lbm::PressureTensorVTKWriter<LatticeModel_T, float>>(
              m_pdf_field_id, "PressureTensorFromPDF"));
    }
    */

    // register object
    if (delta_N) {
      m_vtk_auto[vtk_uid] = {pdf_field_vtk, true};
    } else {
      m_vtk_manual[vtk_uid] = pdf_field_vtk;
    }
  }

  /** Manually call a VTK callback */
  void write_vtk(std::string const &vtk_uid) override {
    if (m_vtk_auto.find(vtk_uid) != m_vtk_auto.end()) {
      throw std::runtime_error("VTKOutput object " + vtk_uid +
                               " is an automatic observable");
    }
    if (m_vtk_manual.find(vtk_uid) == m_vtk_manual.end()) {
      throw std::runtime_error("VTKOutput object " + vtk_uid +
                               " doesn't exist");
    }
    vtk::writeFiles(m_vtk_manual[vtk_uid])();
  }

  /** Activate or deactivate a VTK callback */
  void switch_vtk(std::string const &vtk_uid, int status) override {
    if (m_vtk_manual.find(vtk_uid) != m_vtk_manual.end()) {
      throw std::runtime_error("VTKOutput object " + vtk_uid +
                               " is a manual observable");
    }
    if (m_vtk_auto.find(vtk_uid) == m_vtk_auto.end()) {
      throw std::runtime_error("VTKOutput object " + vtk_uid +
                               " doesn't exist");
    }
    m_vtk_auto[vtk_uid].second = status;
  }

  ~LBWalberlaImpl() override = default;
};
} // namespace walberla

#endif // LB_WALBERLA_H
