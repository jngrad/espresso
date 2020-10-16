
#ifndef LB_WALBERLA_H
#define LB_WALBERLA_H

#include "config.hpp"

#ifdef LB_WALBERLA
#include "blockforest/Initialization.h"
#include "blockforest/StructuredBlockForest.h"
#include "blockforest/communication/UniformBufferedScheme.h"
#include "boost/optional.hpp"
#include "boost/tuple/tuple.hpp"
#include "boundary/BoundaryHandling.h"
#include "core/mpi/Environment.h"
#include "field/FlagField.h"
#include "field/GhostLayerField.h"
#include "field/adaptors/GhostLayerFieldAdaptor.h"
#include "field/vtk/FlagFieldCellFilter.h"
#include "field/vtk/VTKWriter.h"
#include "lbm/boundary/UBB.h"
#include "lbm/field/Adaptors.h"
#include "lbm/field/PdfField.h"
#include "lbm/lattice_model/CollisionModel.h"
#include "lbm/lattice_model/D3Q19.h"
#include "lbm/vtk/all.h"
#include "timeloop/SweepTimeloop.h"
#include "utils/Vector.hpp"
#include "utils/interpolation/bspline_3d.hpp"
#include "utils/math/make_lin_space.hpp"

#include "boundary/BoundaryHandling.h"

#include "core/mpi/Environment.h"
#include "core/mpi/MPIManager.h"
#include "core/mpi/MPITextFile.h"
#include "core/mpi/Reduce.h"

#include "domain_decomposition/SharedSweep.h"

#include "field/AddToStorage.h"
#include "field/FlagField.h"
#include "field/adaptors/AdaptorCreators.h"
#include "field/communication/PackInfo.h"
#include "lbm/boundary/NoSlip.h"
#include "lbm/boundary/UBB.h"
#include "lbm/communication/PdfFieldPackInfo.h"
#include "lbm/field/AddToStorage.h"
#include "lbm/field/PdfField.h"
#include "lbm/lattice_model/D3Q19.h"
#include "lbm/sweeps/CellwiseSweep.h"

#include "stencil/D3Q27.h"

#include "timeloop/SweepTimeloop.h"

#include "LbWalberlaBase.hpp"

namespace walberla {

// Flags marking fluid and boundaries
const FlagUID Fluid_flag("fluid");
const FlagUID UBB_flag("velocity bounce back");

// Vector conversion helpers
inline Utils::Vector3d to_vector3d(const Vector3<real_t> v) {
  return Utils::Vector3d{v[0], v[1], v[2]};
}
inline Vector3<real_t> to_vector3(const Utils::Vector3d v) {
  return Vector3<real_t>{v[0], v[1], v[2]};
}
inline Utils::Vector6d to_vector6d(const Matrix3<real_t> m) {
  return Utils::Vector6d{m[0], m[3], m[4], m[6], m[7], m[8]};
}
inline Utils::Vector3i to_vector3i(const std::array<int, 3> v) {
  return Utils::Vector3i{v[0], v[1], v[2]};
}
/** Sweep that swaps force_toe_be_applied and last_applied_force
and resets force_to_be_applied to the global external force
*/
template <typename PdfField, typename ForceField, typename BoundaryHandling>
class ResetForce {
public:
  ResetForce(const BlockDataID &pdf_field_id,
             const BlockDataID &last_applied_force_field_id,
             const BlockDataID &force_to_be_applied_id,
             const BlockDataID &boundary_handling_id)
      : m_pdf_field_id(pdf_field_id),
        m_last_applied_force_field_id(last_applied_force_field_id),
        m_force_to_be_applied_id(force_to_be_applied_id),
        m_boundary_handling_id(boundary_handling_id),
        m_ext_force(Vector3<real_t>{0, 0, 0}){};

  void set_ext_force(const Utils::Vector3d &ext_force) {
    m_ext_force = to_vector3(ext_force);
  }

  Utils::Vector3d get_ext_force() const { return to_vector3d(m_ext_force); };

  void operator()(IBlock *block) {
    PdfField *pdf_field = block->template getData<PdfField>(m_pdf_field_id);
    ForceField *force_field =
        block->template getData<ForceField>(m_last_applied_force_field_id);
    ForceField *force_to_be_applied =
        block->template getData<ForceField>(m_force_to_be_applied_id);
    BoundaryHandling *boundary_handling =
        block->template getData<BoundaryHandling>(m_boundary_handling_id);

    force_field->swapDataPointers(force_to_be_applied);

    WALBERLA_FOR_ALL_CELLS_INCLUDING_GHOST_LAYER_XYZ(force_field, {
      Cell cell(x, y, z);
      for (int i : {0, 1, 2})
        force_field->get(x, y, z, i) += m_ext_force[i];
    });
    WALBERLA_FOR_ALL_CELLS_INCLUDING_GHOST_LAYER_XYZ(force_to_be_applied, {
      Cell cell(x, y, z);
      for (int i : {0, 1, 2})
        force_to_be_applied->get(cell, i) = real_t{0};
    });
  }

private:
  const BlockDataID m_pdf_field_id;
  const BlockDataID m_last_applied_force_field_id;
  const BlockDataID m_force_to_be_applied_id;
  const BlockDataID m_boundary_handling_id;
  Vector3<real_t> m_ext_force;
};

/** Class that runs and controls the LB on WaLBerla
 */
template <typename LatticeModel> class LbWalberla : public LbWalberlaBase {
protected:
  // Type definitions
  using VectorField = GhostLayerField<real_t, 3>;
  using FlagField = walberla::FlagField<walberla::uint8_t>;
  using PdfField = lbm::PdfField<LatticeModel>;

  using UBB = lbm::UBB<LatticeModel, uint8_t, true, true>;
  using Boundaries =
      BoundaryHandling<FlagField, typename LatticeModel::Stencil, UBB>;

  // Adaptors
  using DensityAdaptor = typename lbm::Adaptor<LatticeModel>::Density;
  using VelocityAdaptor = typename lbm::Adaptor<LatticeModel>::VelocityVector;

  /** VTK writers that are executed automatically */
  std::map<std::string, std::pair<std::shared_ptr<vtk::VTKOutput>, bool>>
      m_vtk_auto;
  /** VTK writers that are executed manually */
  std::map<std::string, std::shared_ptr<vtk::VTKOutput>> m_vtk_manual;

  template <typename VectorField> class force_vector_adaptor_function {
  public:
    typedef VectorField basefield_t;
    typedef typename basefield_t::const_base_iterator basefield_iterator;
    typedef Vector3<real_t> value_type;

    static const uint_t F_SIZE = 1u;

    value_type operator()(const basefield_t &baseField, cell_idx_t x,
                          cell_idx_t y, cell_idx_t z,
                          cell_idx_t /*f*/ = 0) const {
      return baseField.get(x, y, z);
    }

    value_type operator()(const basefield_iterator &it) const {
      auto const *baseFieldPtr =
          dynamic_cast<const basefield_t *>(it.getField());
      return baseFieldPtr->get(it);
    }
  };

  // Member variables
  // For unit conversions
  double m_agrid;
  double m_tau;

  double m_kT;
  double m_density;

  Utils::Vector3i m_grid_dimensions;
  int m_n_ghost_layers;

  // Block data access handles
  BlockDataID m_pdf_field_id;
  BlockDataID m_flag_field_id;

  BlockDataID m_last_applied_force_field_id;
  BlockDataID m_force_to_be_applied_id;

  BlockDataID m_velocity_adaptor_id;

  BlockDataID m_density_adaptor_id;

  BlockDataID m_boundary_handling_id;

  using Communicator = blockforest::communication::UniformBufferedScheme<
      typename stencil::D3Q27>;
  std::shared_ptr<Communicator> m_communication;

  // Block forest
  std::shared_ptr<blockforest::StructuredBlockForest> m_blocks;

  // time loop
  std::shared_ptr<timeloop::SweepTimeloop> m_time_loop;

  // MPI
  std::shared_ptr<mpi::Environment> m_env;

  // Lattice model
  std::shared_ptr<LatticeModel> m_lattice_model;

  // ResetForce sweep + external force handling
  std::shared_ptr<ResetForce<PdfField, VectorField, Boundaries>> m_reset_force;

  // Helpers to retrieve blocks and cells
  struct BlockAndCell {
    IBlock *block;
    Cell cell;
  };
  template <typename PosVector>
  IBlock *get_block_extended(const PosVector &pos) const {
    for (auto b = m_blocks->begin(); b != m_blocks->end(); ++b) {
      if (b->getAABB()
              .getExtended(m_n_ghost_layers)
              .contains(real_c(pos[0]), real_c(pos[1]), real_c(pos[2]))) {
        return &(*b);
      }
    }
    // Cell not in local blocks
    return {};
  }
  boost::optional<BlockAndCell>
  get_block_and_cell(const Utils::Vector3i &node,
                     bool consider_ghost_layers) const {
    // Get block and local cell
    Cell global_cell{uint_c(node[0]), uint_c(node[1]), uint_c(node[2])};
    auto block = m_blocks->getBlock(global_cell, 0);
    // Return if we don't have the cell
    if (consider_ghost_layers and !block) {
      // Try to find a block which has the cell as ghost layer
      block = get_block_extended(node);
    }
    if (!block)
      return {boost::none};

    // Transform coords to block local
    Cell local_cell;
    m_blocks->transformGlobalToBlockLocalCell(local_cell, *block, global_cell);
    return {{block, local_cell}};
  }

  IBlock *get_block(const Utils::Vector3d &pos,
                    bool consider_ghost_layers) const {
    // Get block
    auto block =
        m_blocks->getBlock(real_c(pos[0]), real_c(pos[1]), real_c(pos[2]));
    if (consider_ghost_layers and !block) {
      block = get_block_extended(pos);
    }
    return block;
  };

  // Boundary handling
  class LB_boundary_handling {
  public:
    LB_boundary_handling(const BlockDataID &flag_field_id,
                         const BlockDataID &pdf_field_id)
        : m_flag_field_id(flag_field_id), m_pdf_field_id(pdf_field_id) {}

    Boundaries *operator()(IBlock *const block) {

      FlagField *flag_field =
          block->template getData<FlagField>(m_flag_field_id);
      PdfField *pdf_field = block->template getData<PdfField>(m_pdf_field_id);

      // const uint8_t fluid = flag_field->registerFlag(Fluid_flag);
      const auto fluid = flag_field->flagExists(Fluid_flag)
                             ? flag_field->getFlag(Fluid_flag)
                             : flag_field->registerFlag(Fluid_flag);

      return new Boundaries(
          "boundary handling", flag_field, fluid,
          UBB("velocity bounce back", UBB_flag, pdf_field, nullptr));
    }

  private:
    const BlockDataID m_flag_field_id;
    const BlockDataID m_pdf_field_id;
  };

public:
  LbWalberla(double viscosity, double density, double agrid, double tau,
             const Utils::Vector3d &box_dimensions,
             const Utils::Vector3i &node_grid, int n_ghost_layers) {
    m_agrid = agrid;
    m_tau = tau;
    m_density = density;
    m_n_ghost_layers = n_ghost_layers;
    m_kT = 0;

    if (m_n_ghost_layers <= 0)
      throw std::runtime_error("At least one ghost layer must be used");

    for (int i = 0; i < 3; i++) {
      if (fabs(round(box_dimensions[i] / agrid) * agrid - box_dimensions[i]) /
              box_dimensions[i] >
          std::numeric_limits<double>::epsilon()) {
        throw std::runtime_error(
            "Box length not commensurate with agrid in direction " +
            std::to_string(i));
      }
      m_grid_dimensions[i] = int(std::round(box_dimensions[i] / agrid));
      if (m_grid_dimensions[i] % node_grid[i] != 0) {
        printf("Grid dimension: %d, node grid %d\n", m_grid_dimensions[i],
               node_grid[i]);
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

    m_last_applied_force_field_id = field::addToStorage<VectorField>(
        m_blocks, "force field", real_t{0}, field::fzyx, m_n_ghost_layers);
    m_force_to_be_applied_id = field::addToStorage<VectorField>(
        m_blocks, "force field", real_t{0}, field::fzyx, m_n_ghost_layers);
  };

  void setup_with_valid_lattice_model() {
    m_pdf_field_id = lbm::addPdfFieldToStorage(
        m_blocks, "pdf field", *(m_lattice_model.get()),
        to_vector3(Utils::Vector3d{}), (real_t)1.0, m_n_ghost_layers);

    m_flag_field_id = field::addFlagFieldToStorage<FlagField>(
        m_blocks, "flag field", m_n_ghost_layers);

    m_boundary_handling_id = m_blocks->addBlockData<Boundaries>(
        LB_boundary_handling(m_flag_field_id, m_pdf_field_id),
        "boundary handling");

    clear_boundaries();

    // 1 timestep each time the integrate function is called
    m_time_loop = std::make_shared<timeloop::SweepTimeloop>(
        m_blocks->getBlockStorage(), 1);

    // sets up the communication
    m_communication = std::make_shared<Communicator>(m_blocks);
    m_communication->addPackInfo(
        std::make_shared<field::communication::PackInfo<PdfField>>(
            m_pdf_field_id));
    m_communication->addPackInfo(
        std::make_shared<field::communication::PackInfo<VectorField>>(
            m_last_applied_force_field_id));
    m_reset_force =
        std::make_shared<ResetForce<PdfField, VectorField, Boundaries>>(
            m_pdf_field_id, m_last_applied_force_field_id,
            m_force_to_be_applied_id, m_boundary_handling_id);

    //    m_time_loop->add() << timeloop::BeforeFunction(communication,
    //                                                   "communication")
    m_time_loop->add() << timeloop::Sweep(
        Boundaries::getBlockSweep(m_boundary_handling_id), "boundary handling");
    m_time_loop->add() << timeloop::Sweep(makeSharedSweep(m_reset_force),
                                          "Reset force fields");
    m_time_loop->add() << timeloop::Sweep(
                              typename LatticeModel::Sweep(m_pdf_field_id),
                              "LB stream & collide")
                       << timeloop::AfterFunction(*m_communication,
                                                  "communication");

    m_velocity_adaptor_id = field::addFieldAdaptor<VelocityAdaptor>(
        m_blocks, m_pdf_field_id, "velocity adaptor");

    m_density_adaptor_id = field::addFieldAdaptor<DensityAdaptor>(
        m_blocks, m_pdf_field_id, "density adaptor");

    // Synchronize ghost layers
    (*m_communication)();
  };
  std::shared_ptr<LatticeModel> get_lattice_model() { return m_lattice_model; };

  void integrate() override {
    m_time_loop->singleStep();
    for (auto it = m_vtk_auto.begin(); it != m_vtk_auto.end(); ++it) {
      if (it->second.second)
        vtk::writeFiles(it->second.first)();
    }
  };

  void ghost_communication() override { (*m_communication)(); }

  template <typename Function>
  void interpolate_bspline_at_pos(Utils::Vector3d pos, Function f) const {
    Utils::Interpolation::bspline_3d<2>(
        pos, f, Utils::Vector3d{1.0, 1.0, 1.0}, // grid spacing
        Utils::Vector3d::broadcast(.5));        // offset
  }

  // Velocity
  boost::optional<Utils::Vector3d>
  get_node_velocity(const Utils::Vector3i &node,
                    bool consider_ghosts = false) const override {
    boost::optional<bool> is_boundary =
        get_node_is_boundary(node, consider_ghosts);
    if (is_boundary)    // is info available locally
      if (*is_boundary) // is the node a boundary
        return get_node_velocity_at_boundary(node);
    auto const bc = get_block_and_cell(node, consider_ghosts);
    if (!bc)
      return {};
    auto const &vel_adaptor =
        (*bc).block->template getData<VelocityAdaptor>(m_velocity_adaptor_id);
    return {to_vector3d(vel_adaptor->get((*bc).cell))};
  };
  bool set_node_velocity(const Utils::Vector3i &node,
                         const Utils::Vector3d &v) override {
    auto bc = get_block_and_cell(node, false);
    if (!bc)
      return false;
    auto pdf_field = (*bc).block->template getData<PdfField>(m_pdf_field_id);
    const real_t density = pdf_field->getDensity((*bc).cell);
    pdf_field->setDensityAndVelocity(
        (*bc).cell, Vector3<double>{v[0], v[1], v[2]}, density);
    return true;
  };

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
  };

  // Local force
  bool add_force_at_pos(const Utils::Vector3d &pos,
                        const Utils::Vector3d &force) override {
    if (!pos_in_local_halo(pos))
      return false;
    auto force_at_node = [this, force](const std::array<int, 3> node,
                                       double weight) {
      auto const bc = get_block_and_cell(to_vector3i(node), true);
      if (bc) {
        auto force_field = (*bc).block->template getData<VectorField>(
            m_force_to_be_applied_id);
        for (int i : {0, 1, 2})
          force_field->get((*bc).cell, i) +=
              real_t{force[i] * weight / m_density};
      }
    };
    interpolate_bspline_at_pos(pos, force_at_node);
    return true;
  };

  boost::optional<Utils::Vector3d>
  get_node_force_to_be_applied(const Utils::Vector3i &node) const override {
    auto const bc = get_block_and_cell(node, true);
    if (!bc)
      return {};

    auto const &force_field =
        (*bc).block->template getData<VectorField>(m_force_to_be_applied_id);
    return Utils::Vector3d{{force_field->get((*bc).cell, 0),
                            force_field->get((*bc).cell, 1),
                            force_field->get((*bc).cell, 2)}} *
           m_density;
  };

  boost::optional<Utils::Vector3d>
  get_node_last_applied_force(const Utils::Vector3i &node,
                              bool consider_ghosts = false) const override {
    auto const bc = get_block_and_cell(node, consider_ghosts);
    if (!bc)
      return {};

    auto const &force_field = (*bc).block->template getData<VectorField>(
        m_last_applied_force_field_id);
    return Utils::Vector3d{{force_field->get((*bc).cell, 0),
                            force_field->get((*bc).cell, 1),
                            force_field->get((*bc).cell, 2)}} *
           m_density;
  };

  // Density
  bool set_node_density(const Utils::Vector3i &node, double density) override {
    auto bc = get_block_and_cell(node, false);
    if (!bc)
      return false;

    auto pdf_field = (*bc).block->template getData<PdfField>(m_pdf_field_id);
    auto const &vel_adaptor =
        (*bc).block->template getData<VelocityAdaptor>(m_velocity_adaptor_id);
    Vector3<double> v = vel_adaptor->get((*bc).cell);

    pdf_field->setDensityAndVelocity(
        (*bc).cell, Vector3<double>{v[0], v[1], v[2]}, density / m_density);

    return true;
  };

  boost::optional<double>
  get_node_density(const Utils::Vector3i &node) const override {
    auto bc = get_block_and_cell(node, false);
    if (!bc)
      return {boost::none};

    auto pdf_field = (*bc).block->template getData<PdfField>(m_pdf_field_id);

    return {m_density * pdf_field->getDensity((*bc).cell)};
  };

  // Boundary related
  boost::optional<Utils::Vector3d>
  get_node_velocity_at_boundary(const Utils::Vector3i &node) const override {
    auto bc = get_block_and_cell(node, true);
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
  };
  bool set_node_velocity_at_boundary(const Utils::Vector3i &node,
                                     const Utils::Vector3d &v) override {
    auto bc = get_block_and_cell(node, true);
    if (!bc)
      return false;

    const typename UBB::Velocity velocity(real_c(v[0]), real_c(v[1]),
                                          real_c(v[2]));

    Boundaries *boundary_handling =
        (*bc).block->template getData<Boundaries>(m_boundary_handling_id);
    boundary_handling->forceBoundary(UBB_flag, bc->cell[0], bc->cell[1],
                                     bc->cell[2], velocity);
    return true;
  };
  boost::optional<Utils::Vector3d>
  get_node_boundary_force(const Utils::Vector3i &node) const override {
    auto bc = get_block_and_cell(node, true); // including ghosts
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
    return {m_density * to_vector3d(ubb.getForce((*bc).cell.x(), (*bc).cell.y(),
                                                 (*bc).cell.z()))};
  };
  bool remove_node_from_boundary(const Utils::Vector3i &node) override {
    auto bc = get_block_and_cell(node, true);
    if (!bc)
      return false;
    Boundaries *boundary_handling =
        (*bc).block->template getData<Boundaries>(m_boundary_handling_id);
    boundary_handling->removeBoundary((*bc).cell[0], (*bc).cell[1],
                                      (*bc).cell[2]);
    return true;
  };
  boost::optional<bool>
  get_node_is_boundary(const Utils::Vector3i &node,
                       bool consider_ghosts = false) const override {
    auto bc = get_block_and_cell(node, consider_ghosts);
    if (!bc)
      return {boost::none};

    Boundaries *boundary_handling =
        (*bc).block->template getData<Boundaries>(m_boundary_handling_id);
    return {boundary_handling->isBoundary((*bc).cell)};
  };
  void clear_boundaries() override {
    const CellInterval &domain_bb_in_global_cell_coordinates =
        m_blocks->getCellBBFromAABB(
            m_blocks->begin()->getAABB().getExtended(1.0 * n_ghost_layers()));
    for (auto block = m_blocks->begin(); block != m_blocks->end(); ++block) {

      Boundaries *boundary_handling =
          block->template getData<Boundaries>(m_boundary_handling_id);

      CellInterval domain_bb(domain_bb_in_global_cell_coordinates);
      m_blocks->transformGlobalToBlockLocalCellInterval(domain_bb, *block);

      boundary_handling->fillWithDomain(domain_bb);
    }
  };

  // Pressure tensor
  boost::optional<Utils::Vector6d>
  get_node_pressure_tensor(const Utils::Vector3i &node) const override {
    auto bc = get_block_and_cell(node, false);
    if (!bc)
      return {boost::none};
    auto pdf_field = (*bc).block->template getData<PdfField>(m_pdf_field_id);
    return to_vector6d(pdf_field->getPressureTensor((*bc).cell));
  };

  // Global momentum
  Utils::Vector3d get_momentum() const override {
    Vector3<real_t> mom;
    for (auto block_it = m_blocks->begin(); block_it != m_blocks->end();
         ++block_it) {
      auto pdf_field = block_it->template getData<PdfField>(m_pdf_field_id);
      Vector3<real_t> local_v;
      WALBERLA_FOR_ALL_CELLS_XYZ(pdf_field, {
        double local_dens = pdf_field->getDensityAndVelocity(local_v, x, y, z);
        mom += local_dens * local_v;
      });
    }
    return m_density * to_vector3d(mom);
  };
  // Global external force
  void set_external_force(const Utils::Vector3d &ext_force) override {
    m_reset_force->set_ext_force(ext_force / m_density);
  };
  Utils::Vector3d get_external_force() const override {
    return m_reset_force->get_ext_force() * m_density;
  };

  // Global parameters
  //  double get_viscosity() const override {
  //    return m_lattice_model->collisionModel().viscosity();
  //  };

  double get_tau() const override { return m_tau; };
  double get_kT() const override { return m_kT; };

  // Grid, domain, halo
  int n_ghost_layers() const override { return m_n_ghost_layers; };
  Utils::Vector3i get_grid_dimensions() const override {
    return m_grid_dimensions;
  }
  double get_grid_spacing() const override { return m_agrid; };
  std::pair<Utils::Vector3d, Utils::Vector3d>
  get_local_domain() const override {
    // We only have one block per mpi rank
    assert(++(m_blocks->begin()) == m_blocks->end());

    auto const ab = m_blocks->begin()->getAABB();
    return {to_vector3d(ab.min()), to_vector3d(ab.max())};
  };

  bool node_in_local_domain(const Utils::Vector3i &node) const override {
    // Note: Lattice constant =1, cell centers offset by .5
    return get_block_and_cell(node, false) != boost::none;
  };
  bool node_in_local_halo(const Utils::Vector3i &node) const override {
    return get_block_and_cell(node, true) != boost::none;
  };
  bool pos_in_local_domain(const Utils::Vector3d &pos) const override {
    return get_block(pos, false) != nullptr;
  };
  bool pos_in_local_halo(const Utils::Vector3d &pos) const override {
    return get_block(pos, true) != nullptr;
  };

  std::vector<std::pair<Utils::Vector3i, Utils::Vector3d>>
  node_indices_positions(bool include_ghosts = false) const override {
    int ghost_offset = 0;
    if (include_ghosts)
      ghost_offset = m_n_ghost_layers;
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
  };

  void create_vtk(unsigned delta_N, unsigned initial_count,
                  unsigned flag_observables, std::string const &identifier,
                  std::string const &base_folder, std::string const &prefix) {
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
          make_shared<lbm::DensityVTKWriter<LatticeModel, float>>(
              m_pdf_field_id, "DensityFromPDF"));
    }
    if (static_cast<unsigned>(OutputVTK::velocity_vector) & flag_observables) {
      pdf_field_vtk->addCellDataWriter(
          make_shared<lbm::VelocityVTKWriter<LatticeModel, float>>(
              m_pdf_field_id, "VelocityFromPDF"));
    }
    if (static_cast<unsigned>(OutputVTK::pressure_tensor) & flag_observables) {
      pdf_field_vtk->addCellDataWriter(
          make_shared<lbm::PressureTensorVTKWriter<LatticeModel, float>>(
              m_pdf_field_id, "PressureTensorFromPDF"));
    }

    // register object
    if (delta_N) {
      m_vtk_auto[vtk_uid] = {pdf_field_vtk, true};
    } else {
      m_vtk_manual[vtk_uid] = pdf_field_vtk;
    }
  }

  /** Manually call a VTK callback */
  void write_vtk(std::string const &vtk_uid) {
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
  virtual void switch_vtk(std::string const &vtk_uid, int status) {
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

  /** @brief call, if the lattice model was changed */
  void on_lattice_model_change() {
    for (auto b = m_blocks->begin(); b != m_blocks->end(); ++b) {
      auto pdf_field = b->template getData<PdfField>(m_pdf_field_id);
      pdf_field->resetLatticeModel(*m_lattice_model);
      pdf_field->latticeModel().configure(*b, *m_blocks);
    }
  }

  ~LbWalberla() override = default;
};
} // namespace walberla
#endif // LB_WALBERLA

#endif // LB_WALBERLA_H
