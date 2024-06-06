/*
 * Copyright (C) 2010-2022 The ESPResSo project
 * Copyright (C) 2002,2003,2004,2005,2006,2007,2008,2009,2010
 *   Max-Planck-Institute for Polymer Research, Theory Group
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

#include "cell_system/ParticleDecomposition.hpp"

#include "cell_system/Cell.hpp"

#include "BoxGeometry.hpp"
#include "Particle.hpp"
#include "ghosts.hpp"

#include <utils/Vector.hpp>

#include <boost/mpi/communicator.hpp>

#include <optional>
#include <span>
#include <utility>
#include <vector>

/**
 * @brief Atom decomposition cell system.
 *
 * This implements a distributed particle storage
 * by just evenly distributing the particles over
 * all part-taking nodes. Pairs are found by just
 * considering all pairs independent of logical or
 * physical location, it has therefore quadratic time
 * complexity in the number of particles.
 *
 * For a more detailed discussion please see @cite plimpton95a.
 */
class AtomDecomposition : public ParticleDecomposition {
  boost::mpi::communicator m_comm;
  std::vector<Cell> cells;

  std::vector<Cell *> m_local_cells;
  std::vector<Cell *> m_ghost_cells;

  GhostCommunicator m_exchange_ghosts_comm;
  GhostCommunicator m_collect_ghost_force_comm;

  BoxGeometry const &m_box;

public:
  AtomDecomposition(BoxGeometry const &m_box);
  AtomDecomposition(boost::mpi::communicator comm, BoxGeometry const &box_geo);

  void resort(bool global_flag, std::vector<ParticleChange> &diff) override;

  GhostCommunicator const &exchange_ghosts_comm() const override {
    return m_exchange_ghosts_comm;
  }
  GhostCommunicator const &collect_ghost_force_comm() const override {
    return m_collect_ghost_force_comm;
  }

  std::span<Cell *const> local_cells() const override { return m_local_cells; }
  std::span<Cell *const> ghost_cells() const override { return m_ghost_cells; }

  /* Getter needed for HybridDecomposition */
  auto const &get_local_cells() const { return m_local_cells; }
  auto const &get_ghost_cells() const { return m_ghost_cells; }

  /**
   * @brief Determine which cell a particle id belongs to.
   *
   * Since there is only one local cell this is trivial.
   *
   * @param p Particle to find cell for.
   * @return Pointer to cell or nullptr if not local.
   */
  Cell *particle_to_cell(Particle const &p) override {
    return id_to_cell(p.id());
  }
  Cell const *particle_to_cell(Particle const &p) const override {
    return id_to_cell(p.id());
  }

  Utils::Vector3d max_cutoff() const override;
  Utils::Vector3d max_range() const override;

  /* Return true if minimum image convention is
   * needed for distance calculation. */
  std::optional<BoxGeometry> minimum_image_distance() const override {
    return m_box;
  }

  BoxGeometry const &box() const override { return m_box; }

private:
  /**
   * @brief Find cell for id.
   * @param id to find cell for.
   * @return Cell for id.
   */
  Cell *id_to_cell(int id) {
    return has_id(id) ? std::addressof(local()) : nullptr;
  }

  Cell const *id_to_cell(int id) const {
    return has_id(id) ? std::addressof(local()) : nullptr;
  }

  /**
   * @brief Get the local cell.
   */
  Cell &local() { return cells.at(static_cast<unsigned int>(m_comm.rank())); }
  Cell const &local() const {
    return cells.at(static_cast<unsigned int>(m_comm.rank()));
  }

  void configure_neighbors();
  GhostCommunicator prepare_comm();

  /**
   * @brief Setup ghost communicators.
   */
  void configure_comms();

  /**
   * @brief Mark local and ghost cells.
   */
  void mark_cells();

  /**
   * @brief Determine which rank owns a particle id.
   */
  int id_to_rank(int id) const { return id % m_comm.size(); }

  /**
   * @brief Determine if this rank owns a particle id.
   */
  bool has_id(int id) const { return id_to_rank(id) == m_comm.rank(); }
};
