# Copyright (C) 2019 The ESPResSo project
#
# This file is part of ESPResSo.
#
# ESPResSo is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# ESPResSo is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <http://www.gnu.org/licenses/>.
import espressomd  # pylint: disable=import-error
import unittest


class CaseSystem(object):

    """
    Utility class that initializes the espresso system and sets the random seed.
    """
    system = espressomd.System(box_l=3 * [1.])
    system.seed = system.cell_system.get_state()['n_nodes'] * [1234]


class TestCaseSystem(unittest.TestCase, CaseSystem):
    pass
