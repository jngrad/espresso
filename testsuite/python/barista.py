#
# Copyright (C) 2023 The ESPResSo project
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
#

import unittest as ut

import espressomd
import espressomd.barista


class Test(ut.TestCase):
    system = espressomd.System(box_l=[1., 1., 1.])

    def test_serving(self):
        beverage = espressomd.barista.FancyBrew(base="lungo")
        self.assertEqual(beverage.base, "lungo")
        self.assertEqual(beverage.toppings, [])
        beverage.set_toppings(toppings=["milk", "sugar"])
        self.assertEqual(beverage.toppings, ["milk", "sugar"])
        beverage.set_toppings(toppings=["sugar"])
        self.assertEqual(beverage.toppings, ["sugar"])
        out = beverage.brew_coffee()
        self.assertEqual(out, "Here is your lungo with sugar")

    def test_cleaning(self):
        n_nodes = self.system.cell_system.get_state()["n_nodes"]
        beverage = espressomd.barista.FancyBrew(base="lungo")
        beverage.set_toppings(toppings=["sugar"])
        water_used = beverage.run_cleaning_cycle()
        self.assertAlmostEqual(water_used, n_nodes * 0.1, delta=1e-5)
        beverage.set_toppings(toppings=["sugar", "milk"])
        water_used = beverage.run_cleaning_cycle()
        self.assertAlmostEqual(water_used, n_nodes * 0.3, delta=1e-5)

    def test_exceptions(self):
        with self.assertRaisesRegex(ValueError, "Base not recognized: tea"):
            espressomd.barista.FancyBrew(base="tea")
        with self.assertRaisesRegex(ValueError, "Topping not recognized: chocolate"):
            espressomd.barista.FancyBrew(base="lungo", toppings=["chocolate"])

if __name__ == "__main__":
    ut.main()
