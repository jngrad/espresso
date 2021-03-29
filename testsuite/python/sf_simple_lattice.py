#
# Copyright (C) 2017-2021 The ESPResSo project
#
# This file is part of ESPResSo.
#
# ESPResSo is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published byss
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
import numpy as np
import itertools


class StructureFactorTest(ut.TestCase):

    box_l = 16
    part_ty = 0
    sf_order = 16
    system = espressomd.System(box_l=[box_l, box_l, box_l])

    def tearDown(self):
        self.system.part.clear()

    def test_sc(self):
        """Check simple cubic lattice"""
        xen = range(0, self.box_l, 4)
        for i, j, k in itertools.product(xen, repeat=3):
            self.system.part.add(type=self.part_ty, pos=(i, j, k))
        wavevectors, intensities = self.system.analysis.structure_factor(
            sf_types=[self.part_ty], sf_order=self.sf_order)
        intensities_int = np.around(intensities).astype(int)
        np.testing.assert_array_almost_equal(intensities, intensities_int)
        intensities_ref = np.zeros(intensities.shape)
        intensities_ref[np.nonzero(intensities_int)] = len(self.system.part)
        np.testing.assert_array_equal(intensities_int, intensities_ref)
        # F = f for all planes
        peaks_ref = list(range(1, self.box_l + 1))
        peaks_ref.remove(7)
        peaks_ref.remove(15)
        peaks = (wavevectors[np.nonzero(intensities_int)] * 2 / np.pi)**2
        np.testing.assert_array_almost_equal(peaks, peaks_ref)

    def test_cco(self):
        """Check c-centered orthorhombic lattice"""
        xen = range(0, self.box_l, 4)
        for i, j, k in itertools.product(xen, repeat=3):
            self.system.part.add(type=self.part_ty, pos=(i, j, k))
        for i, j, k in itertools.product(xen, repeat=3):
            self.system.part.add(type=self.part_ty, pos=(i + 2, j + 2, k))
        wavevectors, intensities = self.system.analysis.structure_factor(
            sf_types=[self.part_ty], sf_order=self.sf_order)
        intensities_int = np.around(intensities).astype(int)
        print(intensities)
        print((wavevectors[np.nonzero(intensities_int)] * 2 / np.pi)**2)
        np.testing.assert_array_almost_equal(intensities, intensities_int)
        intensities_ref = np.zeros(intensities.shape)
        intensities_ref[np.nonzero(intensities_int)] = len(self.system.part)
        np.testing.assert_array_equal(intensities_int, intensities_ref)
        # (h+k) even => F = 2f, otherwise F = 0
        peaks_ref = np.arange(2, self.box_l + 1, 2)
        peaks = (wavevectors[np.nonzero(intensities_int)] * 2 / np.pi)**2
        np.testing.assert_array_almost_equal(peaks, peaks_ref)

    def test_bcc(self):
        """Check body-centered cubic lattice"""
        xen = range(0, self.box_l, 4)
        for i, j, k in itertools.product(xen, repeat=3):
            self.system.part.add(type=self.part_ty, pos=(i, j, k))
            self.system.part.add(type=self.part_ty, pos=(i + 2, j + 2, k + 2))
        wavevectors, intensities = self.system.analysis.structure_factor(
            sf_types=[self.part_ty], sf_order=self.sf_order)
        intensities_int = np.around(intensities).astype(int)
        np.testing.assert_array_almost_equal(intensities, intensities_int)
        intensities_ref = np.zeros(intensities.shape)
        intensities_ref[np.nonzero(intensities_int)] = len(self.system.part)
        np.testing.assert_array_equal(intensities_int, intensities_ref)
        # (h+k+l) even => F = 2f, otherwise F = 0
        peaks_ref = np.arange(2, self.box_l + 1, 2)
        peaks = (wavevectors[np.nonzero(intensities_int)] * 2 / np.pi)**2
        np.testing.assert_array_almost_equal(peaks, peaks_ref)

    def test_fcc(self):
        """Check face-centered cubic lattice"""
        xen = range(0, self.box_l, 4)
        for i, j, k in itertools.product(xen, repeat=3):
            self.system.part.add(type=self.part_ty, pos=(i, j, k))
            self.system.part.add(type=self.part_ty, pos=(i + 2, j + 2, k))
            self.system.part.add(type=self.part_ty, pos=(i + 2, j, k + 2))
            self.system.part.add(type=self.part_ty, pos=(i, j + 2, k + 2))
        wavevectors, intensities = self.system.analysis.structure_factor(
            sf_types=[self.part_ty], sf_order=self.sf_order)
        intensities_int = np.around(intensities).astype(int)
        np.testing.assert_array_almost_equal(intensities, intensities_int)
        intensities_ref = np.zeros(intensities.shape)
        intensities_ref[np.nonzero(intensities_int)] = len(self.system.part)
        np.testing.assert_array_equal(intensities_int, intensities_ref)
        # (h,k,l) all even or odd => F = 4f, otherwise F = 0
        peaks_ref = [3, 4, 8, 11, 12, 16]
        peaks = (wavevectors[np.nonzero(intensities_int)] * 2 / np.pi)**2
        np.testing.assert_array_almost_equal(peaks, peaks_ref)

    def test_tetragonal(self):
        """Check tetragonal lattice"""
        xen = range(0, self.box_l, 1)
        yen = range(0, self.box_l, 2)
        zen = range(0, self.box_l, 4)
        for i, j, k in itertools.product(xen, yen, zen):
            self.system.part.add(type=self.part_ty, pos=(i, j, k))
        wavevectors, intensities = self.system.analysis.structure_factor(
            sf_types=[self.part_ty], sf_order=self.sf_order)
        intensities_int = np.around(intensities).astype(int)
        print(intensities_int)
        print(intensities[np.nonzero(intensities_int)])
        print((wavevectors[np.nonzero(intensities_int)] * 2 / np.pi)**2)

    def test_exceptions(self):
        with self.assertRaisesRegex(ValueError, 'order has to be a strictly positive number'):
            self.system.analysis.structure_factor(sf_types=[0], sf_order=0)


if __name__ == "__main__":
    ut.main()
