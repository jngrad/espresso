#
# Copyright (C) 2024-2025 The ESPResSo project
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
"""
P3M method implemented in Python.
Limitations:

- only 1 MPI rank
- cubic box
- energy kernels are incorrect
- pressure kernels are not implemented
- tuning algorithm not implemented
- non-metallic boundaries not implemented
"""

import math
import numpy as np


class MpiComm:
    def size(self):
        return 1

    def rank(self):
        return 0


USE_ERFC_APPROXIMATION = True
P3M_EPSILON_METALLIC = 0.
P3M_MESHOFF = 0.5
COLUMN_MAJOR = 0
ROW_MAJOR = 1

comm_cart = MpiComm()


def AS_erfc_part(d: float):
    a1 = 0.254829592
    a2 = -0.284496736
    a3 = 1.421413741
    a4 = -1.453152027
    a5 = 1.061405429
    p = 0.3275911
    t = 1.0 / (1.0 + p * d)
    return t * (a1 + t * (a2 + t * (a3 + t * (a4 + t * a5))))


def bspline(order: int, i: int, x: float):
    if order == 1:
        return (1.)
    if order == 2:
        if i == 0:
            return (0.5) - x
        if i == 1:
            return (0.5) + x
    if order == 3:
        if i == 0:
            return (0.5) * ((0.5) - x)**2
        if i == 1:
            return (0.75) - (x)**2
        if i == 2:
            return (0.5) * ((0.5) + x)**2
    if order == 4:
        if i == 0:
            return ((1.) + x * ((-6.) + x * ((12.) - x * (8.)))) / (48.)
        if i == 1:
            return ((23.) + x * ((-30.) + x * ((-12.) + x * (24.)))) / (48.)
        if i == 2:
            return ((23.) + x * ((30.) + x * ((-12.) - x * (24.)))) / (48.)
        if i == 3:
            return ((1.) + x * ((6.) + x * ((12.) + x * (8.)))) / (48.)
    if order == 5:
        if i == 0:
            return ((1.) + x * ((-8.) + x * ((24.) + x * ((-32.) + x * (16.))))) / (384.)
        if i == 1:
            return ((19.) + x * ((-44.) + x * ((24.) + x * ((16.) - x * (16.))))) / (96.)
        if i == 2:
            return ((115.) + x * x * ((-120.) + x * x * (48.))) / (192.)
        if i == 3:
            return ((19.) +  x * ((44.) + x * ((24.) + x * ((-16.) - x * (16.))))) / (96.)
        if i == 4:
            return ((1.) + x * ((8.) + x * ((24.) + x * ((32.) + x * (16.))))) / (384.)
    if order == 6:
        if i == 0:
            return ((1.) + x * ((-10.) + x * ((40.) + x * ((-80.) + x * ((80.) - x * (32.)))))) / (3840.)
        if i == 1:
            return ((237.) + x * ((-750.) + x * ((840.) + x * ((-240.) + x * ((-240.) + x * (160.)))))) / (3840.)
        if i == 2:
            return ((841.) + x * ((-770.) + x * ((-440.) + x * ((560.) + x * ((80.) - x * (160.)))))) / (1920.)
        if i == 3:
            return ((841.) + x * (+(770.) + x * ((-440.) + x * ((-560.) + x * ((80.) + x * (160.)))))) / (1920.)
        if i == 4:
            return ((237.) + x * ((750.) + x * ((840.) + x * ((240.) + x * ((-240.) - x * (160.)))))) / (3840.)
        if i == 5:
            return ((1.) + x * ((10.) + x * ((40.) + x * ((80.) + x * ((80.) + x * (32.)))))) / (3840.)
    if order == 7:
        if i == 0:
            return ((1.) + x * ((-12.) + x * ((60.) + x * ((-160.) + x * ((240.) + x * ((-192.) + x * (64.))))))) / (46080.)
        if i == 1:
            return ((361.) + x * ((-1416.) + x * ((2220.) + x * ((-1600.) + x * ((240.) + x * ((384.) - x * (192.))))))) / (23040.)
        if i == 2:
            return ((10543.) + x * ((-17340.) + x * ((4740.) + x * ((6880.) + x * ((-4080.) + x * ((-960.) + x * (960.))))))) / (46080.)
        if i == 3:
            return ((5887.) + x * x * ((-4620.) + x * x * ((1680.) - x * x * (320.)))) / (11520.)
        if i == 4:
            return ((10543.) + x * ((17340.) + x * ((4740.) + x * ((-6880.) + x * ((-4080.) + x * ((960.) + x * (960.))))))) / (46080.)
        if i == 5:
            return ((361.) + x * ((1416.) + x * ((2220.) + x * ((1600.) + x * ((240.) + x * ((-384.) - x * (192.))))))) / (23040.)
        if i == 6:
            return ((1.) + x * ((12.) + x * ((60.) + x * ((160.) + x * ((240.) + x * ((192.) + x * (64.))))))) / (46080.)

    raise RuntimeError("Internal interpolation error.")


def bspline_d(order: int, i: int, x: float):
    if order == 1:
        return 0.
    if order == 2:
        if i == 0:
            return (-1.)
        if i == 1:
            return (1.)
    if order == 3:
        if i == 0:
            return x - (0.5)
        if i == 1:
            return (-2.) * x
        if i == 2:
            return x + (0.5)
    if order == 4:
        if i == 0:
            return ((-1.) + x * ((4.) + x * (-4.))) / (8.)
        if i == 1:
            return ((-5.) + x * ((-4.) + x * (12.))) / (8.)
        if i == 2:
            return ((5.) + x * ((-4.) + x * (-12.))) / (8.)
        if i == 3:
            return ((1.) + x * ((4.) + x * (4.))) / (8.)
    if order == 5:
        if i == 0:
            return ((-1.) + x * ((6.) + x * ((-12.) + x * (8.)))) / (48.)
        if i == 1:
            return ((-11.) + x * ((12.) + x * ((12.) + x * (-16.)))) / (24.)
        if i == 2:
            return (x * ((-5.) + x * x * (4.))) / (4.)
        if i == 3:
            return ((11.) + x * ((12.) + x * ((-12.) + x * (-16.)))) / (24.)
        if i == 4:
            return ((1.) + x * ((6.) + x * ((12.) + x * (8.)))) / (48.)
    if order == 6:
        if i == 0:
          return ((-1.) + x * ((8.) + x * ((-24.) + x * ((32.) + x * (-16.))))) / (384.)
        if i == 1:
          return ((-75.) + x * ((168.) + x * ((-72.) + x * ((-96.) + x * ((80.)))))) / (384.)
        if i == 2:
          return ((-77.) + x * ((-88.) + x * ((168.) + x * ((32.) + x * (-80.))))) / (192.)
        if i == 3:
          return ((77.) + x * ((-88.) + x * ((-168.) + x * ((32.) + x * (80.))))) / (192.)
        if i == 4:
          return ((75.) + x * ((168.) + x * ((72.) + x * ((-96.) + x * (-80.))))) / (384.)
        if i == 5:
          return ((1.) + x * ((8.) + x * ((24.) + x * ((32.) + x * (16.))))) / (384.)
    if order == 7:
        if i == 0:
            return ((-1.) + x * ((10.) + x * ((-40.) + x * ((80.) + x * ((-80.) + x * (32.)))))) / (3840.)
        if i == 1:
            return ((-59.) + x * ((185.) + x * ((-200.) + x * ((40.) + x * ((80.) - x * (48.)))))) / (960.)
        if i == 2:
            return ((-289.) + x * ((158.) + x * ((344.) + x * ((-272.) + x * ((-80.) + x * (96.)))))) / (768.)
        if i == 3:
            return (x * ((-77.) + x * x * ((56.) - x * x * (16.)))) / (96.);
        if i == 4:
            return ((289.) + x * ((158.) + x * ((-344.) + x * ((-272.) + x * ((80.) + x * (96.)))))) / (768.)
        if i == 5:
            return ((59.) + x * ((185.) + x * ((200.) + x * ((40.) + x * ((-80.) - x * (48.)))))) / (960.)
        if i == 6:
            return ((1.) + x * ((10.) + x * ((40.) + x * ((80.) + x * ((80.) + x * (32.)))))) / (3840.)

    raise RuntimeError("Internal interpolation error.")


def multiply_complex_by_imaginary(z, k):
    # Perform the multiplication manually: (re + i*imag) * (i*k)
    return complex(-z.imag * k, z.real * k)


def multiply_complex_by_real(z, k):
    # Perform the multiplication manually: (re + i*imag) * (i*k)
    return complex(z.real * k, z.imag * k)


def complex_norm2(z):
    return z.real**2 + z.imag**2


def get_linear_index(pos, adim, memory_order=COLUMN_MAJOR):
    assert (pos[0] >= 0) and (pos[0] < adim[0])
    assert (pos[1] >= 0) and (pos[1] < adim[1])
    assert (pos[2] >= 0) and (pos[2] < adim[2])
    if memory_order == COLUMN_MAJOR:
        return pos[0] + adim[0] * (pos[1] + adim[1] * pos[2])
    return adim[1] * adim[2] * pos[0] + adim[2] * pos[1] + pos[2]


def calc_meshift(mesh_size, zero_out_midpoint=False):
    ret = []
    for i in range(3):
        length = mesh_size[i]
        shift = np.zeros(mesh_size[i], dtype=int)
        for j in range(1, length // 2 + 1):
            shift[j] = j
            shift[-j] = -j
        if zero_out_midpoint:
            shift[length // 2] = 0
        ret.append(shift)
    return ret


def get_mi_coord(a: float, b: float, box_length: float,
                 box_length_inv: float, box_length_half: float, periodic: bool):
    dx = a - b
    if periodic and (abs(dx) > box_length_half):
        return dx - round(dx * box_length_inv) * box_length
    return dx


def periodic_fold_img(x, i, l):
    while x < 0.:
        x += l
        i -= 1
    while x >= l:
        x -= l
        i += 1
    return (x, i)


def periodic_fold(x, l):
    while x < 0.:
        x += l
    while x >= l:
        x -= l
    return x


def image_shift(image_box, box):
    return image_box * box


def unfolded_position(pos, image_box, box):
    return pos + image_shift(image_box, box)


class BoxGeometry:
    def __init__(self):
        self.m_periodic = np.array([True, True, True])
        self.m_length = np.array([1., 1., 1.])
        self.m_length_inv = np.array([1., 1., 1.])
        self.m_length_half = np.array([0.5, 0.5, 0.5])
        self.set_length(np.array([1., 1., 1.]))
        self.set_periodic(0, True)
        self.set_periodic(1, True)
        self.set_periodic(2, True)

    def set_periodic(self, coord: int, val: bool):
        self.m_periodic[coord] = val

    def periodic(self, coord: int):
        return self.m_periodic[coord]

    def length(self):
        return self.m_length

    def length_inv(self):
        return self.m_length_inv

    def length_half(self):
        return self.m_length_half

    def set_length(self, box_l):
        self.m_length = box_l
        self.m_length_inv = 1. / box_l
        self.m_length_half = 0.5 * box_l

    def volume(self):
        return np.prod(self.m_length)

    def get_mi_coord(self, a, b, coord: int):
        return get_mi_coord(a, b, self.m_length[coord], self.m_length_inv[coord],
                            self.m_length_half[coord], self.m_periodic[coord])

    def get_mi_vector(self, a, b):
        return np.array([self.get_mi_coord(a[0], b[0], 0),
                         self.get_mi_coord(a[1], b[1], 1),
                         self.get_mi_coord(a[2], b[2], 2)])

    def fold_position(self, pos, image_box):
        for i in range(3):
            if self.m_periodic[i]:
                pos[i], image_box[i] = periodic_fold_img(pos[i], image_box[i], self.m_length[i])

    def folded_position(self, pos):
        pos_folded = np.copy(pos)
        for i in range(3):
            if self.m_periodic[i]:
                pos_folded[i] = periodic_fold_img(pos[i], self.m_length[i])
        return pos_folded

    def folded_image_box(self, pos, image_box):
        image_box_folded = np.copy(image_box)
        for i in range(3):
            if self.m_periodic[i]:
                image_box_folded[i] = periodic_fold_img(pos[i], image_box[i], self.m_length[i])[1]
        return image_box_folded

    def image_shift(self, image_box):
        return image_shift(image_box, self.m_length)

    def unfolded_position(self, pos, image_box):
        return unfolded_position(pos, image_box, self.m_length)


class LocalBox:
    def __init__(self, lower_corner=(0, 0, 0), local_box_length=(1, 1, 1),
                 boundaries=(0, 0, 0, 0, 0, 0)):
        self.m_local_box_l = np.array(local_box_length, dtype=float)
        self.m_lower_corner = np.array(lower_corner, dtype=float)
        self.m_upper_corner = self.m_lower_corner + self.m_local_box_l
        self.m_boundaries = np.array(boundaries, dtype=int)

    def my_left(self):
        return self.m_lower_corner

    def my_right(self):
        return self.m_upper_corner

    def length(self):
        return self.m_local_box_l

    def boundary(self):
        return self.m_boundaries

    @classmethod
    def make_regular_decomposition(cls, box_l, node_index, node_grid):
        local_length = box_l / node_grid
        my_left = node_index * local_length
        boundaries = np.zeros(6, dtype=int)
        for i in range(3):
            boundaries[2 * i] = (node_index[i] == 0)
            boundaries[2 * i + 1] = -(node_index[i] + 1 == node_grid[i])
        return LocalBox(my_left, local_length, boundaries)


class P3MParameters:
    def __init__(self, tuning, epsilon: float, accuracy: float,
                 r_cut: float, mesh, mesh_off, cao: int, alpha: float):
        self.tuning = bool(tuning)
        self.epsilon = float(epsilon)
        self.accuracy = float(accuracy)
        self.r_cut = float(r_cut)
        self.r_cut_iL = 0.
        self.alpha = float(alpha)
        self.alpha_L = 0.
        self.cao = int(cao)
        self.cao3 = -1
        self.mesh = np.array(mesh, dtype=int)
        self.mesh_off = np.array(mesh_off, dtype=float)
        self.cao_cut = np.zeros(3, dtype=float)
        self.a = np.zeros(3, dtype=float)
        self.ai = np.zeros(3, dtype=float)
        value_to_tune = -1.
        if epsilon < 0.:
            raise ValueError("Parameter 'epsilon' must be >= 0")
        if accuracy <= 0.:
            raise ValueError("Parameter 'accuracy' must be > 0")
        if r_cut <= 0.:
            if tuning and r_cut == value_to_tune:
                self.r_cut = 0.
            else:
                raise ValueError("Parameter 'r_cut' must be > 0")
        if alpha <= 0.:
            if tuning and alpha == value_to_tune:
                self.alpha = 0.
            else:
                raise ValueError("Parameter 'alpha' must be > 0")
        if (not (np.all(mesh >= 1) or
                 ((mesh[0] >= 1) and (mesh == [mesh[0], -1, -1])) or
                 (tuning and np.all(mesh == 1)))):
            raise ValueError("Parameter 'mesh' must be > 0")
        if (not (np.all(mesh_off >= 0.) and
                 np.all(mesh_off <= 1.))):
            if (np.all(mesh_off == -1.)):
                self.mesh_off = 3 * [P3M_MESHOFF]
            else:
                raise ValueError("Parameter 'mesh_off' must be >= 0 and <= 1")
        if (cao < 1 or cao > 7) and (not tuning or cao != -1):
            raise ValueError("Parameter 'cao' must be >= 1 and <= 7")
        if not tuning and np.all(cao > mesh):
            raise ValueError("Parameter 'cao' cannot be larger than 'mesh'")

    def recalc_a_ai_cao_cut(self, box_l):
        self.ai = self.mesh / box_l
        self.a = 1. / self.ai
        self.cao_cut = (self.cao / 2.) * self.a

    def calc_grid_pos(self, pos):
        return (pos * self. ai) - self.mesh_off


class P3MLocalMesh:
    def __init__(self):
        self.dim = np.zeros(3, dtype=int)
        self.dim_no_halo = np.zeros(3, dtype=int)
        self.size = 0
        self.ld_ind = np.zeros(3, dtype=int)
        self.ld_pos = np.zeros(3, dtype=float)
        self.ld_no_halo = np.zeros(3, dtype=int)
        self.ur_no_halo = np.zeros(3, dtype=int)
        self.inner = np.zeros(3, dtype=int)
        self.in_ld = np.zeros(3, dtype=int)
        self.in_ur = np.zeros(3, dtype=int)
        self.margin = np.zeros(6, dtype=int)
        self.n_halo_ld = np.zeros(3, dtype=int)
        self.n_halo_ur = np.zeros(3, dtype=int)
        self.r_margin = np.zeros(6, dtype=int)
        self.q_2_off = 0
        self.q_21_off = 0

    def recalc_ld_pos(self, params: P3MParameters):
        self.ld_pos = (self.ld_ind + params.mesh_off) * params.a

    def calc_local_ca_mesh(self, params: P3MParameters,
                           local_geo: LocalBox, skin: float,
                           space_layer: float):
        # total skin size
        full_skin = params.cao_cut + 3 * [skin] + [0., 0., space_layer]
        # inner left down corner
        inner_ld_pos = local_geo.my_left()
        # inner up right corner
        inner_ur_pos = local_geo.my_right()
        # outer up right corner
        outer_ur_pos = inner_ur_pos + full_skin
        # outer left down corner
        outer_ld_pos = inner_ld_pos - full_skin
        # convert spatial positions to grid positions
        inner_ld_grid_pos = params.calc_grid_pos(inner_ld_pos)
        inner_ur_grid_pos = params.calc_grid_pos(inner_ur_pos)
        outer_ld_grid_pos = params.calc_grid_pos(outer_ld_pos)
        outer_ur_grid_pos = params.calc_grid_pos(outer_ur_pos)

        # inner left down grid point (global index)
        self.in_ld = np.array(np.ceil(inner_ld_grid_pos), dtype=int)
        # inner up right grid point (global index)
        self.in_ur = np.array(np.floor(inner_ur_grid_pos), dtype=int)

        ROUND_ERROR_PREC = 1e-14
        # correct roundoff errors at boundary
        for i in range(3):
            if (inner_ur_grid_pos[i] - self.in_ur[i] < ROUND_ERROR_PREC):
                self.in_ur[i] -= 1
            if (inner_ld_grid_pos[i] - self.in_ld[i] + 1. < ROUND_ERROR_PREC):
                self.in_ld[i] -= 1

        # inner grid dimensions
        self.inner = self.in_ur - self.in_ld + 1
        # index of left down grid point in global mesh
        self.ld_ind = np.array(np.ceil(outer_ld_grid_pos), dtype=int)
        # left down margin
        for i in range(3):
            self.margin[i * 2] = self.in_ld[i] - self.ld_ind[i]
        # up right grid point
        ind = np.array(np.floor(outer_ur_grid_pos), dtype=int)
        # correct roundoff errors at up right boundary
        for i in range(3):
            if (outer_ur_grid_pos[i] - ind[i] == 0.):
                ind[i] -= 1
        # up right margin
        for i in range(3):
            self.margin[(i * 2) + 1] = ind[i] - self.in_ur[i]

        # grid dimension
        size = 1
        for i in range(3):
            self.dim[i] = ind[i] - self.ld_ind[i] + 1
            size *= self.dim[i]

        # reduce inner grid indices from global to local
        for i in range(3):
            self.in_ld[i] = self.margin[i * 2]
        for i in range(3):
            self.in_ur[i] = self.margin[i * 2] + self.inner[i]

        self.q_2_off = self.dim[2] - params.cao
        self.q_21_off = self.dim[2] * (self.dim[1] - params.cao)

        self.n_halo_ld = np.array([self.margin[0], self.margin[2], self.margin[4]], dtype=int)
        self.n_halo_ur = np.array([self.margin[1], self.margin[3], self.margin[5]], dtype=int)
        self.ld_no_halo = self.ld_ind + self.n_halo_ld
        self.ur_no_halo = self.ld_no_halo + self.dim - self.n_halo_ld - self.n_halo_ur
        self.dim_no_halo = self.ur_no_halo - self.ld_no_halo


class InterpolationWeights:
    def __init__(self, cao: int):
        self.cao = cao
        self.ind = None
        self.w_x = np.zeros(cao, dtype=float)
        self.w_y = np.zeros(cao, dtype=float)
        self.w_z = np.zeros(cao, dtype=float)


class p3m_interpolation_cache:
    def __init__(self, cao=0):
        self.m_cao = cao
        self.ca_frac = []
        self.ca_fmp = []

    def size(self):
        return len(self.ca_fmp)

    def cao(self):
        return self.m_cao

    def store(self, cao: int, w: InterpolationWeights):
        assert cao == self.m_cao
        self.ca_fmp.append(w.ind)
        self.ca_frac.extend(w.w_x)
        self.ca_frac.extend(w.w_y)
        self.ca_frac.extend(w.w_z)

    def load(self, cao: int, i: int):
        assert cao == self.m_cao
        assert i < self.size()
        ret = InterpolationWeights(cao)
        ret.ind = self.ca_fmp[i]
        view = self.ca_frac
        offset = 3 * i * cao
        ret.w_x = view[offset + 0 * cao:offset + 1 * cao]
        ret.w_y = view[offset + 1 * cao:offset + 2 * cao]
        ret.w_z = view[offset + 2 * cao:offset + 3 * cao]
        return ret

    def reset(self, cao: int):
        self.m_cao = cao
        self.ca_frac.clear()
        self.ca_fmp.clear()


def p3m_calculate_interpolation_weights(cao: int, position, ai, local_mesh: P3MLocalMesh):
    pos_shift = math.floor((cao - 1) / 2.0) - (cao % 2) / 2.0
    dist = np.zeros(3, dtype=float)
    nmp = np.zeros(3, dtype=int)
    for i in range(3):
        pos = ((position[i] - local_mesh.ld_pos[i]) * ai[i]) - pos_shift
        nmp[i] = int(pos)
        dist[i] = (pos - nmp[i]) - 0.5
    ret = InterpolationWeights(cao)
    ret.ind = get_linear_index(nmp, local_mesh.dim, ROW_MAJOR)
    assert np.all((nmp + cao) <= local_mesh.dim)
    for i in range(cao):
        ret.w_x[i] = bspline(cao, i, dist[0])
        ret.w_y[i] = bspline(cao, i, dist[1])
        ret.w_z[i] = bspline(cao, i, dist[2])
    return ret


def p3m_interpolate(cao: int, local_mesh: P3MLocalMesh,
                    weights: InterpolationWeights, kernel):
    q_ind = weights.ind
    for i0 in range(cao):
        tmp0 = weights.w_x[i0]
        for i1 in range(cao):
            tmp1 = tmp0 * weights.w_y[i1]
            for i2 in range(cao):
                kernel(q_ind, tmp1 * weights.w_z[i2])
                q_ind += 1
            q_ind += local_mesh.q_2_off
        q_ind += local_mesh.q_21_off


def mesh_sendrecv(sendbuf, scount: int, dest: int, recvbuf,
                  rcount: int, source: int, comm, tag: int):
    raise NotImplementedError("MPI communication not implemented")


def fft_pack_block(inp, out, start, size, dim, element: int):
    assert out.dtype is inp.dtype, f"{out.dtype} is not {inp.dtype}"
    # offsets for indices in input grid
    m_in_offset = element * dim[2]
    s_in_offset = element * (dim[2] * (dim[1] - size[1]))
    # offsets for indices in output grid
    m_out_offset = element * size[2]
    # linear index of input grid, linear index of output grid
    li_in = element * (start[2] + dim[2] * (start[1] + dim[1] * start[0]))
    li_out = 0
    for s in range(size[0]):
        for m in range(size[1]):
            for j in range(element * size[2]):
                out[li_out + j] = inp[li_in + j]
            li_in += m_in_offset
            li_out += m_out_offset
        li_in += s_in_offset


def fft_unpack_block(inp, out, start, size, dim, element: int):
    assert out.dtype is inp.dtype, f"{out.dtype} is not {inp.dtype}"
    # offsets for indices in output grid
    m_out_offset = element * dim[2]
    s_out_offset = element * (dim[2] * (dim[1] - size[1]))
    # offset for indices in input grid
    m_in_offset = element * size[2]
    # linear index of in grid, linear index of out grid
    li_in = 0
    li_out = element * (start[2] + dim[2] * (start[1] + dim[1] * start[0]))
    for s in range(size[0]):
        for m in range(size[1]):
            for j in range(element * size[2]):
                out[li_out + j] = inp[li_in + j]
            li_in += m_in_offset
            li_out += m_out_offset
        li_out += s_out_offset


def p3m_add_block(inp, out, start, size, dim):
    assert out.dtype is inp.dtype, f"{out.dtype} is not {inp.dtype}"
    # fast, mid and slow changing indices
    f = 0
    m = 0
    s = 0
    # linear index of in grid, linear index of out grid
    li_in = 0
    li_out = start[2] + (dim[2] * (start[1] + (dim[1] * start[0])))
    # offsets for indices in output grid
    m_out_offset = dim[2] - size[2]
    s_out_offset = (dim[2] * (dim[1] - size[1]))

    for s in range(size[0]):
        for m in range(size[1]):
            for f in range(size[2]):
                out[li_out] += inp[li_in]
                li_out += 1
                li_in += 1
            li_out += m_out_offset
        li_out += s_out_offset


def extract_block(in_array, dimensions, start, stop):
    block_dim = np.array(stop) - start
    out_array = np.empty(np.prod(block_dim), dtype=in_array.dtype)
    for x in range(block_dim[0]):
        for y in range(block_dim[1]):
            for z in range(block_dim[2]):
                # Compute indices for input and output arrays
                in_index = get_linear_index([x + start[0], y + start[1], z + start[2]], dimensions, ROW_MAJOR)
                out_index = get_linear_index([x, y, z], block_dim, ROW_MAJOR)
                out_array[out_index] = in_array[in_index]
    return out_array


def pad_with_zeros_discard_imag(cropped_array, cropped_dim, pad_left, pad_right):
    # Calculate dimensions and strides
    padded_dim = cropped_dim + np.array(pad_left) + pad_right
    cropped_xy_stride = cropped_dim[1] * cropped_dim[2]
    padded_xy_stride = padded_dim[1] * padded_dim[2]

    # Output vector to hold the padded field (initialized with zeros)
    padded_array = np.zeros(np.prod(padded_dim), dtype=cropped_array[0].real.dtype)

    # Calculate the starting position in the padded array for the inner field
    padded_start_x = pad_left[0] * padded_xy_stride
    padded_start_y = pad_left[1] * padded_dim[2] + pad_left[2]

    # Fill in the original cropped field into the padded array by chunks
    for x in range(cropped_dim[0]):
        cropped_x_offset = x * cropped_xy_stride
        padded_x_offset = padded_start_x + x * padded_xy_stride
        for y in range(cropped_dim[1]):
            cropped_y_offset = cropped_x_offset + y * cropped_dim[2]
            padded_y_offset = padded_x_offset + y * padded_dim[2] + padded_start_y
            # Copy a contiguous slice of the z-dimension at once
            for i in range(cropped_dim[2]):
                padded_array[padded_y_offset+i] = cropped_array[cropped_y_offset+i].real

    return padded_array


class p3m_send_mesh:
    REQ_P3M_INIT = 200
    REQ_P3M_GATHER = 201
    REQ_P3M_SPREAD = 202

    def __init__(self):
        # dimension of sub meshes to send
        self.s_dim = np.zeros((6, 3), dtype=int)
        # lower left corners of sub meshes to send
        self.s_ld = np.zeros((6, 3), dtype=int)
        # upper right corners of sub meshes to send
        self.s_ur = np.zeros((6, 3), dtype=int)
        # sizes for send buffers
        self.s_size = np.zeros(6, dtype=int)
        # dimension of sub meshes to recv
        self.r_dim = np.zeros((6, 3), dtype=int)
        # lower left corners of sub meshes to recv
        self.r_ld = np.zeros((6, 3), dtype=int)
        # upper right corners of sub meshes to recv
        self.r_ur = np.zeros((6, 3), dtype=int)
        # sizes for recv buffers
        self.r_size = np.zeros(6, dtype=int)
        # maximal size for send/recv buffers
        self.max = 0
        # grid points to send
        self.send_grid = np.array([], dtype=float)
        # grid points to recv
        self.recv_grid = np.array([], dtype=float)

    def resize(self, comm, local_mesh: P3MLocalMesh):
        if comm.size() != 1:
            raise NotImplementedError("Cartesian communicator not implemented")
        node_neighbors = 6 * [0]

        done = [0, 0, 0]
        # send grids
        for i in range(3):
            for j in range(3):
                # left
                self.s_ld[i * 2][j] = 0 + done[j] * local_mesh.margin[j * 2]
                if j == i:
                    self.s_ur[i * 2][j] = local_mesh.margin[j * 2]
                else:
                    self.s_ur[i * 2][j] = local_mesh.dim[j] - done[j] * local_mesh.margin[(j * 2) + 1]
                # right
                if j == i:
                    self.s_ld[(i * 2) + 1][j] = local_mesh.in_ur[j]
                else:
                    self.s_ld[(i * 2) + 1][j] = 0 + done[j] * local_mesh.margin[j * 2]
                self.s_ur[(i * 2) + 1][j] = local_mesh.dim[j] - done[j] * local_mesh.margin[(j * 2) + 1]
            done[i] = 1

        self.max = 0
        for i in range(6):
            self.s_size[i] = 1
            for j in range(3):
                self.s_dim[i][j] = self.s_ur[i][j] - self.s_ld[i][j]
                self.s_size[i] *= self.s_dim[i][j]
            self.max = max(self.max, self.s_size[i])

        r_margin = np.zeros(6, dtype=int)
        for i in range(6):
            j = i + (1 if (i % 2 == 0) else - 1)
            if node_neighbors[i] != comm.rank():
                mesh_sendrecv(r_margin.view()[i:i + 1], 1, node_neighbors[i],
                              r_margin.view()[j:j + 1], 1, node_neighbors[j],
                              comm, self.REQ_P3M_INIT)
            else:
                r_margin[j] = local_mesh.margin[i]

        # recv grids
        for i in range(3):
            for j in range(3):
                if j == i:
                    self.r_ld[i * 2][j] = self.s_ld[i * 2][j] + local_mesh.margin[2 * j]
                    self.r_ur[i * 2][j] = self.s_ur[i * 2][j] + r_margin[2 * j]
                    self.r_ld[(i * 2) + 1][j] = self.s_ld[(i * 2) + 1][j] - r_margin[(2 * j) + 1]
                    self.r_ur[(i * 2) + 1][j] = self.s_ur[(i * 2) + 1][j] - local_mesh.margin[(2 * j) + 1]
                else:
                    self.r_ld[i * 2][j] = self.s_ld[i * 2][j]
                    self.r_ur[i * 2][j] = self.s_ur[i * 2][j]
                    self.r_ld[(i * 2) + 1][j] = self.s_ld[(i * 2) + 1][j]
                    self.r_ur[(i * 2) + 1][j] = self.s_ur[(i * 2) + 1][j]

        for i in range(6):
            self.r_size[i] = 1
            for j in range(3):
                self.r_dim[i][j] = self.r_ur[i][j] - self.r_ld[i][j]
                self.r_size[i] *= self.r_dim[i][j]
            self.max = max(self.max, self.r_size[i])

    def gather_grid(self, comm, meshes, dim):
        if comm.size() != 1:
            raise NotImplementedError("Cartesian communicator not implemented")
        node_neighbors = 6 * [0]

        self.send_grid = np.empty(self.max * len(meshes), dtype=float)
        self.recv_grid = np.empty(self.max * len(meshes), dtype=float)

        # direction loop
        for s_dir in range(6):
            r_dir = s_dir + (1 if (s_dir % 2 == 0) else - 1)

            # pack send block
            if self.s_size[s_dir] > 0:
                for i in range(len(meshes)):
                    fft_pack_block(meshes[i], self.send_grid.view()[i * self.s_size[s_dir]:],
                                   self.s_ld[s_dir], self.s_dim[s_dir], dim, 1)

            # communication
            if node_neighbors[s_dir] != comm.rank():
                send_size = len(meshes) * self.s_size[s_dir]
                recv_size = len(meshes) * self.r_size[r_dir]
                mesh_sendrecv(self.send_grid.view(), send_size, node_neighbors[s_dir],
                              self.recv_grid.view(), recv_size, node_neighbors[r_dir],
                              comm, self.REQ_P3M_SPREAD)
            else:
                self.send_grid, self.recv_grid = self.recv_grid, self.send_grid

            # add recv block
            if self.r_size[r_dir] > 0:
                for i in range(len(meshes)):
                    p3m_add_block(self.recv_grid.view()[i * self.r_size[r_dir]:], meshes[i],
                                  self.r_ld[r_dir], self.r_dim[r_dir], dim)

    def spread_grid(self, comm, meshes, dim):
        if comm.size() != 1:
            raise NotImplementedError("Cartesian communicator not implemented")
        node_neighbors = 6 * [0]

        self.send_grid = np.empty(self.max * len(meshes), dtype=float)
        self.recv_grid = np.empty(self.max * len(meshes), dtype=float)

        # direction loop
        for s_dir in range(5, -1, -1):
            r_dir = s_dir + (1 if (s_dir % 2 == 0) else - 1)

            # pack send block
            if self.r_size[r_dir] > 0:
                for i in range(len(meshes)):
                    fft_pack_block(meshes[i], self.send_grid.view()[i * self.r_size[r_dir]:],
                                   self.r_ld[r_dir], self.r_dim[r_dir], dim, 1)

            # communication
            if node_neighbors[r_dir] != comm.rank():
                send_size = len(meshes) * self.r_size[r_dir]
                recv_size = len(meshes) * self.s_size[s_dir]
                mesh_sendrecv(self.send_grid.view(), send_size, node_neighbors[r_dir],
                              self.recv_grid.view(), recv_size, node_neighbors[s_dir],
                              comm, self.REQ_P3M_SPREAD)
            else:
                self.send_grid, self.recv_grid = self.recv_grid, self.send_grid

            # unpack recv block
            if self.s_size[s_dir] > 0:
                for i in range(len(meshes)):
                    fft_unpack_block(self.recv_grid.view()[i * self.s_size[s_dir]:], meshes[i],
                                     self.s_ld[s_dir], self.s_dim[s_dir], dim, 1)


class P3MStateCommon:
    def __init__(self, parameters: P3MParameters):
        self.params = parameters
        self.local_mesh = P3MLocalMesh()
        # Spatial differential operator in k-space.
        self.d_op = [[], [], []]
        # Force optimised influence function (k-space)
        self.g_force = []
        # Energy optimised influence function (k-space)
        self.g_energy = []

    def calc_differential_operator(self):
        self.d_op = calc_meshift(self.params.mesh, True)


class CoulombP3MState(P3MStateCommon):
    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self.sum_qpart = 0
        self.sum_q2 = 0.
        self.square_sum_q = 0.
        self.inter_weights = p3m_interpolation_cache()

        self.rs_charge_density = np.zeros(0, dtype=float)
        self.ks_charge_density = np.zeros(0, dtype=np.complex_)
        self.rs_E_fields = [np.zeros(0, dtype=float) for _ in range(3)]
        self.ks_E_fields_storage = np.zeros(0, dtype=np.complex_)
        self.rs_E_fields_no_halo = np.zeros(0, dtype=np.complex_)
        self.halo_comm = p3m_send_mesh()
        self.fft = None


class AssignCharge:
    def __init__(self, cao):
        self.cao = cao

    def __call__(self, p3m, particles):
        for p in particles:
            if p.q != 0.:
                self.assign_and_cache(p3m, p.q, p.pos, p3m.inter_weights)

    def assign_and_cache(self, p3m, q: float, real_pos,
                         inter_weights: p3m_interpolation_cache):
        w = p3m_calculate_interpolation_weights(self.cao,
            real_pos, p3m.params.ai, p3m.local_mesh)
        inter_weights.store(self.cao, w)
        self.assign(p3m, q, w)

    def assign(self, p3m, q: float, w):
        def kernel(ind: int, w: float):
            p3m.rs_charge_density[ind] += w * q
        p3m_interpolate(self.cao, p3m.local_mesh, w, kernel)


class AssignForces:
    def __init__(self, cao):
        self.cao = cao

    def __call__(self, p3m, force_prefac: float, particles):
        assert self.cao == p3m.inter_weights.cao()
        p_index = 0
        for p in particles:
            if p.q != 0.:
                pref = p.q * force_prefac
                w = p3m.inter_weights.load(self.cao, p_index)
                force = np.zeros(3, dtype=float)

                def kernel(ind: int, w: float):
                    force[0] += w * p3m.rs_E_fields[0][ind]
                    force[1] += w * p3m.rs_E_fields[1][ind]
                    force[2] += w * p3m.rs_E_fields[2][ind]
                p3m_interpolate(self.cao, p3m.local_mesh, w, kernel)
                p.f -= pref * force
                p_index += 1


class FFTBox:
    def __init__(self, low, high, layout):
        self.low = low
        self.high = high
        self.layout = layout


class P3MFFT:
    def __init__(self, comm, global_mesh_size, rs_local_ld_index,
                 rs_local_ur_index, mem_layout):
        self.comm = comm
        self.memory_layout = np.array(mem_layout, dtype=int)
        self.in_box = FFTBox(rs_local_ld_index, rs_local_ur_index - 1, self.memory_layout)
        self.out_box = FFTBox(rs_local_ld_index, rs_local_ur_index - 1, self.memory_layout)
        self.global_mesh = np.array(global_mesh_size, dtype=int)
        self.init_fft()

    def set_preferred_kspace_decomposition(self):
        n = comm_cart.size()
        proc_grid = np.array([1, 1, 1], dtype=int)
        local_size = np.array(self.global_mesh, dtype=int)
        for d in [2, 1, 0]:
            while (n % 2 == 0 and local_size[d] % 2 == 0):
                proc_grid[d] *= 2
                n /= 2
                local_size[d] /= 2
        if proc_grid[1] == 1 and proc_grid[0] == 1 and proc_grid[2] % 4 == 0 and local_size[1] % 2 == 0:
            proc_grid[2] /= 2
            proc_grid[1] *= 2
        self.init_fft()

    def init_fft(self):
        pass

    def ks_local_ld_index(self):
        return np.copy(self.out_box.low)

    def ks_local_ur_index(self):
        return np.copy(self.out_box.high) + 1

    def ks_local_size(self):
        return self.ks_local_ur_index() - self.ks_local_ld_index()

    def forward(self, inp):
        return np.fft.fftn(inp.reshape(self.global_mesh),
                           norm="backward").reshape([-1])

    def backward(self, inp):
        return np.fft.ifftn(inp.reshape(self.global_mesh),
                            norm="forward").reshape([-1])

    def backward_batch(self, n: int, inps, outs):
        fft_mesh_length = np.prod(self.global_mesh)
        for d in range(n):
            beg = d * fft_mesh_length
            end = (d + 1) * fft_mesh_length
            outs.view()[beg:end] = self.backward(inps.view()[beg:end])

    def mem_layout(self):
        return self.memory_layout


def for_each_3d(start, stop, counters, kernel, projector=lambda x, y: None):
    for nx in range(start[0], stop[0]):
        counters[0] = nx
        projector(0, nx)
        for ny in range(start[1], stop[1]):
            counters[1] = ny
            projector(1, ny)
            for nz in range(start[2], stop[2]):
                counters[2] = nz
                projector(2, nz)
                kernel()


def sinc(x: float):
    epsilon = 0.1
    pix = math.pi * x

    if abs(x) > epsilon:
        return math.sin(pix) / pix

    # Coefficients of the Taylor expansion of sinc
    c0 = +1. / math.factorial(1)
    c2 = -1. / math.factorial(3)
    c4 = +1. / math.factorial(5)
    c6 = -1. / math.factorial(7)
    c8 = +1. / math.factorial(9)

    pix2 = pix * pix
    return c0 + pix2 * (c2 + pix2 * (c4 + pix2 * (c6 + pix2 * c8)))


def G_opt(S: int, m: int, cao: int, alpha: float, k, h):
    k2 = sum(x**2 for x in k)
    if k2 == 0.:
        return 0.

    limit = 30.
    m_start = 3 * [-m]
    m_stop = 3 * [m + 1]
    exponent_prefactor = 1. / (2. * alpha)**2
    wavevector = (2. * math.pi) / h
    wavevector_i = 1. / wavevector
    indices = np.zeros(3, dtype=int)
    km = np.zeros(3, dtype=float)
    fnm = np.zeros(3, dtype=float)
    division = [0., 0.]

    def kernel():
        U2 = np.prod(fnm)**(2 * cao)
        km2 = np.linalg.norm(km)**2
        exponent = exponent_prefactor * km2
        if exponent < limit:
            f3 = math.exp(-exponent) * (4. * math.pi / km2)
            division[0] += U2 * f3 * np.dot(k, km)**S
        division[1] += U2

    def projector(dim: int, n: int):
        km[dim] = k[dim] + n * wavevector[dim]
        fnm[dim] = sinc(km[dim] * wavevector_i[dim])

    for_each_3d(m_start, m_stop, indices, kernel, projector)

    numerator, denominator = division
    if numerator == 0. and k2 != 0. and denominator != 0.:
        return 0.
    return numerator / ((k2**S) * denominator**2)


def grid_influence_function(S: int, params: P3MParameters,
                            n_start, n_stop, KX: int, KY: int, KZ: int, inv_box_l, m=0):
    shifts = calc_meshift(params.mesh)
    size = n_stop - n_start
    g = np.zeros(np.prod(size), dtype=float)
    wavevector = (2. * np.pi) * inv_box_l
    half_mesh = params.mesh // 2
    indices = np.zeros(3, dtype=int)
    index = [0]

    def kernel():
        if ((indices[KX] % half_mesh[0] != 0) or
            (indices[KY] % half_mesh[1] != 0) or
            (indices[KZ] % half_mesh[2] != 0)):
            k = [shifts[0][indices[KX]] * wavevector[0],
                 shifts[1][indices[KY]] * wavevector[1],
                 shifts[2][indices[KZ]] * wavevector[2]]
            g[index[0]] = G_opt(S, m, params.cao, params.alpha, k, params.a)
        index[0] += 1

    for_each_3d(n_start, n_stop, indices, kernel)
    return g


class CoulombP3M:
    def __init__(self, p3m_params: P3MParameters):
        self.p3m_params = p3m_params
        self.prefactor = 0.

    def set_prefactor(self, prefactor: float):
        if prefactor <= 0.:
            raise ValueError("Parameter 'prefactor' must be > 0")
        self.prefactor = prefactor

    def on_boxl_change(self):
        self.scaleby_box_l()

    def on_node_grid_change(self):
        self.sanity_checks_node_grid()

    def on_periodicity_change(self):
        self.sanity_checks_periodicity()

    def on_cell_structure_change(self):
        self.sanity_checks_cell_structure()
        self.init()

    def sanity_checks(self):
        self.sanity_checks_boxl()
        self.sanity_checks_node_grid()
        self.sanity_checks_periodicity()
        self.sanity_checks_cell_structure()
        self.sanity_checks_charge_neutrality()

    def pair_force(self, q1q2: float, d, dist: float):
        if (q1q2 == 0.) or dist >= self.p3m_params.r_cut or dist <= 0.:
            return np.zeros(3, dtype=float)
        alpha = self.p3m_params.alpha
        adist = alpha * dist
        exp_adist_sq = math.exp(-adist * adist)
        dist_sq = dist * dist
        two_a_sqrt_pi_i = 2. * alpha / math.sqrt(math.pi)
        if USE_ERFC_APPROXIMATION:
            erfc_part_ri = AS_erfc_part(adist) / dist
            fac = exp_adist_sq * (erfc_part_ri + two_a_sqrt_pi_i) / dist_sq
        else:
            erfc_part_ri = math.erfc(adist) / dist
            fac = (erfc_part_ri + two_a_sqrt_pi_i * exp_adist_sq) / dist_sq
        return (fac * self.prefactor * q1q2) * d

    def pair_energy(self, q1q2: float, dist: float):
        if (q1q2 == 0.) or dist >= self.p3m_params.r_cut or dist <= 0.:
            return np.zeros(3, dtype=float)
        adist = self.p3m_params.alpha * dist
        if USE_ERFC_APPROXIMATION:
            erfc_part_ri = AS_erfc_part(adist) / dist
            return self.prefactor * q1q2 * erfc_part_ri * math.exp(-adist * adist)
        else:
            erfc_part_ri = math.erfc(adist) / dist
            return self.prefactor * q1q2 * erfc_part_ri

    def sanity_checks_boxl(self):
        system = get_system()
        box_geo = system.box_geo
        local_geo = system.local_geo
        for i in range(3):
            if self.p3m_params.cao_cut[i] >= box_geo.length_half()[i]:
                raise RuntimeError(f"P3M_init: k-space cutoff {self.p3m_params.cao_cut[i]} is larger than half of box dimension {box_geo.length()[i]}")
            if self.p3m_params.cao_cut[i] >= local_geo.length()[i]:
                raise RuntimeError(f"P3M_init: k-space cutoff {self.p3m_params.cao_cut[i]} is larger than local box dimension {local_geo.length()[i]}")
        if self.p3m_params.epsilon != P3M_EPSILON_METALLIC:
            if ((box_geo.length()[0] != box_geo.length()[1]) or
                (box_geo.length()[1] != box_geo.length()[2]) or
                (self.p3m_params.mesh[0] != self.p3m_params.mesh[1]) or
                (self.p3m_params.mesh[1] != self.p3m_params.mesh[2])):
                raise RuntimeError("CoulombP3M: non-metallic epsilon requires cubic box")

    def sanity_checks_node_grid(self):
        node_grid = communicator.node_grid
        if (node_grid[0] < node_grid[1] or node_grid[1] < node_grid[2]):
            raise RuntimeError("CoulombP3M: node grid must be sorted, largest first")

    def sanity_checks_periodicity(self):
        box_geo = get_system().box_geo;
        if (not box_geo.periodic(0) or not box_geo.periodic(1) or not box_geo.periodic(2)):
            raise RuntimeError("CoulombP3M: requires periodicity (True, True, True)")

    def sanity_checks_cell_structure(self):
        pass


class CoulombP3MImpl(CoulombP3M):
    def mem_layout(self):
        return np.array([2, 1, 0], dtype=int)  # column major

    def __init__(self, p3m_state: CoulombP3MState, prefactor: float):
        super().__init__(p3m_state.params)
        self.p3m = p3m_state
        self.p3m_state_ptr = p3m_state
        self.tune_timings = 10
        self.tune_timings = True
        self.check_complex_residuals = True
        if self.tune_timings <= 0:
            raise ValueError("Parameter 'timings' must be > 0")
        self.m_is_tuned = not self.p3m.params.tuning
        self.p3m.params.tuning = False
        self.set_prefactor(prefactor)
        if self.p3m.params.alpha_L == 0. and self.p3m.params.alpha != 0.:
            self.p3m.params.alpha_L = self.p3m.params.alpha * get_system().box_geo.length()[0]
        if self.p3m.params.r_cut_iL == 0. and self.p3m.params.r_cut != 0.:
            self.p3m.params.r_cut_iL = self.p3m.params.r_cut * get_system().box_geo.length_inv()[0]

    def init(self):
        self.init_cpu_kernels()

    def tune(self):
        system = get_system()
        box_geo = system.box_geo
        if self.p3m.params.alpha_L == 0. and self.p3m.params.alpha != 0.:
            self.p3m.params.alpha_L = self.p3m.params.alpha * box_geo.length()[0]
        if self.p3m.params.r_cut_iL == 0. and self.p3m.params.r_cut != 0.:
            self.p3m.params.r_cut_iL = self.p3m.params.r_cut * box_geo.length_inv()[0]
        if not self.is_tuned():
            self.count_charged_particles()
            if self.p3m.sum_qpart == 0:
                raise RuntimeError("CoulombP3M: no charged particles in the system")
            raise NotImplementedError("CoulombTuningAlgorithm not implemented")
        self.init()

    def count_charged_particles(self):
        local_n = 0
        local_q2 = 0.
        local_q = 0.
        for p in get_system().particles:
            if p.q != 0.:
                local_n += 1
                local_q2 += p.q * p.q
                local_q += p.q
        # boost::mpi::all_reduce(comm_cart, local_n, p3m.sum_qpart, std::plus<>());
        # boost::mpi::all_reduce(comm_cart, local_q2, p3m.sum_q2, std::plus<>());
        # boost::mpi::all_reduce(comm_cart, local_q, p3m.square_sum_q, std::plus<>());
        self.p3m.sum_qpart = local_n
        self.p3m.sum_q2 = local_q2
        self.p3m.square_sum_q = local_q
        self.p3m.square_sum_q = self.p3m.square_sum_q

    def count_charged_particles_elc(self, n: int, sum_q2: float, square_sum_q: float):
        self.p3m.sum_qpart = n
        self.p3m.sum_q2 = sum_q2
        self.p3m.square_sum_q = square_sum_q

    def adapt_epsilon_elc(self):
        self.p3m.params.epsilon = P3M_EPSILON_METALLIC

    def is_tuned(self):
        return self.m_is_tuned

    def is_gpu(self):
        return False

    def is_double_precision(self):
        return True

    def on_activation(self):
        self.sanity_checks()
        self.tune()

    def long_range_energy(self, particles):
        return self.long_range_kernel(False, True, particles)

    def add_long_range_forces(self, particles):
        self.long_range_kernel(True, False, particles)

    def long_range_pressure(self, particles):
        return None

    def charge_assign(self, particles):
        self.prepare_fft_mesh(True)
        AssignCharge(self.p3m.params.cao)(self.p3m, particles)

    def assign_charge(self, q: float, real_pos, skip_cache: bool):
        if skip_cache:
            AssignCharge(self.p3m.params.cao).assign(self.p3m, q, real_pos)
        else:
            AssignCharge(self.p3m.params.cao).assign_and_cache(
                self.p3m, q, real_pos, self.p3m.inter_weights)

    def prepare_fft_mesh(self, reset_weights: bool):
        if reset_weights:
            self.p3m.inter_weights.reset(self.p3m.params.cao)
        self.p3m.rs_charge_density = np.zeros(np.prod(self.p3m.local_mesh.dim), dtype=float)

    def long_range_kernel(self, force_flag: bool, energy_flag: bool, particles):
        system = get_system()
        box_geo = system.box_geo
        p3m = self.p3m
        if self.p3m.sum_q2 > 0.:
            self.charge_assign(particles)
        p3m.halo_comm.gather_grid(comm_cart, [p3m.rs_charge_density.view()], p3m.local_mesh.dim)
        charge_density_no_halos = extract_block(p3m.rs_charge_density, p3m.local_mesh.dim, p3m.local_mesh.n_halo_ld,
                                                p3m.local_mesh.dim - p3m.local_mesh.n_halo_ur)
        p3m.ks_charge_density = p3m.fft.forward(charge_density_no_halos)
        box_dipole = None
        if p3m.params.epsilon != P3M_EPSILON_METALLIC:
            box_dipole = calc_dipole_moment(comm_cart, particles)
        volume = box_geo.volume()
        pref = 4. * np.pi / volume / (2. * p3m.params.epsilon + 1.)
        energy = 0.

        if energy_flag:
            node_energy = 0.
            fft_mesh_length = np.prod(p3m.fft.ks_local_size())
            for i in range(fft_mesh_length):
                node_energy += p3m.g_energy[i] * complex_norm2(p3m.ks_charge_density[i])
            node_energy /= 2. * volume;

            if comm_cart.rank() == 0:
                energy -= p3m.sum_q2 * p3m.params.alpha / math.sqrt(math.pi)
                energy -= p3m.square_sum_q * math.pi / (2. * volume * p3m.params.alpha**2)
                if box_dipole is not None:
                    energy += pref * np.linagl.norm(box_dipole)
            energy *= self.prefactor

        if force_flag:
            mesh_start = np.zeros(3, dtype=int)
            mesh_stop = p3m.fft.ks_local_size()
            wavevector = 2. * math.pi * box_geo.length_inv()
            indices = np.zeros(3, dtype=int)
            index = [0]
            fft_mesh_length = np.prod(p3m.fft.ks_local_size())
            # Holds the electric field in k-space
            ks_E_fields = []
            for d in range(3):
                ks_E_fields.append(p3m.ks_E_fields_storage.view()[d*fft_mesh_length:(d+1)*fft_mesh_length])

            # compute electric field
            for i in range(len(p3m.ks_charge_density)):
                p3m.ks_charge_density[i] = multiply_complex_by_real(
                    p3m.ks_charge_density[i], p3m.g_force[i])

            def kernel():
                phi_hat = p3m.ks_charge_density[index[0]]
                global_index = indices + p3m.fft.ks_local_ld_index()
                for d in range(3):
                    k = p3m.d_op[d][global_index[d]] * wavevector[d]
                    ks_E_fields[d][index[0]] = multiply_complex_by_imaginary(phi_hat, k)
                index[0] += 1

            for_each_3d(mesh_start, mesh_stop, indices, kernel)

            # Back-transform the k-space electric field to real space
            rs_mesh_length_no_halo = np.prod(p3m.local_mesh.dim_no_halo)
            p3m.fft.backward_batch(3, p3m.ks_E_fields_storage, p3m.rs_E_fields_no_halo)

            # add zeros around the E-field in real space to make room for the ghost layers
            for d in range(3):
                start = d * rs_mesh_length_no_halo
                stop = (d + 1) * rs_mesh_length_no_halo
                p3m.rs_E_fields[d] = pad_with_zeros_discard_imag(
                    p3m.rs_E_fields_no_halo.view()[start:stop],
                    p3m.local_mesh.dim_no_halo,
                    p3m.local_mesh.n_halo_ld,
                    p3m.local_mesh.n_halo_ur)

            # Ghost communicate the boundary layers of the E-field in real space
            field_pointers = [p3m.rs_E_fields[0], p3m.rs_E_fields[1], p3m.rs_E_fields[2]]
            p3m.halo_comm.spread_grid(comm_cart, field_pointers, p3m.local_mesh.dim)

            # Assign particle forces
            force_prefac = self.prefactor / volume
            AssignForces(p3m.params.cao)(p3m, force_prefac, particles)

            # add dipole forces
            if box_dipole is not None:
                dm = self.prefactor * pref * box_dipole
                for p in particles:
                    p.f -= p.q * dm

        return energy

    def calc_influence_function_force(self):
        kxyz = self.p3m.fft.mem_layout()
        self.p3m.g_force = grid_influence_function(
            1, self.p3m.params, self.p3m.fft.ks_local_ld_index(), self.p3m.fft.ks_local_ur_index(),
            kxyz[0], kxyz[1], kxyz[2], get_system().box_geo.length_inv())

    def calc_influence_function_energy(self):
        kxyz = self.p3m.fft.mem_layout()
        self.p3m.g_energy = grid_influence_function(
            0, self.p3m.params, self.p3m.fft.ks_local_ld_index(), self.p3m.fft.ks_local_ur_index(),
            kxyz[0], kxyz[1], kxyz[2], get_system().box_geo.length_inv())

    def scaleby_box_l(self):
        box_geo = get_system().box_geo
        self.p3m.params.r_cut = self.p3m.params.r_cut_iL * box_geo.length()[0]
        self.p3m.params.alpha = self.p3m.params.alpha_L * box_geo.length_inv()[0]
        self.p3m.params.recalc_a_ai_cao_cut(box_geo.length())
        self.p3m.local_mesh.recalc_ld_pos(self.p3m.params)
        self.sanity_checks_boxl()
        self.calc_influence_function_force()
        self.calc_influence_function_energy()
        self.p3m.halo_comm.resize(comm_cart, self.p3m.local_mesh)

    def init_cpu_kernels(self):
        p3m = self.p3m
        assert np.all(p3m.params.mesh >= 1)
        assert p3m.params.cao >= 1 and p3m.params.cao <= 7
        assert p3m.params.alpha > 0.
        system = get_system()
        box_geo = system.box_geo
        local_geo = system.local_geo
        skin = system.skin
        self.p3m.params.cao3 = self.p3m.params.cao**3
        self.p3m.params.recalc_a_ai_cao_cut(box_geo.length())
        elc_layer = 0.
        self.p3m.local_mesh.calc_local_ca_mesh(
            self.p3m.params, local_geo, skin, elc_layer)
        self.p3m.fft = P3MFFT(
            comm_cart,
            self.p3m.params.mesh,
            self.p3m.local_mesh.ld_no_halo,
            self.p3m.local_mesh.ur_no_halo,
            self.mem_layout())
        self.p3m.fft.set_preferred_kspace_decomposition()
        rs_array_length = np.prod(self.p3m.local_mesh.dim)
        rs_array_length_no_halo = np.prod(self.p3m.local_mesh.dim_no_halo)
        fft_mesh_length = np.prod(self.p3m.fft.ks_local_size())
        self.p3m.rs_charge_density = np.zeros(rs_array_length, dtype=float)
        self.p3m.ks_charge_density = np.zeros(fft_mesh_length, dtype=np.complex_)
        for i in range(3):
            self.p3m.rs_E_fields[i] = np.zeros(rs_array_length, dtype=float)
        self.p3m.ks_E_fields_storage = np.zeros(3 * fft_mesh_length, dtype=np.complex_)
        self.p3m.rs_E_fields_no_halo = np.zeros(3 * rs_array_length_no_halo, dtype=np.complex_)
        self.p3m.calc_differential_operator()
        self.scaleby_box_l()
        self.count_charged_particles()


class P3M:
    """
    prefactor : :obj:`float`
        Electrostatics prefactor (see :eq:`coulomb_prefactor`).
    accuracy : :obj:`float`
        P3M tunes its parameters to provide this target accuracy.
    alpha : :obj:`float`, optional
        The Ewald parameter.
    cao : :obj:`float`, optional
        The charge-assignment order, an integer between 1 and 7.
    mesh : :obj:`int` or (3,) array_like of :obj:`int`, optional
        The number of mesh points in x, y and z direction. Use a single
        value for cubic boxes.
    r_cut : :obj:`float`, optional
        The real space cutoff.
    """

    def __init__(self, prefactor: float, accuracy: float,
                 r_cut: float, mesh, cao: int, alpha: float, tune: bool):
        self.prefactor = float(prefactor)
        self.accuracy = float(accuracy)
        self.r_cut = float(r_cut)
        self.alpha = float(alpha)
        self.cao = int(cao)
        self.mesh = np.array(mesh, dtype=int)
        self.epsilon = -1.
        self.params = P3MParameters(
            tuning=tune, epsilon=P3M_EPSILON_METALLIC, r_cut=r_cut, cao=cao, alpha=alpha,
            accuracy=accuracy, mesh=np.array(mesh, dtype=int))
        self.state = CoulombP3MState(self.params)
        self.actor = CoulombP3MImpl(p3m_state=self.state, prefactor=prefactor)

    def activate(self):
        self.actor.init()


class Particle:
    def __init__(self, pos, q):
        self.pos = np.array(pos, dtype=float)
        self.q = float(q)
        self.f = np.zeros(3, dtype=float)


class System:
    def __init__(self, box_l):
        self.box_geo = BoxGeometry()
        self.box_geo.set_length(np.array(box_l, dtype=float))
        self.box_l = self.box_geo.m_length
        self.local_geo = LocalBox(local_box_length=self.box_l)
        self.skin = 0.05 * np.min(self.box_l)
        self.particles = []
        self.solver = None


def short_range_loop(system):
    if system.solver is not None:
        for i, p1 in enumerate(system.particles[:-1]):
            for p2 in system.particles[i + 1:]:
                q1q2 = p1.q * p2.q
                if q1q2 != 0.:
                    vec21 = system.box_geo.get_mi_vector(p1.pos, p2.pos)
                    force = system.solver.actor.pair_force(
                        q1q2, vec21, np.linalg.norm(vec21))
                    p1.f += force
                    p2.f -= force


system = System(box_l=[8, 8, 8])


def get_system():
    return system


p1 = Particle(pos=[0, 0, 0], q=+1)
p2 = Particle(pos=[0, 1, 0], q=-1)
ref_forces = [
    [1.87803086e-18, 9.87485249e-01, 2.43971959e-18],
    [-9.89493301e-19, -9.87485249e-01, 2.52230530e-18],
]
system.particles.append(p1)
system.particles.append(p2)
solver = P3M(prefactor=1, accuracy=1e-3, cao=4,
             alpha=5.2859 / np.max(system.box_l),
             r_cut=0.405176 * np.max(system.box_l),
             mesh=3 * [8], tune=False)
system.solver = solver
solver.activate()
energy = solver.actor.long_range_energy(system.particles)
print(energy)
solver.actor.add_long_range_forces(system.particles)
short_range_loop(system)
for p in system.particles:
    print(p.f)

for i, p in enumerate(system.particles):
    np.testing.assert_allclose(p.f, ref_forces[i], rtol=1e-12, atol=1e-9)
