// kernel generated with pystencils v1.0, lbmpy v1.0,
// lbmpy_walberla/pystencils_walberla from commit
// 01a28162ae1aacf7b96152c9f886ce54cc7f53ff

#include "DensityPackInfo_single_precision.h"
#include "core/DataTypes.h"
#include "core/cell/CellInterval.h"
#include "stencil/Directions.h"

#if (defined WALBERLA_CXX_COMPILER_IS_GNU) ||                                  \
    (defined WALBERLA_CXX_COMPILER_IS_CLANG)
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wfloat-equal"
#pragma GCC diagnostic ignored "-Wshadow"
#pragma GCC diagnostic ignored "-Wconversion"
#pragma GCC diagnostic ignored "-Wunused-variable"
#endif

namespace walberla {
namespace pystencils {

using walberla::cell::CellInterval;
using walberla::stencil::Direction;

namespace internal_pack_BSW {
static FUNC_PREFIX void
pack_BSW(float *RESTRICT _data_buffer, float *RESTRICT const _data_j,
         int64_t const _size_j_0, int64_t const _size_j_1,
         int64_t const _size_j_2, int64_t const _stride_j_0,
         int64_t const _stride_j_1, int64_t const _stride_j_2,
         int64_t const _stride_j_3) {
  for (int64_t ctr_2 = 0; ctr_2 < _size_j_2; ctr_2 += 1) {
    float *RESTRICT _data_j_20_39 =
        _data_j + _stride_j_2 * ctr_2 + 9 * _stride_j_3;
    for (int64_t ctr_1 = 0; ctr_1 < _size_j_1; ctr_1 += 1) {
      float *RESTRICT _data_j_20_39_10 = _stride_j_1 * ctr_1 + _data_j_20_39;
      for (int64_t ctr_0 = 0; ctr_0 < _size_j_0; ctr_0 += 1) {
        _data_buffer[_size_j_0 * _size_j_1 * ctr_2 + _size_j_0 * ctr_1 +
                     ctr_0] = _data_j_20_39_10[_stride_j_0 * ctr_0];
      }
    }
  }
}
} // namespace internal_pack_BSW

namespace internal_pack_SW {
static FUNC_PREFIX void
pack_SW(float *RESTRICT _data_buffer, float *RESTRICT const _data_j,
        int64_t const _size_j_0, int64_t const _size_j_1,
        int64_t const _size_j_2, int64_t const _stride_j_0,
        int64_t const _stride_j_1, int64_t const _stride_j_2,
        int64_t const _stride_j_3) {
  for (int64_t ctr_2 = 0; ctr_2 < _size_j_2; ctr_2 += 1) {
    float *RESTRICT _data_j_20_39 =
        _data_j + _stride_j_2 * ctr_2 + 9 * _stride_j_3;
    float *RESTRICT _data_j_20_33 =
        _data_j + _stride_j_2 * ctr_2 + 3 * _stride_j_3;
    float *RESTRICT _data_j_20_310 =
        _data_j + _stride_j_2 * ctr_2 + 10 * _stride_j_3;
    for (int64_t ctr_1 = 0; ctr_1 < _size_j_1; ctr_1 += 1) {
      float *RESTRICT _data_j_20_39_10 = _stride_j_1 * ctr_1 + _data_j_20_39;
      float *RESTRICT _data_j_20_33_10 = _stride_j_1 * ctr_1 + _data_j_20_33;
      float *RESTRICT _data_j_20_310_10 = _stride_j_1 * ctr_1 + _data_j_20_310;
      for (int64_t ctr_0 = 0; ctr_0 < _size_j_0; ctr_0 += 1) {
        _data_buffer[3 * _size_j_0 * _size_j_1 * ctr_2 + 3 * _size_j_0 * ctr_1 +
                     3 * ctr_0] = _data_j_20_39_10[_stride_j_0 * ctr_0];
        _data_buffer[3 * _size_j_0 * _size_j_1 * ctr_2 + 3 * _size_j_0 * ctr_1 +
                     3 * ctr_0 + 1] = _data_j_20_33_10[_stride_j_0 * ctr_0];
        _data_buffer[3 * _size_j_0 * _size_j_1 * ctr_2 + 3 * _size_j_0 * ctr_1 +
                     3 * ctr_0 + 2] = _data_j_20_310_10[_stride_j_0 * ctr_0];
      }
    }
  }
}
} // namespace internal_pack_SW

namespace internal_pack_TSW {
static FUNC_PREFIX void
pack_TSW(float *RESTRICT _data_buffer, float *RESTRICT const _data_j,
         int64_t const _size_j_0, int64_t const _size_j_1,
         int64_t const _size_j_2, int64_t const _stride_j_0,
         int64_t const _stride_j_1, int64_t const _stride_j_2,
         int64_t const _stride_j_3) {
  for (int64_t ctr_2 = 0; ctr_2 < _size_j_2; ctr_2 += 1) {
    float *RESTRICT _data_j_20_310 =
        _data_j + _stride_j_2 * ctr_2 + 10 * _stride_j_3;
    for (int64_t ctr_1 = 0; ctr_1 < _size_j_1; ctr_1 += 1) {
      float *RESTRICT _data_j_20_310_10 = _stride_j_1 * ctr_1 + _data_j_20_310;
      for (int64_t ctr_0 = 0; ctr_0 < _size_j_0; ctr_0 += 1) {
        _data_buffer[_size_j_0 * _size_j_1 * ctr_2 + _size_j_0 * ctr_1 +
                     ctr_0] = _data_j_20_310_10[_stride_j_0 * ctr_0];
      }
    }
  }
}
} // namespace internal_pack_TSW

namespace internal_pack_BW {
static FUNC_PREFIX void
pack_BW(float *RESTRICT _data_buffer, float *RESTRICT const _data_j,
        int64_t const _size_j_0, int64_t const _size_j_1,
        int64_t const _size_j_2, int64_t const _stride_j_0,
        int64_t const _stride_j_1, int64_t const _stride_j_2,
        int64_t const _stride_j_3) {
  for (int64_t ctr_2 = 0; ctr_2 < _size_j_2; ctr_2 += 1) {
    float *RESTRICT _data_j_20_39 =
        _data_j + _stride_j_2 * ctr_2 + 9 * _stride_j_3;
    float *RESTRICT _data_j_20_35 =
        _data_j + _stride_j_2 * ctr_2 + 5 * _stride_j_3;
    float *RESTRICT _data_j_20_311 =
        _data_j + _stride_j_2 * ctr_2 + 11 * _stride_j_3;
    for (int64_t ctr_1 = 0; ctr_1 < _size_j_1; ctr_1 += 1) {
      float *RESTRICT _data_j_20_39_10 = _stride_j_1 * ctr_1 + _data_j_20_39;
      float *RESTRICT _data_j_20_35_10 = _stride_j_1 * ctr_1 + _data_j_20_35;
      float *RESTRICT _data_j_20_311_10 = _stride_j_1 * ctr_1 + _data_j_20_311;
      for (int64_t ctr_0 = 0; ctr_0 < _size_j_0; ctr_0 += 1) {
        _data_buffer[3 * _size_j_0 * _size_j_1 * ctr_2 + 3 * _size_j_0 * ctr_1 +
                     3 * ctr_0] = _data_j_20_39_10[_stride_j_0 * ctr_0];
        _data_buffer[3 * _size_j_0 * _size_j_1 * ctr_2 + 3 * _size_j_0 * ctr_1 +
                     3 * ctr_0 + 1] = _data_j_20_35_10[_stride_j_0 * ctr_0];
        _data_buffer[3 * _size_j_0 * _size_j_1 * ctr_2 + 3 * _size_j_0 * ctr_1 +
                     3 * ctr_0 + 2] = _data_j_20_311_10[_stride_j_0 * ctr_0];
      }
    }
  }
}
} // namespace internal_pack_BW

namespace internal_pack_W {
static FUNC_PREFIX void
pack_W(float *RESTRICT _data_buffer, float *RESTRICT const _data_j,
       int64_t const _size_j_0, int64_t const _size_j_1,
       int64_t const _size_j_2, int64_t const _stride_j_0,
       int64_t const _stride_j_1, int64_t const _stride_j_2,
       int64_t const _stride_j_3) {
  for (int64_t ctr_2 = 0; ctr_2 < _size_j_2; ctr_2 += 1) {
    float *RESTRICT _data_j_20_39 =
        _data_j + _stride_j_2 * ctr_2 + 9 * _stride_j_3;
    float *RESTRICT _data_j_20_33 =
        _data_j + _stride_j_2 * ctr_2 + 3 * _stride_j_3;
    float *RESTRICT _data_j_20_310 =
        _data_j + _stride_j_2 * ctr_2 + 10 * _stride_j_3;
    float *RESTRICT _data_j_20_35 =
        _data_j + _stride_j_2 * ctr_2 + 5 * _stride_j_3;
    float *RESTRICT _data_j_20_30 = _data_j + _stride_j_2 * ctr_2;
    float *RESTRICT _data_j_20_36 =
        _data_j + _stride_j_2 * ctr_2 + 6 * _stride_j_3;
    float *RESTRICT _data_j_20_311 =
        _data_j + _stride_j_2 * ctr_2 + 11 * _stride_j_3;
    float *RESTRICT _data_j_20_34 =
        _data_j + _stride_j_2 * ctr_2 + 4 * _stride_j_3;
    float *RESTRICT _data_j_20_312 =
        _data_j + _stride_j_2 * ctr_2 + 12 * _stride_j_3;
    for (int64_t ctr_1 = 0; ctr_1 < _size_j_1; ctr_1 += 1) {
      float *RESTRICT _data_j_20_39_10 = _stride_j_1 * ctr_1 + _data_j_20_39;
      float *RESTRICT _data_j_20_33_10 = _stride_j_1 * ctr_1 + _data_j_20_33;
      float *RESTRICT _data_j_20_310_10 = _stride_j_1 * ctr_1 + _data_j_20_310;
      float *RESTRICT _data_j_20_35_10 = _stride_j_1 * ctr_1 + _data_j_20_35;
      float *RESTRICT _data_j_20_30_10 = _stride_j_1 * ctr_1 + _data_j_20_30;
      float *RESTRICT _data_j_20_36_10 = _stride_j_1 * ctr_1 + _data_j_20_36;
      float *RESTRICT _data_j_20_311_10 = _stride_j_1 * ctr_1 + _data_j_20_311;
      float *RESTRICT _data_j_20_34_10 = _stride_j_1 * ctr_1 + _data_j_20_34;
      float *RESTRICT _data_j_20_312_10 = _stride_j_1 * ctr_1 + _data_j_20_312;
      for (int64_t ctr_0 = 0; ctr_0 < _size_j_0; ctr_0 += 1) {
        _data_buffer[9 * _size_j_0 * _size_j_1 * ctr_2 + 9 * _size_j_0 * ctr_1 +
                     9 * ctr_0] = _data_j_20_39_10[_stride_j_0 * ctr_0];
        _data_buffer[9 * _size_j_0 * _size_j_1 * ctr_2 + 9 * _size_j_0 * ctr_1 +
                     9 * ctr_0 + 1] = _data_j_20_33_10[_stride_j_0 * ctr_0];
        _data_buffer[9 * _size_j_0 * _size_j_1 * ctr_2 + 9 * _size_j_0 * ctr_1 +
                     9 * ctr_0 + 2] = _data_j_20_310_10[_stride_j_0 * ctr_0];
        _data_buffer[9 * _size_j_0 * _size_j_1 * ctr_2 + 9 * _size_j_0 * ctr_1 +
                     9 * ctr_0 + 3] = _data_j_20_35_10[_stride_j_0 * ctr_0];
        _data_buffer[9 * _size_j_0 * _size_j_1 * ctr_2 + 9 * _size_j_0 * ctr_1 +
                     9 * ctr_0 + 4] = _data_j_20_30_10[_stride_j_0 * ctr_0];
        _data_buffer[9 * _size_j_0 * _size_j_1 * ctr_2 + 9 * _size_j_0 * ctr_1 +
                     9 * ctr_0 + 5] = _data_j_20_36_10[_stride_j_0 * ctr_0];
        _data_buffer[9 * _size_j_0 * _size_j_1 * ctr_2 + 9 * _size_j_0 * ctr_1 +
                     9 * ctr_0 + 6] = _data_j_20_311_10[_stride_j_0 * ctr_0];
        _data_buffer[9 * _size_j_0 * _size_j_1 * ctr_2 + 9 * _size_j_0 * ctr_1 +
                     9 * ctr_0 + 7] = _data_j_20_34_10[_stride_j_0 * ctr_0];
        _data_buffer[9 * _size_j_0 * _size_j_1 * ctr_2 + 9 * _size_j_0 * ctr_1 +
                     9 * ctr_0 + 8] = _data_j_20_312_10[_stride_j_0 * ctr_0];
      }
    }
  }
}
} // namespace internal_pack_W

namespace internal_pack_TW {
static FUNC_PREFIX void
pack_TW(float *RESTRICT _data_buffer, float *RESTRICT const _data_j,
        int64_t const _size_j_0, int64_t const _size_j_1,
        int64_t const _size_j_2, int64_t const _stride_j_0,
        int64_t const _stride_j_1, int64_t const _stride_j_2,
        int64_t const _stride_j_3) {
  for (int64_t ctr_2 = 0; ctr_2 < _size_j_2; ctr_2 += 1) {
    float *RESTRICT _data_j_20_310 =
        _data_j + _stride_j_2 * ctr_2 + 10 * _stride_j_3;
    float *RESTRICT _data_j_20_36 =
        _data_j + _stride_j_2 * ctr_2 + 6 * _stride_j_3;
    float *RESTRICT _data_j_20_312 =
        _data_j + _stride_j_2 * ctr_2 + 12 * _stride_j_3;
    for (int64_t ctr_1 = 0; ctr_1 < _size_j_1; ctr_1 += 1) {
      float *RESTRICT _data_j_20_310_10 = _stride_j_1 * ctr_1 + _data_j_20_310;
      float *RESTRICT _data_j_20_36_10 = _stride_j_1 * ctr_1 + _data_j_20_36;
      float *RESTRICT _data_j_20_312_10 = _stride_j_1 * ctr_1 + _data_j_20_312;
      for (int64_t ctr_0 = 0; ctr_0 < _size_j_0; ctr_0 += 1) {
        _data_buffer[3 * _size_j_0 * _size_j_1 * ctr_2 + 3 * _size_j_0 * ctr_1 +
                     3 * ctr_0] = _data_j_20_310_10[_stride_j_0 * ctr_0];
        _data_buffer[3 * _size_j_0 * _size_j_1 * ctr_2 + 3 * _size_j_0 * ctr_1 +
                     3 * ctr_0 + 1] = _data_j_20_36_10[_stride_j_0 * ctr_0];
        _data_buffer[3 * _size_j_0 * _size_j_1 * ctr_2 + 3 * _size_j_0 * ctr_1 +
                     3 * ctr_0 + 2] = _data_j_20_312_10[_stride_j_0 * ctr_0];
      }
    }
  }
}
} // namespace internal_pack_TW

namespace internal_pack_BNW {
static FUNC_PREFIX void
pack_BNW(float *RESTRICT _data_buffer, float *RESTRICT const _data_j,
         int64_t const _size_j_0, int64_t const _size_j_1,
         int64_t const _size_j_2, int64_t const _stride_j_0,
         int64_t const _stride_j_1, int64_t const _stride_j_2,
         int64_t const _stride_j_3) {
  for (int64_t ctr_2 = 0; ctr_2 < _size_j_2; ctr_2 += 1) {
    float *RESTRICT _data_j_20_311 =
        _data_j + _stride_j_2 * ctr_2 + 11 * _stride_j_3;
    for (int64_t ctr_1 = 0; ctr_1 < _size_j_1; ctr_1 += 1) {
      float *RESTRICT _data_j_20_311_10 = _stride_j_1 * ctr_1 + _data_j_20_311;
      for (int64_t ctr_0 = 0; ctr_0 < _size_j_0; ctr_0 += 1) {
        _data_buffer[_size_j_0 * _size_j_1 * ctr_2 + _size_j_0 * ctr_1 +
                     ctr_0] = _data_j_20_311_10[_stride_j_0 * ctr_0];
      }
    }
  }
}
} // namespace internal_pack_BNW

namespace internal_pack_NW {
static FUNC_PREFIX void
pack_NW(float *RESTRICT _data_buffer, float *RESTRICT const _data_j,
        int64_t const _size_j_0, int64_t const _size_j_1,
        int64_t const _size_j_2, int64_t const _stride_j_0,
        int64_t const _stride_j_1, int64_t const _stride_j_2,
        int64_t const _stride_j_3) {
  for (int64_t ctr_2 = 0; ctr_2 < _size_j_2; ctr_2 += 1) {
    float *RESTRICT _data_j_20_311 =
        _data_j + _stride_j_2 * ctr_2 + 11 * _stride_j_3;
    float *RESTRICT _data_j_20_34 =
        _data_j + _stride_j_2 * ctr_2 + 4 * _stride_j_3;
    float *RESTRICT _data_j_20_312 =
        _data_j + _stride_j_2 * ctr_2 + 12 * _stride_j_3;
    for (int64_t ctr_1 = 0; ctr_1 < _size_j_1; ctr_1 += 1) {
      float *RESTRICT _data_j_20_311_10 = _stride_j_1 * ctr_1 + _data_j_20_311;
      float *RESTRICT _data_j_20_34_10 = _stride_j_1 * ctr_1 + _data_j_20_34;
      float *RESTRICT _data_j_20_312_10 = _stride_j_1 * ctr_1 + _data_j_20_312;
      for (int64_t ctr_0 = 0; ctr_0 < _size_j_0; ctr_0 += 1) {
        _data_buffer[3 * _size_j_0 * _size_j_1 * ctr_2 + 3 * _size_j_0 * ctr_1 +
                     3 * ctr_0] = _data_j_20_311_10[_stride_j_0 * ctr_0];
        _data_buffer[3 * _size_j_0 * _size_j_1 * ctr_2 + 3 * _size_j_0 * ctr_1 +
                     3 * ctr_0 + 1] = _data_j_20_34_10[_stride_j_0 * ctr_0];
        _data_buffer[3 * _size_j_0 * _size_j_1 * ctr_2 + 3 * _size_j_0 * ctr_1 +
                     3 * ctr_0 + 2] = _data_j_20_312_10[_stride_j_0 * ctr_0];
      }
    }
  }
}
} // namespace internal_pack_NW

namespace internal_pack_TNW {
static FUNC_PREFIX void
pack_TNW(float *RESTRICT _data_buffer, float *RESTRICT const _data_j,
         int64_t const _size_j_0, int64_t const _size_j_1,
         int64_t const _size_j_2, int64_t const _stride_j_0,
         int64_t const _stride_j_1, int64_t const _stride_j_2,
         int64_t const _stride_j_3) {
  for (int64_t ctr_2 = 0; ctr_2 < _size_j_2; ctr_2 += 1) {
    float *RESTRICT _data_j_20_312 =
        _data_j + _stride_j_2 * ctr_2 + 12 * _stride_j_3;
    for (int64_t ctr_1 = 0; ctr_1 < _size_j_1; ctr_1 += 1) {
      float *RESTRICT _data_j_20_312_10 = _stride_j_1 * ctr_1 + _data_j_20_312;
      for (int64_t ctr_0 = 0; ctr_0 < _size_j_0; ctr_0 += 1) {
        _data_buffer[_size_j_0 * _size_j_1 * ctr_2 + _size_j_0 * ctr_1 +
                     ctr_0] = _data_j_20_312_10[_stride_j_0 * ctr_0];
      }
    }
  }
}
} // namespace internal_pack_TNW

namespace internal_pack_BS {
static FUNC_PREFIX void
pack_BS(float *RESTRICT _data_buffer, float *RESTRICT const _data_j,
        int64_t const _size_j_0, int64_t const _size_j_1,
        int64_t const _size_j_2, int64_t const _stride_j_0,
        int64_t const _stride_j_1, int64_t const _stride_j_2,
        int64_t const _stride_j_3) {
  for (int64_t ctr_2 = 0; ctr_2 < _size_j_2; ctr_2 += 1) {
    float *RESTRICT _data_j_20_39 =
        _data_j + _stride_j_2 * ctr_2 + 9 * _stride_j_3;
    float *RESTRICT _data_j_20_37 =
        _data_j + _stride_j_2 * ctr_2 + 7 * _stride_j_3;
    for (int64_t ctr_1 = 0; ctr_1 < _size_j_1; ctr_1 += 1) {
      float *RESTRICT _data_j_20_39_10 = _stride_j_1 * ctr_1 + _data_j_20_39;
      float *RESTRICT _data_j_20_37_10 = _stride_j_1 * ctr_1 + _data_j_20_37;
      for (int64_t ctr_0 = 0; ctr_0 < _size_j_0; ctr_0 += 1) {
        _data_buffer[2 * _size_j_0 * _size_j_1 * ctr_2 + 2 * _size_j_0 * ctr_1 +
                     2 * ctr_0] = _data_j_20_39_10[_stride_j_0 * ctr_0];
        _data_buffer[2 * _size_j_0 * _size_j_1 * ctr_2 + 2 * _size_j_0 * ctr_1 +
                     2 * ctr_0 + 1] = _data_j_20_37_10[_stride_j_0 * ctr_0];
      }
    }
  }
}
} // namespace internal_pack_BS

namespace internal_pack_S {
static FUNC_PREFIX void
pack_S(float *RESTRICT _data_buffer, float *RESTRICT const _data_j,
       int64_t const _size_j_0, int64_t const _size_j_1,
       int64_t const _size_j_2, int64_t const _stride_j_0,
       int64_t const _stride_j_1, int64_t const _stride_j_2,
       int64_t const _stride_j_3) {
  for (int64_t ctr_2 = 0; ctr_2 < _size_j_2; ctr_2 += 1) {
    float *RESTRICT _data_j_20_39 =
        _data_j + _stride_j_2 * ctr_2 + 9 * _stride_j_3;
    float *RESTRICT _data_j_20_33 =
        _data_j + _stride_j_2 * ctr_2 + 3 * _stride_j_3;
    float *RESTRICT _data_j_20_310 =
        _data_j + _stride_j_2 * ctr_2 + 10 * _stride_j_3;
    float *RESTRICT _data_j_20_37 =
        _data_j + _stride_j_2 * ctr_2 + 7 * _stride_j_3;
    float *RESTRICT _data_j_20_31 = _data_j + _stride_j_2 * ctr_2 + _stride_j_3;
    float *RESTRICT _data_j_20_38 =
        _data_j + _stride_j_2 * ctr_2 + 8 * _stride_j_3;
    for (int64_t ctr_1 = 0; ctr_1 < _size_j_1; ctr_1 += 1) {
      float *RESTRICT _data_j_20_39_10 = _stride_j_1 * ctr_1 + _data_j_20_39;
      float *RESTRICT _data_j_20_33_10 = _stride_j_1 * ctr_1 + _data_j_20_33;
      float *RESTRICT _data_j_20_310_10 = _stride_j_1 * ctr_1 + _data_j_20_310;
      float *RESTRICT _data_j_20_37_10 = _stride_j_1 * ctr_1 + _data_j_20_37;
      float *RESTRICT _data_j_20_31_10 = _stride_j_1 * ctr_1 + _data_j_20_31;
      float *RESTRICT _data_j_20_38_10 = _stride_j_1 * ctr_1 + _data_j_20_38;
      for (int64_t ctr_0 = 0; ctr_0 < _size_j_0; ctr_0 += 1) {
        _data_buffer[6 * _size_j_0 * _size_j_1 * ctr_2 + 6 * _size_j_0 * ctr_1 +
                     6 * ctr_0] = _data_j_20_39_10[_stride_j_0 * ctr_0];
        _data_buffer[6 * _size_j_0 * _size_j_1 * ctr_2 + 6 * _size_j_0 * ctr_1 +
                     6 * ctr_0 + 1] = _data_j_20_33_10[_stride_j_0 * ctr_0];
        _data_buffer[6 * _size_j_0 * _size_j_1 * ctr_2 + 6 * _size_j_0 * ctr_1 +
                     6 * ctr_0 + 2] = _data_j_20_310_10[_stride_j_0 * ctr_0];
        _data_buffer[6 * _size_j_0 * _size_j_1 * ctr_2 + 6 * _size_j_0 * ctr_1 +
                     6 * ctr_0 + 3] = _data_j_20_37_10[_stride_j_0 * ctr_0];
        _data_buffer[6 * _size_j_0 * _size_j_1 * ctr_2 + 6 * _size_j_0 * ctr_1 +
                     6 * ctr_0 + 4] = _data_j_20_31_10[_stride_j_0 * ctr_0];
        _data_buffer[6 * _size_j_0 * _size_j_1 * ctr_2 + 6 * _size_j_0 * ctr_1 +
                     6 * ctr_0 + 5] = _data_j_20_38_10[_stride_j_0 * ctr_0];
      }
    }
  }
}
} // namespace internal_pack_S

namespace internal_pack_TS {
static FUNC_PREFIX void
pack_TS(float *RESTRICT _data_buffer, float *RESTRICT const _data_j,
        int64_t const _size_j_0, int64_t const _size_j_1,
        int64_t const _size_j_2, int64_t const _stride_j_0,
        int64_t const _stride_j_1, int64_t const _stride_j_2,
        int64_t const _stride_j_3) {
  for (int64_t ctr_2 = 0; ctr_2 < _size_j_2; ctr_2 += 1) {
    float *RESTRICT _data_j_20_310 =
        _data_j + _stride_j_2 * ctr_2 + 10 * _stride_j_3;
    float *RESTRICT _data_j_20_38 =
        _data_j + _stride_j_2 * ctr_2 + 8 * _stride_j_3;
    for (int64_t ctr_1 = 0; ctr_1 < _size_j_1; ctr_1 += 1) {
      float *RESTRICT _data_j_20_310_10 = _stride_j_1 * ctr_1 + _data_j_20_310;
      float *RESTRICT _data_j_20_38_10 = _stride_j_1 * ctr_1 + _data_j_20_38;
      for (int64_t ctr_0 = 0; ctr_0 < _size_j_0; ctr_0 += 1) {
        _data_buffer[2 * _size_j_0 * _size_j_1 * ctr_2 + 2 * _size_j_0 * ctr_1 +
                     2 * ctr_0] = _data_j_20_310_10[_stride_j_0 * ctr_0];
        _data_buffer[2 * _size_j_0 * _size_j_1 * ctr_2 + 2 * _size_j_0 * ctr_1 +
                     2 * ctr_0 + 1] = _data_j_20_38_10[_stride_j_0 * ctr_0];
      }
    }
  }
}
} // namespace internal_pack_TS

namespace internal_pack_B {
static FUNC_PREFIX void
pack_B(float *RESTRICT _data_buffer, float *RESTRICT const _data_j,
       int64_t const _size_j_0, int64_t const _size_j_1,
       int64_t const _size_j_2, int64_t const _stride_j_0,
       int64_t const _stride_j_1, int64_t const _stride_j_2,
       int64_t const _stride_j_3) {
  for (int64_t ctr_2 = 0; ctr_2 < _size_j_2; ctr_2 += 1) {
    float *RESTRICT _data_j_20_39 =
        _data_j + _stride_j_2 * ctr_2 + 9 * _stride_j_3;
    float *RESTRICT _data_j_20_35 =
        _data_j + _stride_j_2 * ctr_2 + 5 * _stride_j_3;
    float *RESTRICT _data_j_20_311 =
        _data_j + _stride_j_2 * ctr_2 + 11 * _stride_j_3;
    float *RESTRICT _data_j_20_37 =
        _data_j + _stride_j_2 * ctr_2 + 7 * _stride_j_3;
    float *RESTRICT _data_j_20_32 =
        _data_j + _stride_j_2 * ctr_2 + 2 * _stride_j_3;
    for (int64_t ctr_1 = 0; ctr_1 < _size_j_1; ctr_1 += 1) {
      float *RESTRICT _data_j_20_39_10 = _stride_j_1 * ctr_1 + _data_j_20_39;
      float *RESTRICT _data_j_20_35_10 = _stride_j_1 * ctr_1 + _data_j_20_35;
      float *RESTRICT _data_j_20_311_10 = _stride_j_1 * ctr_1 + _data_j_20_311;
      float *RESTRICT _data_j_20_37_10 = _stride_j_1 * ctr_1 + _data_j_20_37;
      float *RESTRICT _data_j_20_32_10 = _stride_j_1 * ctr_1 + _data_j_20_32;
      for (int64_t ctr_0 = 0; ctr_0 < _size_j_0; ctr_0 += 1) {
        _data_buffer[5 * _size_j_0 * _size_j_1 * ctr_2 + 5 * _size_j_0 * ctr_1 +
                     5 * ctr_0] = _data_j_20_39_10[_stride_j_0 * ctr_0];
        _data_buffer[5 * _size_j_0 * _size_j_1 * ctr_2 + 5 * _size_j_0 * ctr_1 +
                     5 * ctr_0 + 1] = _data_j_20_35_10[_stride_j_0 * ctr_0];
        _data_buffer[5 * _size_j_0 * _size_j_1 * ctr_2 + 5 * _size_j_0 * ctr_1 +
                     5 * ctr_0 + 2] = _data_j_20_311_10[_stride_j_0 * ctr_0];
        _data_buffer[5 * _size_j_0 * _size_j_1 * ctr_2 + 5 * _size_j_0 * ctr_1 +
                     5 * ctr_0 + 3] = _data_j_20_37_10[_stride_j_0 * ctr_0];
        _data_buffer[5 * _size_j_0 * _size_j_1 * ctr_2 + 5 * _size_j_0 * ctr_1 +
                     5 * ctr_0 + 4] = _data_j_20_32_10[_stride_j_0 * ctr_0];
      }
    }
  }
}
} // namespace internal_pack_B

namespace internal_pack_T {
static FUNC_PREFIX void
pack_T(float *RESTRICT _data_buffer, float *RESTRICT const _data_j,
       int64_t const _size_j_0, int64_t const _size_j_1,
       int64_t const _size_j_2, int64_t const _stride_j_0,
       int64_t const _stride_j_1, int64_t const _stride_j_2,
       int64_t const _stride_j_3) {
  for (int64_t ctr_2 = 0; ctr_2 < _size_j_2; ctr_2 += 1) {
    float *RESTRICT _data_j_20_310 =
        _data_j + _stride_j_2 * ctr_2 + 10 * _stride_j_3;
    float *RESTRICT _data_j_20_36 =
        _data_j + _stride_j_2 * ctr_2 + 6 * _stride_j_3;
    float *RESTRICT _data_j_20_312 =
        _data_j + _stride_j_2 * ctr_2 + 12 * _stride_j_3;
    float *RESTRICT _data_j_20_38 =
        _data_j + _stride_j_2 * ctr_2 + 8 * _stride_j_3;
    for (int64_t ctr_1 = 0; ctr_1 < _size_j_1; ctr_1 += 1) {
      float *RESTRICT _data_j_20_310_10 = _stride_j_1 * ctr_1 + _data_j_20_310;
      float *RESTRICT _data_j_20_36_10 = _stride_j_1 * ctr_1 + _data_j_20_36;
      float *RESTRICT _data_j_20_312_10 = _stride_j_1 * ctr_1 + _data_j_20_312;
      float *RESTRICT _data_j_20_38_10 = _stride_j_1 * ctr_1 + _data_j_20_38;
      for (int64_t ctr_0 = 0; ctr_0 < _size_j_0; ctr_0 += 1) {
        _data_buffer[4 * _size_j_0 * _size_j_1 * ctr_2 + 4 * _size_j_0 * ctr_1 +
                     4 * ctr_0] = _data_j_20_310_10[_stride_j_0 * ctr_0];
        _data_buffer[4 * _size_j_0 * _size_j_1 * ctr_2 + 4 * _size_j_0 * ctr_1 +
                     4 * ctr_0 + 1] = _data_j_20_36_10[_stride_j_0 * ctr_0];
        _data_buffer[4 * _size_j_0 * _size_j_1 * ctr_2 + 4 * _size_j_0 * ctr_1 +
                     4 * ctr_0 + 2] = _data_j_20_312_10[_stride_j_0 * ctr_0];
        _data_buffer[4 * _size_j_0 * _size_j_1 * ctr_2 + 4 * _size_j_0 * ctr_1 +
                     4 * ctr_0 + 3] = _data_j_20_38_10[_stride_j_0 * ctr_0];
      }
    }
  }
}
} // namespace internal_pack_T

namespace internal_pack_BN {
static FUNC_PREFIX void
pack_BN(float *RESTRICT _data_buffer, float *RESTRICT const _data_j,
        int64_t const _size_j_0, int64_t const _size_j_1,
        int64_t const _size_j_2, int64_t const _stride_j_0,
        int64_t const _stride_j_1, int64_t const _stride_j_2,
        int64_t const _stride_j_3) {
  for (int64_t ctr_2 = 0; ctr_2 < _size_j_2; ctr_2 += 1) {
    float *RESTRICT _data_j_20_311 =
        _data_j + _stride_j_2 * ctr_2 + 11 * _stride_j_3;
    for (int64_t ctr_1 = 0; ctr_1 < _size_j_1; ctr_1 += 1) {
      float *RESTRICT _data_j_20_311_10 = _stride_j_1 * ctr_1 + _data_j_20_311;
      for (int64_t ctr_0 = 0; ctr_0 < _size_j_0; ctr_0 += 1) {
        _data_buffer[_size_j_0 * _size_j_1 * ctr_2 + _size_j_0 * ctr_1 +
                     ctr_0] = _data_j_20_311_10[_stride_j_0 * ctr_0];
      }
    }
  }
}
} // namespace internal_pack_BN

namespace internal_pack_N {
static FUNC_PREFIX void
pack_N(float *RESTRICT _data_buffer, float *RESTRICT const _data_j,
       int64_t const _size_j_0, int64_t const _size_j_1,
       int64_t const _size_j_2, int64_t const _stride_j_0,
       int64_t const _stride_j_1, int64_t const _stride_j_2,
       int64_t const _stride_j_3) {
  for (int64_t ctr_2 = 0; ctr_2 < _size_j_2; ctr_2 += 1) {
    float *RESTRICT _data_j_20_311 =
        _data_j + _stride_j_2 * ctr_2 + 11 * _stride_j_3;
    float *RESTRICT _data_j_20_34 =
        _data_j + _stride_j_2 * ctr_2 + 4 * _stride_j_3;
    float *RESTRICT _data_j_20_312 =
        _data_j + _stride_j_2 * ctr_2 + 12 * _stride_j_3;
    for (int64_t ctr_1 = 0; ctr_1 < _size_j_1; ctr_1 += 1) {
      float *RESTRICT _data_j_20_311_10 = _stride_j_1 * ctr_1 + _data_j_20_311;
      float *RESTRICT _data_j_20_34_10 = _stride_j_1 * ctr_1 + _data_j_20_34;
      float *RESTRICT _data_j_20_312_10 = _stride_j_1 * ctr_1 + _data_j_20_312;
      for (int64_t ctr_0 = 0; ctr_0 < _size_j_0; ctr_0 += 1) {
        _data_buffer[3 * _size_j_0 * _size_j_1 * ctr_2 + 3 * _size_j_0 * ctr_1 +
                     3 * ctr_0] = _data_j_20_311_10[_stride_j_0 * ctr_0];
        _data_buffer[3 * _size_j_0 * _size_j_1 * ctr_2 + 3 * _size_j_0 * ctr_1 +
                     3 * ctr_0 + 1] = _data_j_20_34_10[_stride_j_0 * ctr_0];
        _data_buffer[3 * _size_j_0 * _size_j_1 * ctr_2 + 3 * _size_j_0 * ctr_1 +
                     3 * ctr_0 + 2] = _data_j_20_312_10[_stride_j_0 * ctr_0];
      }
    }
  }
}
} // namespace internal_pack_N

namespace internal_pack_TN {
static FUNC_PREFIX void
pack_TN(float *RESTRICT _data_buffer, float *RESTRICT const _data_j,
        int64_t const _size_j_0, int64_t const _size_j_1,
        int64_t const _size_j_2, int64_t const _stride_j_0,
        int64_t const _stride_j_1, int64_t const _stride_j_2,
        int64_t const _stride_j_3) {
  for (int64_t ctr_2 = 0; ctr_2 < _size_j_2; ctr_2 += 1) {
    float *RESTRICT _data_j_20_312 =
        _data_j + _stride_j_2 * ctr_2 + 12 * _stride_j_3;
    for (int64_t ctr_1 = 0; ctr_1 < _size_j_1; ctr_1 += 1) {
      float *RESTRICT _data_j_20_312_10 = _stride_j_1 * ctr_1 + _data_j_20_312;
      for (int64_t ctr_0 = 0; ctr_0 < _size_j_0; ctr_0 += 1) {
        _data_buffer[_size_j_0 * _size_j_1 * ctr_2 + _size_j_0 * ctr_1 +
                     ctr_0] = _data_j_20_312_10[_stride_j_0 * ctr_0];
      }
    }
  }
}
} // namespace internal_pack_TN

namespace internal_unpack_BSW {
static FUNC_PREFIX void
unpack_BSW(float *RESTRICT const _data_buffer, float *RESTRICT _data_j,
           int64_t const _size_j_0, int64_t const _size_j_1,
           int64_t const _size_j_2, int64_t const _stride_j_0,
           int64_t const _stride_j_1, int64_t const _stride_j_2,
           int64_t const _stride_j_3) {
  for (int64_t ctr_2 = 0; ctr_2 < _size_j_2; ctr_2 += 1) {
    float *RESTRICT _data_j_20_39 =
        _data_j + _stride_j_2 * ctr_2 + 9 * _stride_j_3;
    for (int64_t ctr_1 = 0; ctr_1 < _size_j_1; ctr_1 += 1) {
      float *RESTRICT _data_j_20_39_10 = _stride_j_1 * ctr_1 + _data_j_20_39;
      for (int64_t ctr_0 = 0; ctr_0 < _size_j_0; ctr_0 += 1) {
        _data_j_20_39_10[_stride_j_0 * ctr_0] =
            _data_buffer[_size_j_0 * _size_j_1 * ctr_2 + _size_j_0 * ctr_1 +
                         ctr_0];
      }
    }
  }
}
} // namespace internal_unpack_BSW

namespace internal_unpack_SW {
static FUNC_PREFIX void
unpack_SW(float *RESTRICT const _data_buffer, float *RESTRICT _data_j,
          int64_t const _size_j_0, int64_t const _size_j_1,
          int64_t const _size_j_2, int64_t const _stride_j_0,
          int64_t const _stride_j_1, int64_t const _stride_j_2,
          int64_t const _stride_j_3) {
  for (int64_t ctr_2 = 0; ctr_2 < _size_j_2; ctr_2 += 1) {
    float *RESTRICT _data_j_20_39 =
        _data_j + _stride_j_2 * ctr_2 + 9 * _stride_j_3;
    float *RESTRICT _data_j_20_33 =
        _data_j + _stride_j_2 * ctr_2 + 3 * _stride_j_3;
    float *RESTRICT _data_j_20_310 =
        _data_j + _stride_j_2 * ctr_2 + 10 * _stride_j_3;
    for (int64_t ctr_1 = 0; ctr_1 < _size_j_1; ctr_1 += 1) {
      float *RESTRICT _data_j_20_39_10 = _stride_j_1 * ctr_1 + _data_j_20_39;
      float *RESTRICT _data_j_20_33_10 = _stride_j_1 * ctr_1 + _data_j_20_33;
      float *RESTRICT _data_j_20_310_10 = _stride_j_1 * ctr_1 + _data_j_20_310;
      for (int64_t ctr_0 = 0; ctr_0 < _size_j_0; ctr_0 += 1) {
        _data_j_20_39_10[_stride_j_0 * ctr_0] =
            _data_buffer[3 * _size_j_0 * _size_j_1 * ctr_2 +
                         3 * _size_j_0 * ctr_1 + 3 * ctr_0];
        _data_j_20_33_10[_stride_j_0 * ctr_0] =
            _data_buffer[3 * _size_j_0 * _size_j_1 * ctr_2 +
                         3 * _size_j_0 * ctr_1 + 3 * ctr_0 + 1];
        _data_j_20_310_10[_stride_j_0 * ctr_0] =
            _data_buffer[3 * _size_j_0 * _size_j_1 * ctr_2 +
                         3 * _size_j_0 * ctr_1 + 3 * ctr_0 + 2];
      }
    }
  }
}
} // namespace internal_unpack_SW

namespace internal_unpack_TSW {
static FUNC_PREFIX void
unpack_TSW(float *RESTRICT const _data_buffer, float *RESTRICT _data_j,
           int64_t const _size_j_0, int64_t const _size_j_1,
           int64_t const _size_j_2, int64_t const _stride_j_0,
           int64_t const _stride_j_1, int64_t const _stride_j_2,
           int64_t const _stride_j_3) {
  for (int64_t ctr_2 = 0; ctr_2 < _size_j_2; ctr_2 += 1) {
    float *RESTRICT _data_j_20_310 =
        _data_j + _stride_j_2 * ctr_2 + 10 * _stride_j_3;
    for (int64_t ctr_1 = 0; ctr_1 < _size_j_1; ctr_1 += 1) {
      float *RESTRICT _data_j_20_310_10 = _stride_j_1 * ctr_1 + _data_j_20_310;
      for (int64_t ctr_0 = 0; ctr_0 < _size_j_0; ctr_0 += 1) {
        _data_j_20_310_10[_stride_j_0 * ctr_0] =
            _data_buffer[_size_j_0 * _size_j_1 * ctr_2 + _size_j_0 * ctr_1 +
                         ctr_0];
      }
    }
  }
}
} // namespace internal_unpack_TSW

namespace internal_unpack_BW {
static FUNC_PREFIX void
unpack_BW(float *RESTRICT const _data_buffer, float *RESTRICT _data_j,
          int64_t const _size_j_0, int64_t const _size_j_1,
          int64_t const _size_j_2, int64_t const _stride_j_0,
          int64_t const _stride_j_1, int64_t const _stride_j_2,
          int64_t const _stride_j_3) {
  for (int64_t ctr_2 = 0; ctr_2 < _size_j_2; ctr_2 += 1) {
    float *RESTRICT _data_j_20_39 =
        _data_j + _stride_j_2 * ctr_2 + 9 * _stride_j_3;
    float *RESTRICT _data_j_20_35 =
        _data_j + _stride_j_2 * ctr_2 + 5 * _stride_j_3;
    float *RESTRICT _data_j_20_311 =
        _data_j + _stride_j_2 * ctr_2 + 11 * _stride_j_3;
    for (int64_t ctr_1 = 0; ctr_1 < _size_j_1; ctr_1 += 1) {
      float *RESTRICT _data_j_20_39_10 = _stride_j_1 * ctr_1 + _data_j_20_39;
      float *RESTRICT _data_j_20_35_10 = _stride_j_1 * ctr_1 + _data_j_20_35;
      float *RESTRICT _data_j_20_311_10 = _stride_j_1 * ctr_1 + _data_j_20_311;
      for (int64_t ctr_0 = 0; ctr_0 < _size_j_0; ctr_0 += 1) {
        _data_j_20_39_10[_stride_j_0 * ctr_0] =
            _data_buffer[3 * _size_j_0 * _size_j_1 * ctr_2 +
                         3 * _size_j_0 * ctr_1 + 3 * ctr_0];
        _data_j_20_35_10[_stride_j_0 * ctr_0] =
            _data_buffer[3 * _size_j_0 * _size_j_1 * ctr_2 +
                         3 * _size_j_0 * ctr_1 + 3 * ctr_0 + 1];
        _data_j_20_311_10[_stride_j_0 * ctr_0] =
            _data_buffer[3 * _size_j_0 * _size_j_1 * ctr_2 +
                         3 * _size_j_0 * ctr_1 + 3 * ctr_0 + 2];
      }
    }
  }
}
} // namespace internal_unpack_BW

namespace internal_unpack_W {
static FUNC_PREFIX void
unpack_W(float *RESTRICT const _data_buffer, float *RESTRICT _data_j,
         int64_t const _size_j_0, int64_t const _size_j_1,
         int64_t const _size_j_2, int64_t const _stride_j_0,
         int64_t const _stride_j_1, int64_t const _stride_j_2,
         int64_t const _stride_j_3) {
  for (int64_t ctr_2 = 0; ctr_2 < _size_j_2; ctr_2 += 1) {
    float *RESTRICT _data_j_20_39 =
        _data_j + _stride_j_2 * ctr_2 + 9 * _stride_j_3;
    float *RESTRICT _data_j_20_33 =
        _data_j + _stride_j_2 * ctr_2 + 3 * _stride_j_3;
    float *RESTRICT _data_j_20_310 =
        _data_j + _stride_j_2 * ctr_2 + 10 * _stride_j_3;
    float *RESTRICT _data_j_20_35 =
        _data_j + _stride_j_2 * ctr_2 + 5 * _stride_j_3;
    float *RESTRICT _data_j_20_30 = _data_j + _stride_j_2 * ctr_2;
    float *RESTRICT _data_j_20_36 =
        _data_j + _stride_j_2 * ctr_2 + 6 * _stride_j_3;
    float *RESTRICT _data_j_20_311 =
        _data_j + _stride_j_2 * ctr_2 + 11 * _stride_j_3;
    float *RESTRICT _data_j_20_34 =
        _data_j + _stride_j_2 * ctr_2 + 4 * _stride_j_3;
    float *RESTRICT _data_j_20_312 =
        _data_j + _stride_j_2 * ctr_2 + 12 * _stride_j_3;
    for (int64_t ctr_1 = 0; ctr_1 < _size_j_1; ctr_1 += 1) {
      float *RESTRICT _data_j_20_39_10 = _stride_j_1 * ctr_1 + _data_j_20_39;
      float *RESTRICT _data_j_20_33_10 = _stride_j_1 * ctr_1 + _data_j_20_33;
      float *RESTRICT _data_j_20_310_10 = _stride_j_1 * ctr_1 + _data_j_20_310;
      float *RESTRICT _data_j_20_35_10 = _stride_j_1 * ctr_1 + _data_j_20_35;
      float *RESTRICT _data_j_20_30_10 = _stride_j_1 * ctr_1 + _data_j_20_30;
      float *RESTRICT _data_j_20_36_10 = _stride_j_1 * ctr_1 + _data_j_20_36;
      float *RESTRICT _data_j_20_311_10 = _stride_j_1 * ctr_1 + _data_j_20_311;
      float *RESTRICT _data_j_20_34_10 = _stride_j_1 * ctr_1 + _data_j_20_34;
      float *RESTRICT _data_j_20_312_10 = _stride_j_1 * ctr_1 + _data_j_20_312;
      for (int64_t ctr_0 = 0; ctr_0 < _size_j_0; ctr_0 += 1) {
        _data_j_20_39_10[_stride_j_0 * ctr_0] =
            _data_buffer[9 * _size_j_0 * _size_j_1 * ctr_2 +
                         9 * _size_j_0 * ctr_1 + 9 * ctr_0];
        _data_j_20_33_10[_stride_j_0 * ctr_0] =
            _data_buffer[9 * _size_j_0 * _size_j_1 * ctr_2 +
                         9 * _size_j_0 * ctr_1 + 9 * ctr_0 + 1];
        _data_j_20_310_10[_stride_j_0 * ctr_0] =
            _data_buffer[9 * _size_j_0 * _size_j_1 * ctr_2 +
                         9 * _size_j_0 * ctr_1 + 9 * ctr_0 + 2];
        _data_j_20_35_10[_stride_j_0 * ctr_0] =
            _data_buffer[9 * _size_j_0 * _size_j_1 * ctr_2 +
                         9 * _size_j_0 * ctr_1 + 9 * ctr_0 + 3];
        _data_j_20_30_10[_stride_j_0 * ctr_0] =
            _data_buffer[9 * _size_j_0 * _size_j_1 * ctr_2 +
                         9 * _size_j_0 * ctr_1 + 9 * ctr_0 + 4];
        _data_j_20_36_10[_stride_j_0 * ctr_0] =
            _data_buffer[9 * _size_j_0 * _size_j_1 * ctr_2 +
                         9 * _size_j_0 * ctr_1 + 9 * ctr_0 + 5];
        _data_j_20_311_10[_stride_j_0 * ctr_0] =
            _data_buffer[9 * _size_j_0 * _size_j_1 * ctr_2 +
                         9 * _size_j_0 * ctr_1 + 9 * ctr_0 + 6];
        _data_j_20_34_10[_stride_j_0 * ctr_0] =
            _data_buffer[9 * _size_j_0 * _size_j_1 * ctr_2 +
                         9 * _size_j_0 * ctr_1 + 9 * ctr_0 + 7];
        _data_j_20_312_10[_stride_j_0 * ctr_0] =
            _data_buffer[9 * _size_j_0 * _size_j_1 * ctr_2 +
                         9 * _size_j_0 * ctr_1 + 9 * ctr_0 + 8];
      }
    }
  }
}
} // namespace internal_unpack_W

namespace internal_unpack_TW {
static FUNC_PREFIX void
unpack_TW(float *RESTRICT const _data_buffer, float *RESTRICT _data_j,
          int64_t const _size_j_0, int64_t const _size_j_1,
          int64_t const _size_j_2, int64_t const _stride_j_0,
          int64_t const _stride_j_1, int64_t const _stride_j_2,
          int64_t const _stride_j_3) {
  for (int64_t ctr_2 = 0; ctr_2 < _size_j_2; ctr_2 += 1) {
    float *RESTRICT _data_j_20_310 =
        _data_j + _stride_j_2 * ctr_2 + 10 * _stride_j_3;
    float *RESTRICT _data_j_20_36 =
        _data_j + _stride_j_2 * ctr_2 + 6 * _stride_j_3;
    float *RESTRICT _data_j_20_312 =
        _data_j + _stride_j_2 * ctr_2 + 12 * _stride_j_3;
    for (int64_t ctr_1 = 0; ctr_1 < _size_j_1; ctr_1 += 1) {
      float *RESTRICT _data_j_20_310_10 = _stride_j_1 * ctr_1 + _data_j_20_310;
      float *RESTRICT _data_j_20_36_10 = _stride_j_1 * ctr_1 + _data_j_20_36;
      float *RESTRICT _data_j_20_312_10 = _stride_j_1 * ctr_1 + _data_j_20_312;
      for (int64_t ctr_0 = 0; ctr_0 < _size_j_0; ctr_0 += 1) {
        _data_j_20_310_10[_stride_j_0 * ctr_0] =
            _data_buffer[3 * _size_j_0 * _size_j_1 * ctr_2 +
                         3 * _size_j_0 * ctr_1 + 3 * ctr_0];
        _data_j_20_36_10[_stride_j_0 * ctr_0] =
            _data_buffer[3 * _size_j_0 * _size_j_1 * ctr_2 +
                         3 * _size_j_0 * ctr_1 + 3 * ctr_0 + 1];
        _data_j_20_312_10[_stride_j_0 * ctr_0] =
            _data_buffer[3 * _size_j_0 * _size_j_1 * ctr_2 +
                         3 * _size_j_0 * ctr_1 + 3 * ctr_0 + 2];
      }
    }
  }
}
} // namespace internal_unpack_TW

namespace internal_unpack_BNW {
static FUNC_PREFIX void
unpack_BNW(float *RESTRICT const _data_buffer, float *RESTRICT _data_j,
           int64_t const _size_j_0, int64_t const _size_j_1,
           int64_t const _size_j_2, int64_t const _stride_j_0,
           int64_t const _stride_j_1, int64_t const _stride_j_2,
           int64_t const _stride_j_3) {
  for (int64_t ctr_2 = 0; ctr_2 < _size_j_2; ctr_2 += 1) {
    float *RESTRICT _data_j_20_311 =
        _data_j + _stride_j_2 * ctr_2 + 11 * _stride_j_3;
    for (int64_t ctr_1 = 0; ctr_1 < _size_j_1; ctr_1 += 1) {
      float *RESTRICT _data_j_20_311_10 = _stride_j_1 * ctr_1 + _data_j_20_311;
      for (int64_t ctr_0 = 0; ctr_0 < _size_j_0; ctr_0 += 1) {
        _data_j_20_311_10[_stride_j_0 * ctr_0] =
            _data_buffer[_size_j_0 * _size_j_1 * ctr_2 + _size_j_0 * ctr_1 +
                         ctr_0];
      }
    }
  }
}
} // namespace internal_unpack_BNW

namespace internal_unpack_NW {
static FUNC_PREFIX void
unpack_NW(float *RESTRICT const _data_buffer, float *RESTRICT _data_j,
          int64_t const _size_j_0, int64_t const _size_j_1,
          int64_t const _size_j_2, int64_t const _stride_j_0,
          int64_t const _stride_j_1, int64_t const _stride_j_2,
          int64_t const _stride_j_3) {
  for (int64_t ctr_2 = 0; ctr_2 < _size_j_2; ctr_2 += 1) {
    float *RESTRICT _data_j_20_311 =
        _data_j + _stride_j_2 * ctr_2 + 11 * _stride_j_3;
    float *RESTRICT _data_j_20_34 =
        _data_j + _stride_j_2 * ctr_2 + 4 * _stride_j_3;
    float *RESTRICT _data_j_20_312 =
        _data_j + _stride_j_2 * ctr_2 + 12 * _stride_j_3;
    for (int64_t ctr_1 = 0; ctr_1 < _size_j_1; ctr_1 += 1) {
      float *RESTRICT _data_j_20_311_10 = _stride_j_1 * ctr_1 + _data_j_20_311;
      float *RESTRICT _data_j_20_34_10 = _stride_j_1 * ctr_1 + _data_j_20_34;
      float *RESTRICT _data_j_20_312_10 = _stride_j_1 * ctr_1 + _data_j_20_312;
      for (int64_t ctr_0 = 0; ctr_0 < _size_j_0; ctr_0 += 1) {
        _data_j_20_311_10[_stride_j_0 * ctr_0] =
            _data_buffer[3 * _size_j_0 * _size_j_1 * ctr_2 +
                         3 * _size_j_0 * ctr_1 + 3 * ctr_0];
        _data_j_20_34_10[_stride_j_0 * ctr_0] =
            _data_buffer[3 * _size_j_0 * _size_j_1 * ctr_2 +
                         3 * _size_j_0 * ctr_1 + 3 * ctr_0 + 1];
        _data_j_20_312_10[_stride_j_0 * ctr_0] =
            _data_buffer[3 * _size_j_0 * _size_j_1 * ctr_2 +
                         3 * _size_j_0 * ctr_1 + 3 * ctr_0 + 2];
      }
    }
  }
}
} // namespace internal_unpack_NW

namespace internal_unpack_TNW {
static FUNC_PREFIX void
unpack_TNW(float *RESTRICT const _data_buffer, float *RESTRICT _data_j,
           int64_t const _size_j_0, int64_t const _size_j_1,
           int64_t const _size_j_2, int64_t const _stride_j_0,
           int64_t const _stride_j_1, int64_t const _stride_j_2,
           int64_t const _stride_j_3) {
  for (int64_t ctr_2 = 0; ctr_2 < _size_j_2; ctr_2 += 1) {
    float *RESTRICT _data_j_20_312 =
        _data_j + _stride_j_2 * ctr_2 + 12 * _stride_j_3;
    for (int64_t ctr_1 = 0; ctr_1 < _size_j_1; ctr_1 += 1) {
      float *RESTRICT _data_j_20_312_10 = _stride_j_1 * ctr_1 + _data_j_20_312;
      for (int64_t ctr_0 = 0; ctr_0 < _size_j_0; ctr_0 += 1) {
        _data_j_20_312_10[_stride_j_0 * ctr_0] =
            _data_buffer[_size_j_0 * _size_j_1 * ctr_2 + _size_j_0 * ctr_1 +
                         ctr_0];
      }
    }
  }
}
} // namespace internal_unpack_TNW

namespace internal_unpack_BS {
static FUNC_PREFIX void
unpack_BS(float *RESTRICT const _data_buffer, float *RESTRICT _data_j,
          int64_t const _size_j_0, int64_t const _size_j_1,
          int64_t const _size_j_2, int64_t const _stride_j_0,
          int64_t const _stride_j_1, int64_t const _stride_j_2,
          int64_t const _stride_j_3) {
  for (int64_t ctr_2 = 0; ctr_2 < _size_j_2; ctr_2 += 1) {
    float *RESTRICT _data_j_20_39 =
        _data_j + _stride_j_2 * ctr_2 + 9 * _stride_j_3;
    float *RESTRICT _data_j_20_37 =
        _data_j + _stride_j_2 * ctr_2 + 7 * _stride_j_3;
    for (int64_t ctr_1 = 0; ctr_1 < _size_j_1; ctr_1 += 1) {
      float *RESTRICT _data_j_20_39_10 = _stride_j_1 * ctr_1 + _data_j_20_39;
      float *RESTRICT _data_j_20_37_10 = _stride_j_1 * ctr_1 + _data_j_20_37;
      for (int64_t ctr_0 = 0; ctr_0 < _size_j_0; ctr_0 += 1) {
        _data_j_20_39_10[_stride_j_0 * ctr_0] =
            _data_buffer[2 * _size_j_0 * _size_j_1 * ctr_2 +
                         2 * _size_j_0 * ctr_1 + 2 * ctr_0];
        _data_j_20_37_10[_stride_j_0 * ctr_0] =
            _data_buffer[2 * _size_j_0 * _size_j_1 * ctr_2 +
                         2 * _size_j_0 * ctr_1 + 2 * ctr_0 + 1];
      }
    }
  }
}
} // namespace internal_unpack_BS

namespace internal_unpack_S {
static FUNC_PREFIX void
unpack_S(float *RESTRICT const _data_buffer, float *RESTRICT _data_j,
         int64_t const _size_j_0, int64_t const _size_j_1,
         int64_t const _size_j_2, int64_t const _stride_j_0,
         int64_t const _stride_j_1, int64_t const _stride_j_2,
         int64_t const _stride_j_3) {
  for (int64_t ctr_2 = 0; ctr_2 < _size_j_2; ctr_2 += 1) {
    float *RESTRICT _data_j_20_39 =
        _data_j + _stride_j_2 * ctr_2 + 9 * _stride_j_3;
    float *RESTRICT _data_j_20_33 =
        _data_j + _stride_j_2 * ctr_2 + 3 * _stride_j_3;
    float *RESTRICT _data_j_20_310 =
        _data_j + _stride_j_2 * ctr_2 + 10 * _stride_j_3;
    float *RESTRICT _data_j_20_37 =
        _data_j + _stride_j_2 * ctr_2 + 7 * _stride_j_3;
    float *RESTRICT _data_j_20_31 = _data_j + _stride_j_2 * ctr_2 + _stride_j_3;
    float *RESTRICT _data_j_20_38 =
        _data_j + _stride_j_2 * ctr_2 + 8 * _stride_j_3;
    for (int64_t ctr_1 = 0; ctr_1 < _size_j_1; ctr_1 += 1) {
      float *RESTRICT _data_j_20_39_10 = _stride_j_1 * ctr_1 + _data_j_20_39;
      float *RESTRICT _data_j_20_33_10 = _stride_j_1 * ctr_1 + _data_j_20_33;
      float *RESTRICT _data_j_20_310_10 = _stride_j_1 * ctr_1 + _data_j_20_310;
      float *RESTRICT _data_j_20_37_10 = _stride_j_1 * ctr_1 + _data_j_20_37;
      float *RESTRICT _data_j_20_31_10 = _stride_j_1 * ctr_1 + _data_j_20_31;
      float *RESTRICT _data_j_20_38_10 = _stride_j_1 * ctr_1 + _data_j_20_38;
      for (int64_t ctr_0 = 0; ctr_0 < _size_j_0; ctr_0 += 1) {
        _data_j_20_39_10[_stride_j_0 * ctr_0] =
            _data_buffer[6 * _size_j_0 * _size_j_1 * ctr_2 +
                         6 * _size_j_0 * ctr_1 + 6 * ctr_0];
        _data_j_20_33_10[_stride_j_0 * ctr_0] =
            _data_buffer[6 * _size_j_0 * _size_j_1 * ctr_2 +
                         6 * _size_j_0 * ctr_1 + 6 * ctr_0 + 1];
        _data_j_20_310_10[_stride_j_0 * ctr_0] =
            _data_buffer[6 * _size_j_0 * _size_j_1 * ctr_2 +
                         6 * _size_j_0 * ctr_1 + 6 * ctr_0 + 2];
        _data_j_20_37_10[_stride_j_0 * ctr_0] =
            _data_buffer[6 * _size_j_0 * _size_j_1 * ctr_2 +
                         6 * _size_j_0 * ctr_1 + 6 * ctr_0 + 3];
        _data_j_20_31_10[_stride_j_0 * ctr_0] =
            _data_buffer[6 * _size_j_0 * _size_j_1 * ctr_2 +
                         6 * _size_j_0 * ctr_1 + 6 * ctr_0 + 4];
        _data_j_20_38_10[_stride_j_0 * ctr_0] =
            _data_buffer[6 * _size_j_0 * _size_j_1 * ctr_2 +
                         6 * _size_j_0 * ctr_1 + 6 * ctr_0 + 5];
      }
    }
  }
}
} // namespace internal_unpack_S

namespace internal_unpack_TS {
static FUNC_PREFIX void
unpack_TS(float *RESTRICT const _data_buffer, float *RESTRICT _data_j,
          int64_t const _size_j_0, int64_t const _size_j_1,
          int64_t const _size_j_2, int64_t const _stride_j_0,
          int64_t const _stride_j_1, int64_t const _stride_j_2,
          int64_t const _stride_j_3) {
  for (int64_t ctr_2 = 0; ctr_2 < _size_j_2; ctr_2 += 1) {
    float *RESTRICT _data_j_20_310 =
        _data_j + _stride_j_2 * ctr_2 + 10 * _stride_j_3;
    float *RESTRICT _data_j_20_38 =
        _data_j + _stride_j_2 * ctr_2 + 8 * _stride_j_3;
    for (int64_t ctr_1 = 0; ctr_1 < _size_j_1; ctr_1 += 1) {
      float *RESTRICT _data_j_20_310_10 = _stride_j_1 * ctr_1 + _data_j_20_310;
      float *RESTRICT _data_j_20_38_10 = _stride_j_1 * ctr_1 + _data_j_20_38;
      for (int64_t ctr_0 = 0; ctr_0 < _size_j_0; ctr_0 += 1) {
        _data_j_20_310_10[_stride_j_0 * ctr_0] =
            _data_buffer[2 * _size_j_0 * _size_j_1 * ctr_2 +
                         2 * _size_j_0 * ctr_1 + 2 * ctr_0];
        _data_j_20_38_10[_stride_j_0 * ctr_0] =
            _data_buffer[2 * _size_j_0 * _size_j_1 * ctr_2 +
                         2 * _size_j_0 * ctr_1 + 2 * ctr_0 + 1];
      }
    }
  }
}
} // namespace internal_unpack_TS

namespace internal_unpack_B {
static FUNC_PREFIX void
unpack_B(float *RESTRICT const _data_buffer, float *RESTRICT _data_j,
         int64_t const _size_j_0, int64_t const _size_j_1,
         int64_t const _size_j_2, int64_t const _stride_j_0,
         int64_t const _stride_j_1, int64_t const _stride_j_2,
         int64_t const _stride_j_3) {
  for (int64_t ctr_2 = 0; ctr_2 < _size_j_2; ctr_2 += 1) {
    float *RESTRICT _data_j_20_39 =
        _data_j + _stride_j_2 * ctr_2 + 9 * _stride_j_3;
    float *RESTRICT _data_j_20_35 =
        _data_j + _stride_j_2 * ctr_2 + 5 * _stride_j_3;
    float *RESTRICT _data_j_20_311 =
        _data_j + _stride_j_2 * ctr_2 + 11 * _stride_j_3;
    float *RESTRICT _data_j_20_37 =
        _data_j + _stride_j_2 * ctr_2 + 7 * _stride_j_3;
    float *RESTRICT _data_j_20_32 =
        _data_j + _stride_j_2 * ctr_2 + 2 * _stride_j_3;
    for (int64_t ctr_1 = 0; ctr_1 < _size_j_1; ctr_1 += 1) {
      float *RESTRICT _data_j_20_39_10 = _stride_j_1 * ctr_1 + _data_j_20_39;
      float *RESTRICT _data_j_20_35_10 = _stride_j_1 * ctr_1 + _data_j_20_35;
      float *RESTRICT _data_j_20_311_10 = _stride_j_1 * ctr_1 + _data_j_20_311;
      float *RESTRICT _data_j_20_37_10 = _stride_j_1 * ctr_1 + _data_j_20_37;
      float *RESTRICT _data_j_20_32_10 = _stride_j_1 * ctr_1 + _data_j_20_32;
      for (int64_t ctr_0 = 0; ctr_0 < _size_j_0; ctr_0 += 1) {
        _data_j_20_39_10[_stride_j_0 * ctr_0] =
            _data_buffer[5 * _size_j_0 * _size_j_1 * ctr_2 +
                         5 * _size_j_0 * ctr_1 + 5 * ctr_0];
        _data_j_20_35_10[_stride_j_0 * ctr_0] =
            _data_buffer[5 * _size_j_0 * _size_j_1 * ctr_2 +
                         5 * _size_j_0 * ctr_1 + 5 * ctr_0 + 1];
        _data_j_20_311_10[_stride_j_0 * ctr_0] =
            _data_buffer[5 * _size_j_0 * _size_j_1 * ctr_2 +
                         5 * _size_j_0 * ctr_1 + 5 * ctr_0 + 2];
        _data_j_20_37_10[_stride_j_0 * ctr_0] =
            _data_buffer[5 * _size_j_0 * _size_j_1 * ctr_2 +
                         5 * _size_j_0 * ctr_1 + 5 * ctr_0 + 3];
        _data_j_20_32_10[_stride_j_0 * ctr_0] =
            _data_buffer[5 * _size_j_0 * _size_j_1 * ctr_2 +
                         5 * _size_j_0 * ctr_1 + 5 * ctr_0 + 4];
      }
    }
  }
}
} // namespace internal_unpack_B

namespace internal_unpack_T {
static FUNC_PREFIX void
unpack_T(float *RESTRICT const _data_buffer, float *RESTRICT _data_j,
         int64_t const _size_j_0, int64_t const _size_j_1,
         int64_t const _size_j_2, int64_t const _stride_j_0,
         int64_t const _stride_j_1, int64_t const _stride_j_2,
         int64_t const _stride_j_3) {
  for (int64_t ctr_2 = 0; ctr_2 < _size_j_2; ctr_2 += 1) {
    float *RESTRICT _data_j_20_310 =
        _data_j + _stride_j_2 * ctr_2 + 10 * _stride_j_3;
    float *RESTRICT _data_j_20_36 =
        _data_j + _stride_j_2 * ctr_2 + 6 * _stride_j_3;
    float *RESTRICT _data_j_20_312 =
        _data_j + _stride_j_2 * ctr_2 + 12 * _stride_j_3;
    float *RESTRICT _data_j_20_38 =
        _data_j + _stride_j_2 * ctr_2 + 8 * _stride_j_3;
    for (int64_t ctr_1 = 0; ctr_1 < _size_j_1; ctr_1 += 1) {
      float *RESTRICT _data_j_20_310_10 = _stride_j_1 * ctr_1 + _data_j_20_310;
      float *RESTRICT _data_j_20_36_10 = _stride_j_1 * ctr_1 + _data_j_20_36;
      float *RESTRICT _data_j_20_312_10 = _stride_j_1 * ctr_1 + _data_j_20_312;
      float *RESTRICT _data_j_20_38_10 = _stride_j_1 * ctr_1 + _data_j_20_38;
      for (int64_t ctr_0 = 0; ctr_0 < _size_j_0; ctr_0 += 1) {
        _data_j_20_310_10[_stride_j_0 * ctr_0] =
            _data_buffer[4 * _size_j_0 * _size_j_1 * ctr_2 +
                         4 * _size_j_0 * ctr_1 + 4 * ctr_0];
        _data_j_20_36_10[_stride_j_0 * ctr_0] =
            _data_buffer[4 * _size_j_0 * _size_j_1 * ctr_2 +
                         4 * _size_j_0 * ctr_1 + 4 * ctr_0 + 1];
        _data_j_20_312_10[_stride_j_0 * ctr_0] =
            _data_buffer[4 * _size_j_0 * _size_j_1 * ctr_2 +
                         4 * _size_j_0 * ctr_1 + 4 * ctr_0 + 2];
        _data_j_20_38_10[_stride_j_0 * ctr_0] =
            _data_buffer[4 * _size_j_0 * _size_j_1 * ctr_2 +
                         4 * _size_j_0 * ctr_1 + 4 * ctr_0 + 3];
      }
    }
  }
}
} // namespace internal_unpack_T

namespace internal_unpack_BN {
static FUNC_PREFIX void
unpack_BN(float *RESTRICT const _data_buffer, float *RESTRICT _data_j,
          int64_t const _size_j_0, int64_t const _size_j_1,
          int64_t const _size_j_2, int64_t const _stride_j_0,
          int64_t const _stride_j_1, int64_t const _stride_j_2,
          int64_t const _stride_j_3) {
  for (int64_t ctr_2 = 0; ctr_2 < _size_j_2; ctr_2 += 1) {
    float *RESTRICT _data_j_20_311 =
        _data_j + _stride_j_2 * ctr_2 + 11 * _stride_j_3;
    for (int64_t ctr_1 = 0; ctr_1 < _size_j_1; ctr_1 += 1) {
      float *RESTRICT _data_j_20_311_10 = _stride_j_1 * ctr_1 + _data_j_20_311;
      for (int64_t ctr_0 = 0; ctr_0 < _size_j_0; ctr_0 += 1) {
        _data_j_20_311_10[_stride_j_0 * ctr_0] =
            _data_buffer[_size_j_0 * _size_j_1 * ctr_2 + _size_j_0 * ctr_1 +
                         ctr_0];
      }
    }
  }
}
} // namespace internal_unpack_BN

namespace internal_unpack_N {
static FUNC_PREFIX void
unpack_N(float *RESTRICT const _data_buffer, float *RESTRICT _data_j,
         int64_t const _size_j_0, int64_t const _size_j_1,
         int64_t const _size_j_2, int64_t const _stride_j_0,
         int64_t const _stride_j_1, int64_t const _stride_j_2,
         int64_t const _stride_j_3) {
  for (int64_t ctr_2 = 0; ctr_2 < _size_j_2; ctr_2 += 1) {
    float *RESTRICT _data_j_20_311 =
        _data_j + _stride_j_2 * ctr_2 + 11 * _stride_j_3;
    float *RESTRICT _data_j_20_34 =
        _data_j + _stride_j_2 * ctr_2 + 4 * _stride_j_3;
    float *RESTRICT _data_j_20_312 =
        _data_j + _stride_j_2 * ctr_2 + 12 * _stride_j_3;
    for (int64_t ctr_1 = 0; ctr_1 < _size_j_1; ctr_1 += 1) {
      float *RESTRICT _data_j_20_311_10 = _stride_j_1 * ctr_1 + _data_j_20_311;
      float *RESTRICT _data_j_20_34_10 = _stride_j_1 * ctr_1 + _data_j_20_34;
      float *RESTRICT _data_j_20_312_10 = _stride_j_1 * ctr_1 + _data_j_20_312;
      for (int64_t ctr_0 = 0; ctr_0 < _size_j_0; ctr_0 += 1) {
        _data_j_20_311_10[_stride_j_0 * ctr_0] =
            _data_buffer[3 * _size_j_0 * _size_j_1 * ctr_2 +
                         3 * _size_j_0 * ctr_1 + 3 * ctr_0];
        _data_j_20_34_10[_stride_j_0 * ctr_0] =
            _data_buffer[3 * _size_j_0 * _size_j_1 * ctr_2 +
                         3 * _size_j_0 * ctr_1 + 3 * ctr_0 + 1];
        _data_j_20_312_10[_stride_j_0 * ctr_0] =
            _data_buffer[3 * _size_j_0 * _size_j_1 * ctr_2 +
                         3 * _size_j_0 * ctr_1 + 3 * ctr_0 + 2];
      }
    }
  }
}
} // namespace internal_unpack_N

namespace internal_unpack_TN {
static FUNC_PREFIX void
unpack_TN(float *RESTRICT const _data_buffer, float *RESTRICT _data_j,
          int64_t const _size_j_0, int64_t const _size_j_1,
          int64_t const _size_j_2, int64_t const _stride_j_0,
          int64_t const _stride_j_1, int64_t const _stride_j_2,
          int64_t const _stride_j_3) {
  for (int64_t ctr_2 = 0; ctr_2 < _size_j_2; ctr_2 += 1) {
    float *RESTRICT _data_j_20_312 =
        _data_j + _stride_j_2 * ctr_2 + 12 * _stride_j_3;
    for (int64_t ctr_1 = 0; ctr_1 < _size_j_1; ctr_1 += 1) {
      float *RESTRICT _data_j_20_312_10 = _stride_j_1 * ctr_1 + _data_j_20_312;
      for (int64_t ctr_0 = 0; ctr_0 < _size_j_0; ctr_0 += 1) {
        _data_j_20_312_10[_stride_j_0 * ctr_0] =
            _data_buffer[_size_j_0 * _size_j_1 * ctr_2 + _size_j_0 * ctr_1 +
                         ctr_0];
      }
    }
  }
}
} // namespace internal_unpack_TN

void DensityPackInfo_single_precision::pack(Direction dir,
                                            unsigned char *byte_buffer,
                                            IBlock *block) const {
  float *buffer = reinterpret_cast<float *>(byte_buffer);

  auto j = block->getData<field::GhostLayerField<float, 13>>(jID);

  CellInterval ci;
  j->getSliceBeforeGhostLayer(dir, ci, 1, false);

  switch (dir) {
  case stencil::BSW: {
    float *RESTRICT _data_buffer = buffer;
    WALBERLA_ASSERT_GREATER_EQUAL(ci.xMin(), -int_c(j->nrOfGhostLayers()));
    WALBERLA_ASSERT_GREATER_EQUAL(ci.yMin(), -int_c(j->nrOfGhostLayers()));
    WALBERLA_ASSERT_GREATER_EQUAL(ci.zMin(), -int_c(j->nrOfGhostLayers()));
    float *RESTRICT const _data_j =
        j->dataAt(ci.xMin(), ci.yMin(), ci.zMin(), 0);
    WALBERLA_ASSERT_GREATER_EQUAL(j->xSizeWithGhostLayer(),
                                  int64_t(cell_idx_c(ci.xSize()) + 0));
    const int64_t _size_j_0 = int64_t(cell_idx_c(ci.xSize()) + 0);
    WALBERLA_ASSERT_GREATER_EQUAL(j->ySizeWithGhostLayer(),
                                  int64_t(cell_idx_c(ci.ySize()) + 0));
    const int64_t _size_j_1 = int64_t(cell_idx_c(ci.ySize()) + 0);
    WALBERLA_ASSERT_GREATER_EQUAL(j->zSizeWithGhostLayer(),
                                  int64_t(cell_idx_c(ci.zSize()) + 0));
    const int64_t _size_j_2 = int64_t(cell_idx_c(ci.zSize()) + 0);
    const int64_t _stride_j_0 = int64_t(j->xStride());
    const int64_t _stride_j_1 = int64_t(j->yStride());
    const int64_t _stride_j_2 = int64_t(j->zStride());
    const int64_t _stride_j_3 = int64_t(1 * int64_t(j->fStride()));
    internal_pack_BSW::pack_BSW(_data_buffer, _data_j, _size_j_0, _size_j_1,
                                _size_j_2, _stride_j_0, _stride_j_1,
                                _stride_j_2, _stride_j_3);
    break;
  }

  case stencil::SW: {
    float *RESTRICT _data_buffer = buffer;
    WALBERLA_ASSERT_GREATER_EQUAL(ci.xMin(), -int_c(j->nrOfGhostLayers()));
    WALBERLA_ASSERT_GREATER_EQUAL(ci.yMin(), -int_c(j->nrOfGhostLayers()));
    WALBERLA_ASSERT_GREATER_EQUAL(ci.zMin(), -int_c(j->nrOfGhostLayers()));
    float *RESTRICT const _data_j =
        j->dataAt(ci.xMin(), ci.yMin(), ci.zMin(), 0);
    WALBERLA_ASSERT_GREATER_EQUAL(j->xSizeWithGhostLayer(),
                                  int64_t(cell_idx_c(ci.xSize()) + 0));
    const int64_t _size_j_0 = int64_t(cell_idx_c(ci.xSize()) + 0);
    WALBERLA_ASSERT_GREATER_EQUAL(j->ySizeWithGhostLayer(),
                                  int64_t(cell_idx_c(ci.ySize()) + 0));
    const int64_t _size_j_1 = int64_t(cell_idx_c(ci.ySize()) + 0);
    WALBERLA_ASSERT_GREATER_EQUAL(j->zSizeWithGhostLayer(),
                                  int64_t(cell_idx_c(ci.zSize()) + 0));
    const int64_t _size_j_2 = int64_t(cell_idx_c(ci.zSize()) + 0);
    const int64_t _stride_j_0 = int64_t(j->xStride());
    const int64_t _stride_j_1 = int64_t(j->yStride());
    const int64_t _stride_j_2 = int64_t(j->zStride());
    const int64_t _stride_j_3 = int64_t(1 * int64_t(j->fStride()));
    internal_pack_SW::pack_SW(_data_buffer, _data_j, _size_j_0, _size_j_1,
                              _size_j_2, _stride_j_0, _stride_j_1, _stride_j_2,
                              _stride_j_3);
    break;
  }

  case stencil::TSW: {
    float *RESTRICT _data_buffer = buffer;
    WALBERLA_ASSERT_GREATER_EQUAL(ci.xMin(), -int_c(j->nrOfGhostLayers()));
    WALBERLA_ASSERT_GREATER_EQUAL(ci.yMin(), -int_c(j->nrOfGhostLayers()));
    WALBERLA_ASSERT_GREATER_EQUAL(ci.zMin(), -int_c(j->nrOfGhostLayers()));
    float *RESTRICT const _data_j =
        j->dataAt(ci.xMin(), ci.yMin(), ci.zMin(), 0);
    WALBERLA_ASSERT_GREATER_EQUAL(j->xSizeWithGhostLayer(),
                                  int64_t(cell_idx_c(ci.xSize()) + 0));
    const int64_t _size_j_0 = int64_t(cell_idx_c(ci.xSize()) + 0);
    WALBERLA_ASSERT_GREATER_EQUAL(j->ySizeWithGhostLayer(),
                                  int64_t(cell_idx_c(ci.ySize()) + 0));
    const int64_t _size_j_1 = int64_t(cell_idx_c(ci.ySize()) + 0);
    WALBERLA_ASSERT_GREATER_EQUAL(j->zSizeWithGhostLayer(),
                                  int64_t(cell_idx_c(ci.zSize()) + 0));
    const int64_t _size_j_2 = int64_t(cell_idx_c(ci.zSize()) + 0);
    const int64_t _stride_j_0 = int64_t(j->xStride());
    const int64_t _stride_j_1 = int64_t(j->yStride());
    const int64_t _stride_j_2 = int64_t(j->zStride());
    const int64_t _stride_j_3 = int64_t(1 * int64_t(j->fStride()));
    internal_pack_TSW::pack_TSW(_data_buffer, _data_j, _size_j_0, _size_j_1,
                                _size_j_2, _stride_j_0, _stride_j_1,
                                _stride_j_2, _stride_j_3);
    break;
  }

  case stencil::BW: {
    float *RESTRICT _data_buffer = buffer;
    WALBERLA_ASSERT_GREATER_EQUAL(ci.xMin(), -int_c(j->nrOfGhostLayers()));
    WALBERLA_ASSERT_GREATER_EQUAL(ci.yMin(), -int_c(j->nrOfGhostLayers()));
    WALBERLA_ASSERT_GREATER_EQUAL(ci.zMin(), -int_c(j->nrOfGhostLayers()));
    float *RESTRICT const _data_j =
        j->dataAt(ci.xMin(), ci.yMin(), ci.zMin(), 0);
    WALBERLA_ASSERT_GREATER_EQUAL(j->xSizeWithGhostLayer(),
                                  int64_t(cell_idx_c(ci.xSize()) + 0));
    const int64_t _size_j_0 = int64_t(cell_idx_c(ci.xSize()) + 0);
    WALBERLA_ASSERT_GREATER_EQUAL(j->ySizeWithGhostLayer(),
                                  int64_t(cell_idx_c(ci.ySize()) + 0));
    const int64_t _size_j_1 = int64_t(cell_idx_c(ci.ySize()) + 0);
    WALBERLA_ASSERT_GREATER_EQUAL(j->zSizeWithGhostLayer(),
                                  int64_t(cell_idx_c(ci.zSize()) + 0));
    const int64_t _size_j_2 = int64_t(cell_idx_c(ci.zSize()) + 0);
    const int64_t _stride_j_0 = int64_t(j->xStride());
    const int64_t _stride_j_1 = int64_t(j->yStride());
    const int64_t _stride_j_2 = int64_t(j->zStride());
    const int64_t _stride_j_3 = int64_t(1 * int64_t(j->fStride()));
    internal_pack_BW::pack_BW(_data_buffer, _data_j, _size_j_0, _size_j_1,
                              _size_j_2, _stride_j_0, _stride_j_1, _stride_j_2,
                              _stride_j_3);
    break;
  }

  case stencil::W: {
    float *RESTRICT _data_buffer = buffer;
    WALBERLA_ASSERT_GREATER_EQUAL(ci.xMin(), -int_c(j->nrOfGhostLayers()));
    WALBERLA_ASSERT_GREATER_EQUAL(ci.yMin(), -int_c(j->nrOfGhostLayers()));
    WALBERLA_ASSERT_GREATER_EQUAL(ci.zMin(), -int_c(j->nrOfGhostLayers()));
    float *RESTRICT const _data_j =
        j->dataAt(ci.xMin(), ci.yMin(), ci.zMin(), 0);
    WALBERLA_ASSERT_GREATER_EQUAL(j->xSizeWithGhostLayer(),
                                  int64_t(cell_idx_c(ci.xSize()) + 0));
    const int64_t _size_j_0 = int64_t(cell_idx_c(ci.xSize()) + 0);
    WALBERLA_ASSERT_GREATER_EQUAL(j->ySizeWithGhostLayer(),
                                  int64_t(cell_idx_c(ci.ySize()) + 0));
    const int64_t _size_j_1 = int64_t(cell_idx_c(ci.ySize()) + 0);
    WALBERLA_ASSERT_GREATER_EQUAL(j->zSizeWithGhostLayer(),
                                  int64_t(cell_idx_c(ci.zSize()) + 0));
    const int64_t _size_j_2 = int64_t(cell_idx_c(ci.zSize()) + 0);
    const int64_t _stride_j_0 = int64_t(j->xStride());
    const int64_t _stride_j_1 = int64_t(j->yStride());
    const int64_t _stride_j_2 = int64_t(j->zStride());
    const int64_t _stride_j_3 = int64_t(1 * int64_t(j->fStride()));
    internal_pack_W::pack_W(_data_buffer, _data_j, _size_j_0, _size_j_1,
                            _size_j_2, _stride_j_0, _stride_j_1, _stride_j_2,
                            _stride_j_3);
    break;
  }

  case stencil::TW: {
    float *RESTRICT _data_buffer = buffer;
    WALBERLA_ASSERT_GREATER_EQUAL(ci.xMin(), -int_c(j->nrOfGhostLayers()));
    WALBERLA_ASSERT_GREATER_EQUAL(ci.yMin(), -int_c(j->nrOfGhostLayers()));
    WALBERLA_ASSERT_GREATER_EQUAL(ci.zMin(), -int_c(j->nrOfGhostLayers()));
    float *RESTRICT const _data_j =
        j->dataAt(ci.xMin(), ci.yMin(), ci.zMin(), 0);
    WALBERLA_ASSERT_GREATER_EQUAL(j->xSizeWithGhostLayer(),
                                  int64_t(cell_idx_c(ci.xSize()) + 0));
    const int64_t _size_j_0 = int64_t(cell_idx_c(ci.xSize()) + 0);
    WALBERLA_ASSERT_GREATER_EQUAL(j->ySizeWithGhostLayer(),
                                  int64_t(cell_idx_c(ci.ySize()) + 0));
    const int64_t _size_j_1 = int64_t(cell_idx_c(ci.ySize()) + 0);
    WALBERLA_ASSERT_GREATER_EQUAL(j->zSizeWithGhostLayer(),
                                  int64_t(cell_idx_c(ci.zSize()) + 0));
    const int64_t _size_j_2 = int64_t(cell_idx_c(ci.zSize()) + 0);
    const int64_t _stride_j_0 = int64_t(j->xStride());
    const int64_t _stride_j_1 = int64_t(j->yStride());
    const int64_t _stride_j_2 = int64_t(j->zStride());
    const int64_t _stride_j_3 = int64_t(1 * int64_t(j->fStride()));
    internal_pack_TW::pack_TW(_data_buffer, _data_j, _size_j_0, _size_j_1,
                              _size_j_2, _stride_j_0, _stride_j_1, _stride_j_2,
                              _stride_j_3);
    break;
  }

  case stencil::BNW: {
    float *RESTRICT _data_buffer = buffer;
    WALBERLA_ASSERT_GREATER_EQUAL(ci.xMin(), -int_c(j->nrOfGhostLayers()));
    WALBERLA_ASSERT_GREATER_EQUAL(ci.yMin(), -int_c(j->nrOfGhostLayers()));
    WALBERLA_ASSERT_GREATER_EQUAL(ci.zMin(), -int_c(j->nrOfGhostLayers()));
    float *RESTRICT const _data_j =
        j->dataAt(ci.xMin(), ci.yMin(), ci.zMin(), 0);
    WALBERLA_ASSERT_GREATER_EQUAL(j->xSizeWithGhostLayer(),
                                  int64_t(cell_idx_c(ci.xSize()) + 0));
    const int64_t _size_j_0 = int64_t(cell_idx_c(ci.xSize()) + 0);
    WALBERLA_ASSERT_GREATER_EQUAL(j->ySizeWithGhostLayer(),
                                  int64_t(cell_idx_c(ci.ySize()) + 0));
    const int64_t _size_j_1 = int64_t(cell_idx_c(ci.ySize()) + 0);
    WALBERLA_ASSERT_GREATER_EQUAL(j->zSizeWithGhostLayer(),
                                  int64_t(cell_idx_c(ci.zSize()) + 0));
    const int64_t _size_j_2 = int64_t(cell_idx_c(ci.zSize()) + 0);
    const int64_t _stride_j_0 = int64_t(j->xStride());
    const int64_t _stride_j_1 = int64_t(j->yStride());
    const int64_t _stride_j_2 = int64_t(j->zStride());
    const int64_t _stride_j_3 = int64_t(1 * int64_t(j->fStride()));
    internal_pack_BNW::pack_BNW(_data_buffer, _data_j, _size_j_0, _size_j_1,
                                _size_j_2, _stride_j_0, _stride_j_1,
                                _stride_j_2, _stride_j_3);
    break;
  }

  case stencil::NW: {
    float *RESTRICT _data_buffer = buffer;
    WALBERLA_ASSERT_GREATER_EQUAL(ci.xMin(), -int_c(j->nrOfGhostLayers()));
    WALBERLA_ASSERT_GREATER_EQUAL(ci.yMin(), -int_c(j->nrOfGhostLayers()));
    WALBERLA_ASSERT_GREATER_EQUAL(ci.zMin(), -int_c(j->nrOfGhostLayers()));
    float *RESTRICT const _data_j =
        j->dataAt(ci.xMin(), ci.yMin(), ci.zMin(), 0);
    WALBERLA_ASSERT_GREATER_EQUAL(j->xSizeWithGhostLayer(),
                                  int64_t(cell_idx_c(ci.xSize()) + 0));
    const int64_t _size_j_0 = int64_t(cell_idx_c(ci.xSize()) + 0);
    WALBERLA_ASSERT_GREATER_EQUAL(j->ySizeWithGhostLayer(),
                                  int64_t(cell_idx_c(ci.ySize()) + 0));
    const int64_t _size_j_1 = int64_t(cell_idx_c(ci.ySize()) + 0);
    WALBERLA_ASSERT_GREATER_EQUAL(j->zSizeWithGhostLayer(),
                                  int64_t(cell_idx_c(ci.zSize()) + 0));
    const int64_t _size_j_2 = int64_t(cell_idx_c(ci.zSize()) + 0);
    const int64_t _stride_j_0 = int64_t(j->xStride());
    const int64_t _stride_j_1 = int64_t(j->yStride());
    const int64_t _stride_j_2 = int64_t(j->zStride());
    const int64_t _stride_j_3 = int64_t(1 * int64_t(j->fStride()));
    internal_pack_NW::pack_NW(_data_buffer, _data_j, _size_j_0, _size_j_1,
                              _size_j_2, _stride_j_0, _stride_j_1, _stride_j_2,
                              _stride_j_3);
    break;
  }

  case stencil::TNW: {
    float *RESTRICT _data_buffer = buffer;
    WALBERLA_ASSERT_GREATER_EQUAL(ci.xMin(), -int_c(j->nrOfGhostLayers()));
    WALBERLA_ASSERT_GREATER_EQUAL(ci.yMin(), -int_c(j->nrOfGhostLayers()));
    WALBERLA_ASSERT_GREATER_EQUAL(ci.zMin(), -int_c(j->nrOfGhostLayers()));
    float *RESTRICT const _data_j =
        j->dataAt(ci.xMin(), ci.yMin(), ci.zMin(), 0);
    WALBERLA_ASSERT_GREATER_EQUAL(j->xSizeWithGhostLayer(),
                                  int64_t(cell_idx_c(ci.xSize()) + 0));
    const int64_t _size_j_0 = int64_t(cell_idx_c(ci.xSize()) + 0);
    WALBERLA_ASSERT_GREATER_EQUAL(j->ySizeWithGhostLayer(),
                                  int64_t(cell_idx_c(ci.ySize()) + 0));
    const int64_t _size_j_1 = int64_t(cell_idx_c(ci.ySize()) + 0);
    WALBERLA_ASSERT_GREATER_EQUAL(j->zSizeWithGhostLayer(),
                                  int64_t(cell_idx_c(ci.zSize()) + 0));
    const int64_t _size_j_2 = int64_t(cell_idx_c(ci.zSize()) + 0);
    const int64_t _stride_j_0 = int64_t(j->xStride());
    const int64_t _stride_j_1 = int64_t(j->yStride());
    const int64_t _stride_j_2 = int64_t(j->zStride());
    const int64_t _stride_j_3 = int64_t(1 * int64_t(j->fStride()));
    internal_pack_TNW::pack_TNW(_data_buffer, _data_j, _size_j_0, _size_j_1,
                                _size_j_2, _stride_j_0, _stride_j_1,
                                _stride_j_2, _stride_j_3);
    break;
  }

  case stencil::BS: {
    float *RESTRICT _data_buffer = buffer;
    WALBERLA_ASSERT_GREATER_EQUAL(ci.xMin(), -int_c(j->nrOfGhostLayers()));
    WALBERLA_ASSERT_GREATER_EQUAL(ci.yMin(), -int_c(j->nrOfGhostLayers()));
    WALBERLA_ASSERT_GREATER_EQUAL(ci.zMin(), -int_c(j->nrOfGhostLayers()));
    float *RESTRICT const _data_j =
        j->dataAt(ci.xMin(), ci.yMin(), ci.zMin(), 0);
    WALBERLA_ASSERT_GREATER_EQUAL(j->xSizeWithGhostLayer(),
                                  int64_t(cell_idx_c(ci.xSize()) + 0));
    const int64_t _size_j_0 = int64_t(cell_idx_c(ci.xSize()) + 0);
    WALBERLA_ASSERT_GREATER_EQUAL(j->ySizeWithGhostLayer(),
                                  int64_t(cell_idx_c(ci.ySize()) + 0));
    const int64_t _size_j_1 = int64_t(cell_idx_c(ci.ySize()) + 0);
    WALBERLA_ASSERT_GREATER_EQUAL(j->zSizeWithGhostLayer(),
                                  int64_t(cell_idx_c(ci.zSize()) + 0));
    const int64_t _size_j_2 = int64_t(cell_idx_c(ci.zSize()) + 0);
    const int64_t _stride_j_0 = int64_t(j->xStride());
    const int64_t _stride_j_1 = int64_t(j->yStride());
    const int64_t _stride_j_2 = int64_t(j->zStride());
    const int64_t _stride_j_3 = int64_t(1 * int64_t(j->fStride()));
    internal_pack_BS::pack_BS(_data_buffer, _data_j, _size_j_0, _size_j_1,
                              _size_j_2, _stride_j_0, _stride_j_1, _stride_j_2,
                              _stride_j_3);
    break;
  }

  case stencil::S: {
    float *RESTRICT _data_buffer = buffer;
    WALBERLA_ASSERT_GREATER_EQUAL(ci.xMin(), -int_c(j->nrOfGhostLayers()));
    WALBERLA_ASSERT_GREATER_EQUAL(ci.yMin(), -int_c(j->nrOfGhostLayers()));
    WALBERLA_ASSERT_GREATER_EQUAL(ci.zMin(), -int_c(j->nrOfGhostLayers()));
    float *RESTRICT const _data_j =
        j->dataAt(ci.xMin(), ci.yMin(), ci.zMin(), 0);
    WALBERLA_ASSERT_GREATER_EQUAL(j->xSizeWithGhostLayer(),
                                  int64_t(cell_idx_c(ci.xSize()) + 0));
    const int64_t _size_j_0 = int64_t(cell_idx_c(ci.xSize()) + 0);
    WALBERLA_ASSERT_GREATER_EQUAL(j->ySizeWithGhostLayer(),
                                  int64_t(cell_idx_c(ci.ySize()) + 0));
    const int64_t _size_j_1 = int64_t(cell_idx_c(ci.ySize()) + 0);
    WALBERLA_ASSERT_GREATER_EQUAL(j->zSizeWithGhostLayer(),
                                  int64_t(cell_idx_c(ci.zSize()) + 0));
    const int64_t _size_j_2 = int64_t(cell_idx_c(ci.zSize()) + 0);
    const int64_t _stride_j_0 = int64_t(j->xStride());
    const int64_t _stride_j_1 = int64_t(j->yStride());
    const int64_t _stride_j_2 = int64_t(j->zStride());
    const int64_t _stride_j_3 = int64_t(1 * int64_t(j->fStride()));
    internal_pack_S::pack_S(_data_buffer, _data_j, _size_j_0, _size_j_1,
                            _size_j_2, _stride_j_0, _stride_j_1, _stride_j_2,
                            _stride_j_3);
    break;
  }

  case stencil::TS: {
    float *RESTRICT _data_buffer = buffer;
    WALBERLA_ASSERT_GREATER_EQUAL(ci.xMin(), -int_c(j->nrOfGhostLayers()));
    WALBERLA_ASSERT_GREATER_EQUAL(ci.yMin(), -int_c(j->nrOfGhostLayers()));
    WALBERLA_ASSERT_GREATER_EQUAL(ci.zMin(), -int_c(j->nrOfGhostLayers()));
    float *RESTRICT const _data_j =
        j->dataAt(ci.xMin(), ci.yMin(), ci.zMin(), 0);
    WALBERLA_ASSERT_GREATER_EQUAL(j->xSizeWithGhostLayer(),
                                  int64_t(cell_idx_c(ci.xSize()) + 0));
    const int64_t _size_j_0 = int64_t(cell_idx_c(ci.xSize()) + 0);
    WALBERLA_ASSERT_GREATER_EQUAL(j->ySizeWithGhostLayer(),
                                  int64_t(cell_idx_c(ci.ySize()) + 0));
    const int64_t _size_j_1 = int64_t(cell_idx_c(ci.ySize()) + 0);
    WALBERLA_ASSERT_GREATER_EQUAL(j->zSizeWithGhostLayer(),
                                  int64_t(cell_idx_c(ci.zSize()) + 0));
    const int64_t _size_j_2 = int64_t(cell_idx_c(ci.zSize()) + 0);
    const int64_t _stride_j_0 = int64_t(j->xStride());
    const int64_t _stride_j_1 = int64_t(j->yStride());
    const int64_t _stride_j_2 = int64_t(j->zStride());
    const int64_t _stride_j_3 = int64_t(1 * int64_t(j->fStride()));
    internal_pack_TS::pack_TS(_data_buffer, _data_j, _size_j_0, _size_j_1,
                              _size_j_2, _stride_j_0, _stride_j_1, _stride_j_2,
                              _stride_j_3);
    break;
  }

  case stencil::B: {
    float *RESTRICT _data_buffer = buffer;
    WALBERLA_ASSERT_GREATER_EQUAL(ci.xMin(), -int_c(j->nrOfGhostLayers()));
    WALBERLA_ASSERT_GREATER_EQUAL(ci.yMin(), -int_c(j->nrOfGhostLayers()));
    WALBERLA_ASSERT_GREATER_EQUAL(ci.zMin(), -int_c(j->nrOfGhostLayers()));
    float *RESTRICT const _data_j =
        j->dataAt(ci.xMin(), ci.yMin(), ci.zMin(), 0);
    WALBERLA_ASSERT_GREATER_EQUAL(j->xSizeWithGhostLayer(),
                                  int64_t(cell_idx_c(ci.xSize()) + 0));
    const int64_t _size_j_0 = int64_t(cell_idx_c(ci.xSize()) + 0);
    WALBERLA_ASSERT_GREATER_EQUAL(j->ySizeWithGhostLayer(),
                                  int64_t(cell_idx_c(ci.ySize()) + 0));
    const int64_t _size_j_1 = int64_t(cell_idx_c(ci.ySize()) + 0);
    WALBERLA_ASSERT_GREATER_EQUAL(j->zSizeWithGhostLayer(),
                                  int64_t(cell_idx_c(ci.zSize()) + 0));
    const int64_t _size_j_2 = int64_t(cell_idx_c(ci.zSize()) + 0);
    const int64_t _stride_j_0 = int64_t(j->xStride());
    const int64_t _stride_j_1 = int64_t(j->yStride());
    const int64_t _stride_j_2 = int64_t(j->zStride());
    const int64_t _stride_j_3 = int64_t(1 * int64_t(j->fStride()));
    internal_pack_B::pack_B(_data_buffer, _data_j, _size_j_0, _size_j_1,
                            _size_j_2, _stride_j_0, _stride_j_1, _stride_j_2,
                            _stride_j_3);
    break;
  }

  case stencil::T: {
    float *RESTRICT _data_buffer = buffer;
    WALBERLA_ASSERT_GREATER_EQUAL(ci.xMin(), -int_c(j->nrOfGhostLayers()));
    WALBERLA_ASSERT_GREATER_EQUAL(ci.yMin(), -int_c(j->nrOfGhostLayers()));
    WALBERLA_ASSERT_GREATER_EQUAL(ci.zMin(), -int_c(j->nrOfGhostLayers()));
    float *RESTRICT const _data_j =
        j->dataAt(ci.xMin(), ci.yMin(), ci.zMin(), 0);
    WALBERLA_ASSERT_GREATER_EQUAL(j->xSizeWithGhostLayer(),
                                  int64_t(cell_idx_c(ci.xSize()) + 0));
    const int64_t _size_j_0 = int64_t(cell_idx_c(ci.xSize()) + 0);
    WALBERLA_ASSERT_GREATER_EQUAL(j->ySizeWithGhostLayer(),
                                  int64_t(cell_idx_c(ci.ySize()) + 0));
    const int64_t _size_j_1 = int64_t(cell_idx_c(ci.ySize()) + 0);
    WALBERLA_ASSERT_GREATER_EQUAL(j->zSizeWithGhostLayer(),
                                  int64_t(cell_idx_c(ci.zSize()) + 0));
    const int64_t _size_j_2 = int64_t(cell_idx_c(ci.zSize()) + 0);
    const int64_t _stride_j_0 = int64_t(j->xStride());
    const int64_t _stride_j_1 = int64_t(j->yStride());
    const int64_t _stride_j_2 = int64_t(j->zStride());
    const int64_t _stride_j_3 = int64_t(1 * int64_t(j->fStride()));
    internal_pack_T::pack_T(_data_buffer, _data_j, _size_j_0, _size_j_1,
                            _size_j_2, _stride_j_0, _stride_j_1, _stride_j_2,
                            _stride_j_3);
    break;
  }

  case stencil::BN: {
    float *RESTRICT _data_buffer = buffer;
    WALBERLA_ASSERT_GREATER_EQUAL(ci.xMin(), -int_c(j->nrOfGhostLayers()));
    WALBERLA_ASSERT_GREATER_EQUAL(ci.yMin(), -int_c(j->nrOfGhostLayers()));
    WALBERLA_ASSERT_GREATER_EQUAL(ci.zMin(), -int_c(j->nrOfGhostLayers()));
    float *RESTRICT const _data_j =
        j->dataAt(ci.xMin(), ci.yMin(), ci.zMin(), 0);
    WALBERLA_ASSERT_GREATER_EQUAL(j->xSizeWithGhostLayer(),
                                  int64_t(cell_idx_c(ci.xSize()) + 0));
    const int64_t _size_j_0 = int64_t(cell_idx_c(ci.xSize()) + 0);
    WALBERLA_ASSERT_GREATER_EQUAL(j->ySizeWithGhostLayer(),
                                  int64_t(cell_idx_c(ci.ySize()) + 0));
    const int64_t _size_j_1 = int64_t(cell_idx_c(ci.ySize()) + 0);
    WALBERLA_ASSERT_GREATER_EQUAL(j->zSizeWithGhostLayer(),
                                  int64_t(cell_idx_c(ci.zSize()) + 0));
    const int64_t _size_j_2 = int64_t(cell_idx_c(ci.zSize()) + 0);
    const int64_t _stride_j_0 = int64_t(j->xStride());
    const int64_t _stride_j_1 = int64_t(j->yStride());
    const int64_t _stride_j_2 = int64_t(j->zStride());
    const int64_t _stride_j_3 = int64_t(1 * int64_t(j->fStride()));
    internal_pack_BN::pack_BN(_data_buffer, _data_j, _size_j_0, _size_j_1,
                              _size_j_2, _stride_j_0, _stride_j_1, _stride_j_2,
                              _stride_j_3);
    break;
  }

  case stencil::N: {
    float *RESTRICT _data_buffer = buffer;
    WALBERLA_ASSERT_GREATER_EQUAL(ci.xMin(), -int_c(j->nrOfGhostLayers()));
    WALBERLA_ASSERT_GREATER_EQUAL(ci.yMin(), -int_c(j->nrOfGhostLayers()));
    WALBERLA_ASSERT_GREATER_EQUAL(ci.zMin(), -int_c(j->nrOfGhostLayers()));
    float *RESTRICT const _data_j =
        j->dataAt(ci.xMin(), ci.yMin(), ci.zMin(), 0);
    WALBERLA_ASSERT_GREATER_EQUAL(j->xSizeWithGhostLayer(),
                                  int64_t(cell_idx_c(ci.xSize()) + 0));
    const int64_t _size_j_0 = int64_t(cell_idx_c(ci.xSize()) + 0);
    WALBERLA_ASSERT_GREATER_EQUAL(j->ySizeWithGhostLayer(),
                                  int64_t(cell_idx_c(ci.ySize()) + 0));
    const int64_t _size_j_1 = int64_t(cell_idx_c(ci.ySize()) + 0);
    WALBERLA_ASSERT_GREATER_EQUAL(j->zSizeWithGhostLayer(),
                                  int64_t(cell_idx_c(ci.zSize()) + 0));
    const int64_t _size_j_2 = int64_t(cell_idx_c(ci.zSize()) + 0);
    const int64_t _stride_j_0 = int64_t(j->xStride());
    const int64_t _stride_j_1 = int64_t(j->yStride());
    const int64_t _stride_j_2 = int64_t(j->zStride());
    const int64_t _stride_j_3 = int64_t(1 * int64_t(j->fStride()));
    internal_pack_N::pack_N(_data_buffer, _data_j, _size_j_0, _size_j_1,
                            _size_j_2, _stride_j_0, _stride_j_1, _stride_j_2,
                            _stride_j_3);
    break;
  }

  case stencil::TN: {
    float *RESTRICT _data_buffer = buffer;
    WALBERLA_ASSERT_GREATER_EQUAL(ci.xMin(), -int_c(j->nrOfGhostLayers()));
    WALBERLA_ASSERT_GREATER_EQUAL(ci.yMin(), -int_c(j->nrOfGhostLayers()));
    WALBERLA_ASSERT_GREATER_EQUAL(ci.zMin(), -int_c(j->nrOfGhostLayers()));
    float *RESTRICT const _data_j =
        j->dataAt(ci.xMin(), ci.yMin(), ci.zMin(), 0);
    WALBERLA_ASSERT_GREATER_EQUAL(j->xSizeWithGhostLayer(),
                                  int64_t(cell_idx_c(ci.xSize()) + 0));
    const int64_t _size_j_0 = int64_t(cell_idx_c(ci.xSize()) + 0);
    WALBERLA_ASSERT_GREATER_EQUAL(j->ySizeWithGhostLayer(),
                                  int64_t(cell_idx_c(ci.ySize()) + 0));
    const int64_t _size_j_1 = int64_t(cell_idx_c(ci.ySize()) + 0);
    WALBERLA_ASSERT_GREATER_EQUAL(j->zSizeWithGhostLayer(),
                                  int64_t(cell_idx_c(ci.zSize()) + 0));
    const int64_t _size_j_2 = int64_t(cell_idx_c(ci.zSize()) + 0);
    const int64_t _stride_j_0 = int64_t(j->xStride());
    const int64_t _stride_j_1 = int64_t(j->yStride());
    const int64_t _stride_j_2 = int64_t(j->zStride());
    const int64_t _stride_j_3 = int64_t(1 * int64_t(j->fStride()));
    internal_pack_TN::pack_TN(_data_buffer, _data_j, _size_j_0, _size_j_1,
                              _size_j_2, _stride_j_0, _stride_j_1, _stride_j_2,
                              _stride_j_3);
    break;
  }

  default:
    WALBERLA_ASSERT(false);
  }
}

void DensityPackInfo_single_precision::unpack(Direction dir,
                                              unsigned char *byte_buffer,
                                              IBlock *block) const {
  float *buffer = reinterpret_cast<float *>(byte_buffer);

  auto j = block->getData<field::GhostLayerField<float, 13>>(jID);

  CellInterval ci;
  j->getGhostRegion(dir, ci, 1, false);
  auto communciationDirection = stencil::inverseDir[dir];

  switch (communciationDirection) {
  case stencil::BSW: {
    float *RESTRICT const _data_buffer = buffer;
    WALBERLA_ASSERT_GREATER_EQUAL(ci.xMin(), -int_c(j->nrOfGhostLayers()));
    WALBERLA_ASSERT_GREATER_EQUAL(ci.yMin(), -int_c(j->nrOfGhostLayers()));
    WALBERLA_ASSERT_GREATER_EQUAL(ci.zMin(), -int_c(j->nrOfGhostLayers()));
    float *RESTRICT _data_j = j->dataAt(ci.xMin(), ci.yMin(), ci.zMin(), 0);
    WALBERLA_ASSERT_GREATER_EQUAL(j->xSizeWithGhostLayer(),
                                  int64_t(cell_idx_c(ci.xSize()) + 0));
    const int64_t _size_j_0 = int64_t(cell_idx_c(ci.xSize()) + 0);
    WALBERLA_ASSERT_GREATER_EQUAL(j->ySizeWithGhostLayer(),
                                  int64_t(cell_idx_c(ci.ySize()) + 0));
    const int64_t _size_j_1 = int64_t(cell_idx_c(ci.ySize()) + 0);
    WALBERLA_ASSERT_GREATER_EQUAL(j->zSizeWithGhostLayer(),
                                  int64_t(cell_idx_c(ci.zSize()) + 0));
    const int64_t _size_j_2 = int64_t(cell_idx_c(ci.zSize()) + 0);
    const int64_t _stride_j_0 = int64_t(j->xStride());
    const int64_t _stride_j_1 = int64_t(j->yStride());
    const int64_t _stride_j_2 = int64_t(j->zStride());
    const int64_t _stride_j_3 = int64_t(1 * int64_t(j->fStride()));
    internal_unpack_BSW::unpack_BSW(_data_buffer, _data_j, _size_j_0, _size_j_1,
                                    _size_j_2, _stride_j_0, _stride_j_1,
                                    _stride_j_2, _stride_j_3);
    break;
  }

  case stencil::SW: {
    float *RESTRICT const _data_buffer = buffer;
    WALBERLA_ASSERT_GREATER_EQUAL(ci.xMin(), -int_c(j->nrOfGhostLayers()));
    WALBERLA_ASSERT_GREATER_EQUAL(ci.yMin(), -int_c(j->nrOfGhostLayers()));
    WALBERLA_ASSERT_GREATER_EQUAL(ci.zMin(), -int_c(j->nrOfGhostLayers()));
    float *RESTRICT _data_j = j->dataAt(ci.xMin(), ci.yMin(), ci.zMin(), 0);
    WALBERLA_ASSERT_GREATER_EQUAL(j->xSizeWithGhostLayer(),
                                  int64_t(cell_idx_c(ci.xSize()) + 0));
    const int64_t _size_j_0 = int64_t(cell_idx_c(ci.xSize()) + 0);
    WALBERLA_ASSERT_GREATER_EQUAL(j->ySizeWithGhostLayer(),
                                  int64_t(cell_idx_c(ci.ySize()) + 0));
    const int64_t _size_j_1 = int64_t(cell_idx_c(ci.ySize()) + 0);
    WALBERLA_ASSERT_GREATER_EQUAL(j->zSizeWithGhostLayer(),
                                  int64_t(cell_idx_c(ci.zSize()) + 0));
    const int64_t _size_j_2 = int64_t(cell_idx_c(ci.zSize()) + 0);
    const int64_t _stride_j_0 = int64_t(j->xStride());
    const int64_t _stride_j_1 = int64_t(j->yStride());
    const int64_t _stride_j_2 = int64_t(j->zStride());
    const int64_t _stride_j_3 = int64_t(1 * int64_t(j->fStride()));
    internal_unpack_SW::unpack_SW(_data_buffer, _data_j, _size_j_0, _size_j_1,
                                  _size_j_2, _stride_j_0, _stride_j_1,
                                  _stride_j_2, _stride_j_3);
    break;
  }

  case stencil::TSW: {
    float *RESTRICT const _data_buffer = buffer;
    WALBERLA_ASSERT_GREATER_EQUAL(ci.xMin(), -int_c(j->nrOfGhostLayers()));
    WALBERLA_ASSERT_GREATER_EQUAL(ci.yMin(), -int_c(j->nrOfGhostLayers()));
    WALBERLA_ASSERT_GREATER_EQUAL(ci.zMin(), -int_c(j->nrOfGhostLayers()));
    float *RESTRICT _data_j = j->dataAt(ci.xMin(), ci.yMin(), ci.zMin(), 0);
    WALBERLA_ASSERT_GREATER_EQUAL(j->xSizeWithGhostLayer(),
                                  int64_t(cell_idx_c(ci.xSize()) + 0));
    const int64_t _size_j_0 = int64_t(cell_idx_c(ci.xSize()) + 0);
    WALBERLA_ASSERT_GREATER_EQUAL(j->ySizeWithGhostLayer(),
                                  int64_t(cell_idx_c(ci.ySize()) + 0));
    const int64_t _size_j_1 = int64_t(cell_idx_c(ci.ySize()) + 0);
    WALBERLA_ASSERT_GREATER_EQUAL(j->zSizeWithGhostLayer(),
                                  int64_t(cell_idx_c(ci.zSize()) + 0));
    const int64_t _size_j_2 = int64_t(cell_idx_c(ci.zSize()) + 0);
    const int64_t _stride_j_0 = int64_t(j->xStride());
    const int64_t _stride_j_1 = int64_t(j->yStride());
    const int64_t _stride_j_2 = int64_t(j->zStride());
    const int64_t _stride_j_3 = int64_t(1 * int64_t(j->fStride()));
    internal_unpack_TSW::unpack_TSW(_data_buffer, _data_j, _size_j_0, _size_j_1,
                                    _size_j_2, _stride_j_0, _stride_j_1,
                                    _stride_j_2, _stride_j_3);
    break;
  }

  case stencil::BW: {
    float *RESTRICT const _data_buffer = buffer;
    WALBERLA_ASSERT_GREATER_EQUAL(ci.xMin(), -int_c(j->nrOfGhostLayers()));
    WALBERLA_ASSERT_GREATER_EQUAL(ci.yMin(), -int_c(j->nrOfGhostLayers()));
    WALBERLA_ASSERT_GREATER_EQUAL(ci.zMin(), -int_c(j->nrOfGhostLayers()));
    float *RESTRICT _data_j = j->dataAt(ci.xMin(), ci.yMin(), ci.zMin(), 0);
    WALBERLA_ASSERT_GREATER_EQUAL(j->xSizeWithGhostLayer(),
                                  int64_t(cell_idx_c(ci.xSize()) + 0));
    const int64_t _size_j_0 = int64_t(cell_idx_c(ci.xSize()) + 0);
    WALBERLA_ASSERT_GREATER_EQUAL(j->ySizeWithGhostLayer(),
                                  int64_t(cell_idx_c(ci.ySize()) + 0));
    const int64_t _size_j_1 = int64_t(cell_idx_c(ci.ySize()) + 0);
    WALBERLA_ASSERT_GREATER_EQUAL(j->zSizeWithGhostLayer(),
                                  int64_t(cell_idx_c(ci.zSize()) + 0));
    const int64_t _size_j_2 = int64_t(cell_idx_c(ci.zSize()) + 0);
    const int64_t _stride_j_0 = int64_t(j->xStride());
    const int64_t _stride_j_1 = int64_t(j->yStride());
    const int64_t _stride_j_2 = int64_t(j->zStride());
    const int64_t _stride_j_3 = int64_t(1 * int64_t(j->fStride()));
    internal_unpack_BW::unpack_BW(_data_buffer, _data_j, _size_j_0, _size_j_1,
                                  _size_j_2, _stride_j_0, _stride_j_1,
                                  _stride_j_2, _stride_j_3);
    break;
  }

  case stencil::W: {
    float *RESTRICT const _data_buffer = buffer;
    WALBERLA_ASSERT_GREATER_EQUAL(ci.xMin(), -int_c(j->nrOfGhostLayers()));
    WALBERLA_ASSERT_GREATER_EQUAL(ci.yMin(), -int_c(j->nrOfGhostLayers()));
    WALBERLA_ASSERT_GREATER_EQUAL(ci.zMin(), -int_c(j->nrOfGhostLayers()));
    float *RESTRICT _data_j = j->dataAt(ci.xMin(), ci.yMin(), ci.zMin(), 0);
    WALBERLA_ASSERT_GREATER_EQUAL(j->xSizeWithGhostLayer(),
                                  int64_t(cell_idx_c(ci.xSize()) + 0));
    const int64_t _size_j_0 = int64_t(cell_idx_c(ci.xSize()) + 0);
    WALBERLA_ASSERT_GREATER_EQUAL(j->ySizeWithGhostLayer(),
                                  int64_t(cell_idx_c(ci.ySize()) + 0));
    const int64_t _size_j_1 = int64_t(cell_idx_c(ci.ySize()) + 0);
    WALBERLA_ASSERT_GREATER_EQUAL(j->zSizeWithGhostLayer(),
                                  int64_t(cell_idx_c(ci.zSize()) + 0));
    const int64_t _size_j_2 = int64_t(cell_idx_c(ci.zSize()) + 0);
    const int64_t _stride_j_0 = int64_t(j->xStride());
    const int64_t _stride_j_1 = int64_t(j->yStride());
    const int64_t _stride_j_2 = int64_t(j->zStride());
    const int64_t _stride_j_3 = int64_t(1 * int64_t(j->fStride()));
    internal_unpack_W::unpack_W(_data_buffer, _data_j, _size_j_0, _size_j_1,
                                _size_j_2, _stride_j_0, _stride_j_1,
                                _stride_j_2, _stride_j_3);
    break;
  }

  case stencil::TW: {
    float *RESTRICT const _data_buffer = buffer;
    WALBERLA_ASSERT_GREATER_EQUAL(ci.xMin(), -int_c(j->nrOfGhostLayers()));
    WALBERLA_ASSERT_GREATER_EQUAL(ci.yMin(), -int_c(j->nrOfGhostLayers()));
    WALBERLA_ASSERT_GREATER_EQUAL(ci.zMin(), -int_c(j->nrOfGhostLayers()));
    float *RESTRICT _data_j = j->dataAt(ci.xMin(), ci.yMin(), ci.zMin(), 0);
    WALBERLA_ASSERT_GREATER_EQUAL(j->xSizeWithGhostLayer(),
                                  int64_t(cell_idx_c(ci.xSize()) + 0));
    const int64_t _size_j_0 = int64_t(cell_idx_c(ci.xSize()) + 0);
    WALBERLA_ASSERT_GREATER_EQUAL(j->ySizeWithGhostLayer(),
                                  int64_t(cell_idx_c(ci.ySize()) + 0));
    const int64_t _size_j_1 = int64_t(cell_idx_c(ci.ySize()) + 0);
    WALBERLA_ASSERT_GREATER_EQUAL(j->zSizeWithGhostLayer(),
                                  int64_t(cell_idx_c(ci.zSize()) + 0));
    const int64_t _size_j_2 = int64_t(cell_idx_c(ci.zSize()) + 0);
    const int64_t _stride_j_0 = int64_t(j->xStride());
    const int64_t _stride_j_1 = int64_t(j->yStride());
    const int64_t _stride_j_2 = int64_t(j->zStride());
    const int64_t _stride_j_3 = int64_t(1 * int64_t(j->fStride()));
    internal_unpack_TW::unpack_TW(_data_buffer, _data_j, _size_j_0, _size_j_1,
                                  _size_j_2, _stride_j_0, _stride_j_1,
                                  _stride_j_2, _stride_j_3);
    break;
  }

  case stencil::BNW: {
    float *RESTRICT const _data_buffer = buffer;
    WALBERLA_ASSERT_GREATER_EQUAL(ci.xMin(), -int_c(j->nrOfGhostLayers()));
    WALBERLA_ASSERT_GREATER_EQUAL(ci.yMin(), -int_c(j->nrOfGhostLayers()));
    WALBERLA_ASSERT_GREATER_EQUAL(ci.zMin(), -int_c(j->nrOfGhostLayers()));
    float *RESTRICT _data_j = j->dataAt(ci.xMin(), ci.yMin(), ci.zMin(), 0);
    WALBERLA_ASSERT_GREATER_EQUAL(j->xSizeWithGhostLayer(),
                                  int64_t(cell_idx_c(ci.xSize()) + 0));
    const int64_t _size_j_0 = int64_t(cell_idx_c(ci.xSize()) + 0);
    WALBERLA_ASSERT_GREATER_EQUAL(j->ySizeWithGhostLayer(),
                                  int64_t(cell_idx_c(ci.ySize()) + 0));
    const int64_t _size_j_1 = int64_t(cell_idx_c(ci.ySize()) + 0);
    WALBERLA_ASSERT_GREATER_EQUAL(j->zSizeWithGhostLayer(),
                                  int64_t(cell_idx_c(ci.zSize()) + 0));
    const int64_t _size_j_2 = int64_t(cell_idx_c(ci.zSize()) + 0);
    const int64_t _stride_j_0 = int64_t(j->xStride());
    const int64_t _stride_j_1 = int64_t(j->yStride());
    const int64_t _stride_j_2 = int64_t(j->zStride());
    const int64_t _stride_j_3 = int64_t(1 * int64_t(j->fStride()));
    internal_unpack_BNW::unpack_BNW(_data_buffer, _data_j, _size_j_0, _size_j_1,
                                    _size_j_2, _stride_j_0, _stride_j_1,
                                    _stride_j_2, _stride_j_3);
    break;
  }

  case stencil::NW: {
    float *RESTRICT const _data_buffer = buffer;
    WALBERLA_ASSERT_GREATER_EQUAL(ci.xMin(), -int_c(j->nrOfGhostLayers()));
    WALBERLA_ASSERT_GREATER_EQUAL(ci.yMin(), -int_c(j->nrOfGhostLayers()));
    WALBERLA_ASSERT_GREATER_EQUAL(ci.zMin(), -int_c(j->nrOfGhostLayers()));
    float *RESTRICT _data_j = j->dataAt(ci.xMin(), ci.yMin(), ci.zMin(), 0);
    WALBERLA_ASSERT_GREATER_EQUAL(j->xSizeWithGhostLayer(),
                                  int64_t(cell_idx_c(ci.xSize()) + 0));
    const int64_t _size_j_0 = int64_t(cell_idx_c(ci.xSize()) + 0);
    WALBERLA_ASSERT_GREATER_EQUAL(j->ySizeWithGhostLayer(),
                                  int64_t(cell_idx_c(ci.ySize()) + 0));
    const int64_t _size_j_1 = int64_t(cell_idx_c(ci.ySize()) + 0);
    WALBERLA_ASSERT_GREATER_EQUAL(j->zSizeWithGhostLayer(),
                                  int64_t(cell_idx_c(ci.zSize()) + 0));
    const int64_t _size_j_2 = int64_t(cell_idx_c(ci.zSize()) + 0);
    const int64_t _stride_j_0 = int64_t(j->xStride());
    const int64_t _stride_j_1 = int64_t(j->yStride());
    const int64_t _stride_j_2 = int64_t(j->zStride());
    const int64_t _stride_j_3 = int64_t(1 * int64_t(j->fStride()));
    internal_unpack_NW::unpack_NW(_data_buffer, _data_j, _size_j_0, _size_j_1,
                                  _size_j_2, _stride_j_0, _stride_j_1,
                                  _stride_j_2, _stride_j_3);
    break;
  }

  case stencil::TNW: {
    float *RESTRICT const _data_buffer = buffer;
    WALBERLA_ASSERT_GREATER_EQUAL(ci.xMin(), -int_c(j->nrOfGhostLayers()));
    WALBERLA_ASSERT_GREATER_EQUAL(ci.yMin(), -int_c(j->nrOfGhostLayers()));
    WALBERLA_ASSERT_GREATER_EQUAL(ci.zMin(), -int_c(j->nrOfGhostLayers()));
    float *RESTRICT _data_j = j->dataAt(ci.xMin(), ci.yMin(), ci.zMin(), 0);
    WALBERLA_ASSERT_GREATER_EQUAL(j->xSizeWithGhostLayer(),
                                  int64_t(cell_idx_c(ci.xSize()) + 0));
    const int64_t _size_j_0 = int64_t(cell_idx_c(ci.xSize()) + 0);
    WALBERLA_ASSERT_GREATER_EQUAL(j->ySizeWithGhostLayer(),
                                  int64_t(cell_idx_c(ci.ySize()) + 0));
    const int64_t _size_j_1 = int64_t(cell_idx_c(ci.ySize()) + 0);
    WALBERLA_ASSERT_GREATER_EQUAL(j->zSizeWithGhostLayer(),
                                  int64_t(cell_idx_c(ci.zSize()) + 0));
    const int64_t _size_j_2 = int64_t(cell_idx_c(ci.zSize()) + 0);
    const int64_t _stride_j_0 = int64_t(j->xStride());
    const int64_t _stride_j_1 = int64_t(j->yStride());
    const int64_t _stride_j_2 = int64_t(j->zStride());
    const int64_t _stride_j_3 = int64_t(1 * int64_t(j->fStride()));
    internal_unpack_TNW::unpack_TNW(_data_buffer, _data_j, _size_j_0, _size_j_1,
                                    _size_j_2, _stride_j_0, _stride_j_1,
                                    _stride_j_2, _stride_j_3);
    break;
  }

  case stencil::BS: {
    float *RESTRICT const _data_buffer = buffer;
    WALBERLA_ASSERT_GREATER_EQUAL(ci.xMin(), -int_c(j->nrOfGhostLayers()));
    WALBERLA_ASSERT_GREATER_EQUAL(ci.yMin(), -int_c(j->nrOfGhostLayers()));
    WALBERLA_ASSERT_GREATER_EQUAL(ci.zMin(), -int_c(j->nrOfGhostLayers()));
    float *RESTRICT _data_j = j->dataAt(ci.xMin(), ci.yMin(), ci.zMin(), 0);
    WALBERLA_ASSERT_GREATER_EQUAL(j->xSizeWithGhostLayer(),
                                  int64_t(cell_idx_c(ci.xSize()) + 0));
    const int64_t _size_j_0 = int64_t(cell_idx_c(ci.xSize()) + 0);
    WALBERLA_ASSERT_GREATER_EQUAL(j->ySizeWithGhostLayer(),
                                  int64_t(cell_idx_c(ci.ySize()) + 0));
    const int64_t _size_j_1 = int64_t(cell_idx_c(ci.ySize()) + 0);
    WALBERLA_ASSERT_GREATER_EQUAL(j->zSizeWithGhostLayer(),
                                  int64_t(cell_idx_c(ci.zSize()) + 0));
    const int64_t _size_j_2 = int64_t(cell_idx_c(ci.zSize()) + 0);
    const int64_t _stride_j_0 = int64_t(j->xStride());
    const int64_t _stride_j_1 = int64_t(j->yStride());
    const int64_t _stride_j_2 = int64_t(j->zStride());
    const int64_t _stride_j_3 = int64_t(1 * int64_t(j->fStride()));
    internal_unpack_BS::unpack_BS(_data_buffer, _data_j, _size_j_0, _size_j_1,
                                  _size_j_2, _stride_j_0, _stride_j_1,
                                  _stride_j_2, _stride_j_3);
    break;
  }

  case stencil::S: {
    float *RESTRICT const _data_buffer = buffer;
    WALBERLA_ASSERT_GREATER_EQUAL(ci.xMin(), -int_c(j->nrOfGhostLayers()));
    WALBERLA_ASSERT_GREATER_EQUAL(ci.yMin(), -int_c(j->nrOfGhostLayers()));
    WALBERLA_ASSERT_GREATER_EQUAL(ci.zMin(), -int_c(j->nrOfGhostLayers()));
    float *RESTRICT _data_j = j->dataAt(ci.xMin(), ci.yMin(), ci.zMin(), 0);
    WALBERLA_ASSERT_GREATER_EQUAL(j->xSizeWithGhostLayer(),
                                  int64_t(cell_idx_c(ci.xSize()) + 0));
    const int64_t _size_j_0 = int64_t(cell_idx_c(ci.xSize()) + 0);
    WALBERLA_ASSERT_GREATER_EQUAL(j->ySizeWithGhostLayer(),
                                  int64_t(cell_idx_c(ci.ySize()) + 0));
    const int64_t _size_j_1 = int64_t(cell_idx_c(ci.ySize()) + 0);
    WALBERLA_ASSERT_GREATER_EQUAL(j->zSizeWithGhostLayer(),
                                  int64_t(cell_idx_c(ci.zSize()) + 0));
    const int64_t _size_j_2 = int64_t(cell_idx_c(ci.zSize()) + 0);
    const int64_t _stride_j_0 = int64_t(j->xStride());
    const int64_t _stride_j_1 = int64_t(j->yStride());
    const int64_t _stride_j_2 = int64_t(j->zStride());
    const int64_t _stride_j_3 = int64_t(1 * int64_t(j->fStride()));
    internal_unpack_S::unpack_S(_data_buffer, _data_j, _size_j_0, _size_j_1,
                                _size_j_2, _stride_j_0, _stride_j_1,
                                _stride_j_2, _stride_j_3);
    break;
  }

  case stencil::TS: {
    float *RESTRICT const _data_buffer = buffer;
    WALBERLA_ASSERT_GREATER_EQUAL(ci.xMin(), -int_c(j->nrOfGhostLayers()));
    WALBERLA_ASSERT_GREATER_EQUAL(ci.yMin(), -int_c(j->nrOfGhostLayers()));
    WALBERLA_ASSERT_GREATER_EQUAL(ci.zMin(), -int_c(j->nrOfGhostLayers()));
    float *RESTRICT _data_j = j->dataAt(ci.xMin(), ci.yMin(), ci.zMin(), 0);
    WALBERLA_ASSERT_GREATER_EQUAL(j->xSizeWithGhostLayer(),
                                  int64_t(cell_idx_c(ci.xSize()) + 0));
    const int64_t _size_j_0 = int64_t(cell_idx_c(ci.xSize()) + 0);
    WALBERLA_ASSERT_GREATER_EQUAL(j->ySizeWithGhostLayer(),
                                  int64_t(cell_idx_c(ci.ySize()) + 0));
    const int64_t _size_j_1 = int64_t(cell_idx_c(ci.ySize()) + 0);
    WALBERLA_ASSERT_GREATER_EQUAL(j->zSizeWithGhostLayer(),
                                  int64_t(cell_idx_c(ci.zSize()) + 0));
    const int64_t _size_j_2 = int64_t(cell_idx_c(ci.zSize()) + 0);
    const int64_t _stride_j_0 = int64_t(j->xStride());
    const int64_t _stride_j_1 = int64_t(j->yStride());
    const int64_t _stride_j_2 = int64_t(j->zStride());
    const int64_t _stride_j_3 = int64_t(1 * int64_t(j->fStride()));
    internal_unpack_TS::unpack_TS(_data_buffer, _data_j, _size_j_0, _size_j_1,
                                  _size_j_2, _stride_j_0, _stride_j_1,
                                  _stride_j_2, _stride_j_3);
    break;
  }

  case stencil::B: {
    float *RESTRICT const _data_buffer = buffer;
    WALBERLA_ASSERT_GREATER_EQUAL(ci.xMin(), -int_c(j->nrOfGhostLayers()));
    WALBERLA_ASSERT_GREATER_EQUAL(ci.yMin(), -int_c(j->nrOfGhostLayers()));
    WALBERLA_ASSERT_GREATER_EQUAL(ci.zMin(), -int_c(j->nrOfGhostLayers()));
    float *RESTRICT _data_j = j->dataAt(ci.xMin(), ci.yMin(), ci.zMin(), 0);
    WALBERLA_ASSERT_GREATER_EQUAL(j->xSizeWithGhostLayer(),
                                  int64_t(cell_idx_c(ci.xSize()) + 0));
    const int64_t _size_j_0 = int64_t(cell_idx_c(ci.xSize()) + 0);
    WALBERLA_ASSERT_GREATER_EQUAL(j->ySizeWithGhostLayer(),
                                  int64_t(cell_idx_c(ci.ySize()) + 0));
    const int64_t _size_j_1 = int64_t(cell_idx_c(ci.ySize()) + 0);
    WALBERLA_ASSERT_GREATER_EQUAL(j->zSizeWithGhostLayer(),
                                  int64_t(cell_idx_c(ci.zSize()) + 0));
    const int64_t _size_j_2 = int64_t(cell_idx_c(ci.zSize()) + 0);
    const int64_t _stride_j_0 = int64_t(j->xStride());
    const int64_t _stride_j_1 = int64_t(j->yStride());
    const int64_t _stride_j_2 = int64_t(j->zStride());
    const int64_t _stride_j_3 = int64_t(1 * int64_t(j->fStride()));
    internal_unpack_B::unpack_B(_data_buffer, _data_j, _size_j_0, _size_j_1,
                                _size_j_2, _stride_j_0, _stride_j_1,
                                _stride_j_2, _stride_j_3);
    break;
  }

  case stencil::T: {
    float *RESTRICT const _data_buffer = buffer;
    WALBERLA_ASSERT_GREATER_EQUAL(ci.xMin(), -int_c(j->nrOfGhostLayers()));
    WALBERLA_ASSERT_GREATER_EQUAL(ci.yMin(), -int_c(j->nrOfGhostLayers()));
    WALBERLA_ASSERT_GREATER_EQUAL(ci.zMin(), -int_c(j->nrOfGhostLayers()));
    float *RESTRICT _data_j = j->dataAt(ci.xMin(), ci.yMin(), ci.zMin(), 0);
    WALBERLA_ASSERT_GREATER_EQUAL(j->xSizeWithGhostLayer(),
                                  int64_t(cell_idx_c(ci.xSize()) + 0));
    const int64_t _size_j_0 = int64_t(cell_idx_c(ci.xSize()) + 0);
    WALBERLA_ASSERT_GREATER_EQUAL(j->ySizeWithGhostLayer(),
                                  int64_t(cell_idx_c(ci.ySize()) + 0));
    const int64_t _size_j_1 = int64_t(cell_idx_c(ci.ySize()) + 0);
    WALBERLA_ASSERT_GREATER_EQUAL(j->zSizeWithGhostLayer(),
                                  int64_t(cell_idx_c(ci.zSize()) + 0));
    const int64_t _size_j_2 = int64_t(cell_idx_c(ci.zSize()) + 0);
    const int64_t _stride_j_0 = int64_t(j->xStride());
    const int64_t _stride_j_1 = int64_t(j->yStride());
    const int64_t _stride_j_2 = int64_t(j->zStride());
    const int64_t _stride_j_3 = int64_t(1 * int64_t(j->fStride()));
    internal_unpack_T::unpack_T(_data_buffer, _data_j, _size_j_0, _size_j_1,
                                _size_j_2, _stride_j_0, _stride_j_1,
                                _stride_j_2, _stride_j_3);
    break;
  }

  case stencil::BN: {
    float *RESTRICT const _data_buffer = buffer;
    WALBERLA_ASSERT_GREATER_EQUAL(ci.xMin(), -int_c(j->nrOfGhostLayers()));
    WALBERLA_ASSERT_GREATER_EQUAL(ci.yMin(), -int_c(j->nrOfGhostLayers()));
    WALBERLA_ASSERT_GREATER_EQUAL(ci.zMin(), -int_c(j->nrOfGhostLayers()));
    float *RESTRICT _data_j = j->dataAt(ci.xMin(), ci.yMin(), ci.zMin(), 0);
    WALBERLA_ASSERT_GREATER_EQUAL(j->xSizeWithGhostLayer(),
                                  int64_t(cell_idx_c(ci.xSize()) + 0));
    const int64_t _size_j_0 = int64_t(cell_idx_c(ci.xSize()) + 0);
    WALBERLA_ASSERT_GREATER_EQUAL(j->ySizeWithGhostLayer(),
                                  int64_t(cell_idx_c(ci.ySize()) + 0));
    const int64_t _size_j_1 = int64_t(cell_idx_c(ci.ySize()) + 0);
    WALBERLA_ASSERT_GREATER_EQUAL(j->zSizeWithGhostLayer(),
                                  int64_t(cell_idx_c(ci.zSize()) + 0));
    const int64_t _size_j_2 = int64_t(cell_idx_c(ci.zSize()) + 0);
    const int64_t _stride_j_0 = int64_t(j->xStride());
    const int64_t _stride_j_1 = int64_t(j->yStride());
    const int64_t _stride_j_2 = int64_t(j->zStride());
    const int64_t _stride_j_3 = int64_t(1 * int64_t(j->fStride()));
    internal_unpack_BN::unpack_BN(_data_buffer, _data_j, _size_j_0, _size_j_1,
                                  _size_j_2, _stride_j_0, _stride_j_1,
                                  _stride_j_2, _stride_j_3);
    break;
  }

  case stencil::N: {
    float *RESTRICT const _data_buffer = buffer;
    WALBERLA_ASSERT_GREATER_EQUAL(ci.xMin(), -int_c(j->nrOfGhostLayers()));
    WALBERLA_ASSERT_GREATER_EQUAL(ci.yMin(), -int_c(j->nrOfGhostLayers()));
    WALBERLA_ASSERT_GREATER_EQUAL(ci.zMin(), -int_c(j->nrOfGhostLayers()));
    float *RESTRICT _data_j = j->dataAt(ci.xMin(), ci.yMin(), ci.zMin(), 0);
    WALBERLA_ASSERT_GREATER_EQUAL(j->xSizeWithGhostLayer(),
                                  int64_t(cell_idx_c(ci.xSize()) + 0));
    const int64_t _size_j_0 = int64_t(cell_idx_c(ci.xSize()) + 0);
    WALBERLA_ASSERT_GREATER_EQUAL(j->ySizeWithGhostLayer(),
                                  int64_t(cell_idx_c(ci.ySize()) + 0));
    const int64_t _size_j_1 = int64_t(cell_idx_c(ci.ySize()) + 0);
    WALBERLA_ASSERT_GREATER_EQUAL(j->zSizeWithGhostLayer(),
                                  int64_t(cell_idx_c(ci.zSize()) + 0));
    const int64_t _size_j_2 = int64_t(cell_idx_c(ci.zSize()) + 0);
    const int64_t _stride_j_0 = int64_t(j->xStride());
    const int64_t _stride_j_1 = int64_t(j->yStride());
    const int64_t _stride_j_2 = int64_t(j->zStride());
    const int64_t _stride_j_3 = int64_t(1 * int64_t(j->fStride()));
    internal_unpack_N::unpack_N(_data_buffer, _data_j, _size_j_0, _size_j_1,
                                _size_j_2, _stride_j_0, _stride_j_1,
                                _stride_j_2, _stride_j_3);
    break;
  }

  case stencil::TN: {
    float *RESTRICT const _data_buffer = buffer;
    WALBERLA_ASSERT_GREATER_EQUAL(ci.xMin(), -int_c(j->nrOfGhostLayers()));
    WALBERLA_ASSERT_GREATER_EQUAL(ci.yMin(), -int_c(j->nrOfGhostLayers()));
    WALBERLA_ASSERT_GREATER_EQUAL(ci.zMin(), -int_c(j->nrOfGhostLayers()));
    float *RESTRICT _data_j = j->dataAt(ci.xMin(), ci.yMin(), ci.zMin(), 0);
    WALBERLA_ASSERT_GREATER_EQUAL(j->xSizeWithGhostLayer(),
                                  int64_t(cell_idx_c(ci.xSize()) + 0));
    const int64_t _size_j_0 = int64_t(cell_idx_c(ci.xSize()) + 0);
    WALBERLA_ASSERT_GREATER_EQUAL(j->ySizeWithGhostLayer(),
                                  int64_t(cell_idx_c(ci.ySize()) + 0));
    const int64_t _size_j_1 = int64_t(cell_idx_c(ci.ySize()) + 0);
    WALBERLA_ASSERT_GREATER_EQUAL(j->zSizeWithGhostLayer(),
                                  int64_t(cell_idx_c(ci.zSize()) + 0));
    const int64_t _size_j_2 = int64_t(cell_idx_c(ci.zSize()) + 0);
    const int64_t _stride_j_0 = int64_t(j->xStride());
    const int64_t _stride_j_1 = int64_t(j->yStride());
    const int64_t _stride_j_2 = int64_t(j->zStride());
    const int64_t _stride_j_3 = int64_t(1 * int64_t(j->fStride()));
    internal_unpack_TN::unpack_TN(_data_buffer, _data_j, _size_j_0, _size_j_1,
                                  _size_j_2, _stride_j_0, _stride_j_1,
                                  _stride_j_2, _stride_j_3);
    break;
  }

  default:
    WALBERLA_ASSERT(false);
  }
}

uint_t DensityPackInfo_single_precision::size(stencil::Direction dir,
                                              const IBlock *block) const {
  auto j = block->getData<field::GhostLayerField<float, 13>>(jID);

  CellInterval ci;
  j->getGhostRegion(dir, ci, 1, false);

  uint_t elementsPerCell = 0;

  switch (dir) {
  case stencil::BSW:
    elementsPerCell = 1;
    break;

  case stencil::SW:
    elementsPerCell = 3;
    break;

  case stencil::TSW:
    elementsPerCell = 1;
    break;

  case stencil::BW:
    elementsPerCell = 3;
    break;

  case stencil::W:
    elementsPerCell = 9;
    break;

  case stencil::TW:
    elementsPerCell = 3;
    break;

  case stencil::BNW:
    elementsPerCell = 1;
    break;

  case stencil::NW:
    elementsPerCell = 3;
    break;

  case stencil::TNW:
    elementsPerCell = 1;
    break;

  case stencil::BS:
    elementsPerCell = 2;
    break;

  case stencil::S:
    elementsPerCell = 6;
    break;

  case stencil::TS:
    elementsPerCell = 2;
    break;

  case stencil::B:
    elementsPerCell = 5;
    break;

  case stencil::T:
    elementsPerCell = 4;
    break;

  case stencil::BN:
    elementsPerCell = 1;
    break;

  case stencil::N:
    elementsPerCell = 3;
    break;

  case stencil::TN:
    elementsPerCell = 1;
    break;

  default:
    elementsPerCell = 0;
  }
  return ci.numCells() * elementsPerCell * sizeof(float);
}

} // namespace pystencils
} // namespace walberla