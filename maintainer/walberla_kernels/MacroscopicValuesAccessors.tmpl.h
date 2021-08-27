//======================================================================================================================
//
//  This file is part of waLBerla. waLBerla is free software: you can
//  redistribute it and/or modify it under the terms of the GNU General Public
//  License as published by the Free Software Foundation, either version 3 of
//  the License, or (at your option) any later version.
//
//  waLBerla is distributed in the hope that it will be useful, but WITHOUT
//  ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
//  FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
//  for more details.
//
//  You should have received a copy of the GNU General Public License along
//  with waLBerla (see COPYING.txt). If not, see <http://www.gnu.org/licenses/>.
//
//! \\author Martin Bauer <martin.bauer@fau.de>
//
//======================================================================================================================

#pragma once
#include "core/DataTypes.h"
#include "core/logging/Logging.h"

#include "field/GhostLayerField.h"
#include "field/SwapableCompare.h"
#include "domain_decomposition/BlockDataID.h"
#include "domain_decomposition/IBlock.h"
#include "stencil/{{stencil_name}}.h"

#include "lbm/lattice_model/EquilibriumDistribution.h"
#include "lbm/field/Density.h"
#include "lbm/field/DensityAndMomentumDensity.h"
#include "lbm/field/DensityAndVelocity.h"
#include "lbm/field/PressureTensor.h"
#include "lbm/field/ShearRate.h"

#include <vector>

#ifdef __GNUC__
#define RESTRICT __restrict__
#elif _MSC_VER
#define RESTRICT __restrict
#else
#define RESTRICT
#endif

#ifdef WALBERLA_CXX_COMPILER_IS_GNU
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wunused-variable"
#pragma GCC diagnostic ignored "-Wunused-parameter"
#endif

#ifdef WALBERLA_CXX_COMPILER_IS_CLANG
#pragma clang diagnostic push
#pragma clang diagnostic ignored "-Wunused-variable"
#pragma clang diagnostic ignored "-Wunused-parameter"
#endif

{% set lmIgnores = ('pdfs', 'pdfs_tmp') %}
{% set lmOffsets = ('block_offset_0', 'block_offset_1', 'block_offset_2') %}



namespace walberla {
namespace {{namespace}} {


//======================================================================================================================
//
//  Implementation of macroscopic value backend
//
//======================================================================================================================

template<>
class EquilibriumDistribution< {{class_name}}, void>
{
public:
   typedef stencil::D3Q19 Stencil;

   static real_t get( const stencil::Direction direction,
                      const Vector3< real_t > & u = Vector3< real_t >( real_t(0.0) ),
                      real_t rho = real_t(1.0) )
   {
        {% if not compressible %}
        rho -= real_t(1.0);
        {% endif %}
        {{equilibrium_from_direction}}
   }

   static real_t getSymmetricPart( const stencil::Direction direction,
                                   const Vector3<real_t> & u = Vector3< real_t >(real_t(0.0)),
                                   real_t rho = real_t(1.0) )
   {
        {% if not compressible %}
        rho -= real_t(1.0);
        {% endif %}
        {{symmetric_equilibrium_from_direction}}
   }

   static real_t getAsymmetricPart( const stencil::Direction direction,
                                    const Vector3< real_t > & u = Vector3<real_t>( real_t(0.0) ),
                                    real_t rho = real_t(1.0) )
   {
        {% if not compressible %}
        rho -= real_t(1.0);
        {% endif %}
        {{asymmetric_equilibrium_from_direction}}
   }

   static std::vector< real_t > get( const Vector3< real_t > & u = Vector3<real_t>( real_t(0.0) ),
                                     real_t rho = real_t(1.0) )
   {
      {% if not compressible %}
      rho -= real_t(1.0);
      {% endif %}

      std::vector< real_t > equilibrium( Stencil::Size );
      for( auto d = Stencil::begin(); d != Stencil::end(); ++d )
      {
         equilibrium[d.toIdx()] = get(*d, u, rho);
      }
      return equilibrium;
   }
};


namespace internal {

template<>
struct AdaptVelocityToForce<{{class_name}}, void>
{
   template< typename FieldPtrOrIterator >
   static Vector3<real_t> get( FieldPtrOrIterator & it, const {{class_name}} & lm,
                               const Vector3< real_t > & velocity, const real_t rho )
   {
      auto x = it.x();
      auto y = it.y();
      auto z = it.z();
      {% if macroscopic_velocity_shift %}
      return velocity - Vector3<real_t>({{macroscopic_velocity_shift | join(",") }} {% if D == 2 %}, real_t(0.0) {%endif %} );
      {% else %}
      return velocity;
      {% endif %}
   }

   static Vector3<real_t> get( const cell_idx_t x, const cell_idx_t y, const cell_idx_t z, const {{class_name}} & lm,
                               const Vector3< real_t > & velocity, const real_t rho )
   {
      {% if macroscopic_velocity_shift %}

      return velocity - Vector3<real_t>({{macroscopic_velocity_shift | join(",") }} {% if D == 2 %}, real_t(0.0) {%endif %} );
      {% else %}
      return velocity;
      {% endif %}
   }

   static Vector3<real_t> get( const cell_idx_t x, const cell_idx_t y, const cell_idx_t z,
                               const GhostLayerField<real_t, 3u> * force_field,
                               const Vector3< real_t > & velocity, const real_t rho ) {

      return get(x, y, z, *force_field, velocity, rho);
   }
};
} // namespace internal



template<>
struct Equilibrium< {{class_name}}, void >
{

   template< typename FieldPtrOrIterator >
   static void set( FieldPtrOrIterator & it,
                    const Vector3< real_t > & u = Vector3< real_t >( real_t(0.0) ), real_t rho = real_t(1.0) )
   {
        {%if not compressible %}
        rho -= real_t(1.0);
        {%endif %}

       {% for eqTerm in equilibrium -%}
       it[{{loop.index0 }}] = {{eqTerm}};
       {% endfor -%}
   }

   template< typename PdfField_T >
   static void set( PdfField_T & pdf, const cell_idx_t x, const cell_idx_t y, const cell_idx_t z,
                    const Vector3< real_t > & u = Vector3< real_t >( real_t(0.0) ), real_t rho = real_t(1.0) )
   {
      {%if not compressible %}
      rho -= real_t(1.0);
      {%endif %}

      real_t & xyz0 = pdf(x,y,z,0);
      {% for eqTerm in equilibrium -%}
      pdf.getF( &xyz0, {{loop.index0 }})= {{eqTerm}};
      {% endfor -%}
   }
};


template<>
struct Density<{{class_name}}, void>
{
   template< typename FieldPtrOrIterator >
   static inline real_t get( {{class_name}} const & , const FieldPtrOrIterator & it )
   {
        {% for i in range(Q) -%}
            const real_t f_{{i}} = it[{{i}}];
        {% endfor -%}
        {{density_getters | indent(8)}}
        return rho;
   }

   template< typename PdfField_T >
   static inline real_t get( {{class_name}} const & ,
                             const PdfField_T & pdf, const cell_idx_t x, const cell_idx_t y, const cell_idx_t z )
   {
        const real_t & xyz0 = pdf(x,y,z,0);
        {% for i in range(Q) -%}
            const real_t f_{{i}} = pdf.getF( &xyz0, {{i}});
        {% endfor -%}
        {{density_getters | indent(8)}}
        return rho;
   }
};


template<>
struct DensityAndVelocity<{{class_name}}>
{
    template< typename FieldPtrOrIterator >
    static void set( FieldPtrOrIterator & it, const {{class_name}} & lm,
                     const Vector3< real_t > & u = Vector3< real_t >( real_t(0.0) ), const real_t rho_in = real_t(1.0) )
    {
        auto x = it.x();
        auto y = it.y();
        auto z = it.z();

        {{density_velocity_setter_macroscopic_values | indent(8)}}
        {% if D == 2 -%}
        const real_t u_2(0.0);
        {% endif %}

        Equilibrium<{{class_name}}>::set(it, Vector3<real_t>(u_0, u_1, u_2), rho{%if not compressible %} + real_t(1) {%endif%});
    }

    template< typename PdfField_T >
    static void set( PdfField_T & pdf, const cell_idx_t x, const cell_idx_t y, const cell_idx_t z, const {{class_name}} & lm,
                     const Vector3< real_t > & u = Vector3< real_t >( real_t(0.0) ), const real_t rho_in = real_t(1.0) )
    {
        {{density_velocity_setter_macroscopic_values | indent(8)}}
        {% if D == 2 -%}
        const real_t u_2(0.0);
        {% endif %}

        Equilibrium<{{class_name}}>::set(pdf, x, y, z, Vector3<real_t>(u_0, u_1, u_2), rho {%if not compressible %} + real_t(1) {%endif%});
    }
};


template<typename FieldIteratorXYZ >
struct DensityAndVelocityRange<{{class_name}}, FieldIteratorXYZ>
{

   static void set( FieldIteratorXYZ & begin, const FieldIteratorXYZ & end, const {{class_name}} & lm,
                    const Vector3< real_t > & u = Vector3< real_t >( real_t(0.0) ), const real_t rho_in = real_t(1.0) )
   {
        for( auto cellIt = begin; cellIt != end; ++cellIt )
        {
            const auto x = cellIt.x();
            const auto y = cellIt.y();
            const auto z = cellIt.z();
            {{density_velocity_setter_macroscopic_values | indent(12)}}
            {% if D == 2 -%}
            const real_t u_2(0.0);
            {% endif %}

            Equilibrium<{{class_name}}>::set(cellIt, Vector3<real_t>(u_0, u_1, u_2), rho{%if not compressible %} + real_t(1) {%endif%});
        }
   }
};



template<>
struct DensityAndMomentumDensity<{{class_name}}>
{
   template< typename FieldPtrOrIterator >
   static real_t get( Vector3< real_t > & momentumDensity, const {{class_name}} & lm,
                      const FieldPtrOrIterator & it )
   {
        const auto x = it.x();
        const auto y = it.y();
        const auto z = it.z();

        {% for i in range(Q) -%}
            const real_t f_{{i}} = it[{{i}}];
        {% endfor -%}

        {{momentum_density_getter | indent(8) }}
        {% for i in range(D) -%}
            momentumDensity[{{i}}] = md_{{i}};
        {% endfor %}
        return rho;
   }

   template< typename PdfField_T >
   static real_t get( Vector3< real_t > & momentumDensity, const {{class_name}} & lm, const PdfField_T & pdf,
                      const cell_idx_t x, const cell_idx_t y, const cell_idx_t z )
   {
        const real_t & xyz0 = pdf(x,y,z,0);
        {% for i in range(Q) -%}
            const real_t f_{{i}} = pdf.getF( &xyz0, {{i}});
        {% endfor -%}

        {{momentum_density_getter | indent(8) }}
        {% for i in range(D) -%}
            momentumDensity[{{i}}] = md_{{i}};
        {% endfor %}
       return rho;
   }
};


template<>
struct MomentumDensity< {{class_name}}>
{
   template< typename FieldPtrOrIterator >
   static void get( Vector3< real_t > & momentumDensity, const {{class_name}} & lm, const FieldPtrOrIterator & it )
   {
        const auto x = it.x();
        const auto y = it.y();
        const auto z = it.z();

        {% for i in range(Q) -%}
            const real_t f_{{i}} = it[{{i}}];
        {% endfor -%}

        {{momentum_density_getter | indent(8) }}
        {% for i in range(D) -%}
            momentumDensity[{{i}}] = md_{{i}};
        {% endfor %}
   }

   template< typename PdfField_T >
   static void get( Vector3< real_t > & momentumDensity, const {{class_name}} & lm, const PdfField_T & pdf,
                    const cell_idx_t x, const cell_idx_t y, const cell_idx_t z )
   {
        const real_t & xyz0 = pdf(x,y,z,0);
        {% for i in range(Q) -%}
            const real_t f_{{i}} = pdf.getF( &xyz0, {{i}});
        {% endfor -%}

        {{momentum_density_getter | indent(8) }}
        {% for i in range(D) -%}
            momentumDensity[{{i}}] = md_{{i}};
        {% endfor %}
   }
};


template<>
struct PressureTensor<{{class_name}}>
{
   template< typename FieldPtrOrIterator >
   static void get( Matrix3< real_t > & /* pressureTensor */, {{class_name}} const & /* latticeModel */, const FieldPtrOrIterator & /* it */ )
   {
       WALBERLA_ABORT("Not implemented");
   }

   template< typename PdfField_T >
   static void get( Matrix3< real_t > & /* pressureTensor */, {{class_name}} const & /* latticeModel */, const PdfField_T & /* pdf */,
                    const cell_idx_t /* x */, const cell_idx_t /* y */, const cell_idx_t /* z */ )
   {
       WALBERLA_ABORT("Not implemented");
   }
};


template<>
struct ShearRate<{{class_name}}>
{
   template< typename FieldPtrOrIterator >
   static inline real_t get( {{class_name}} const & /* latticeModel */, const FieldPtrOrIterator & /* it */,
                             const Vector3< real_t > & /* velocity */, const real_t /* rho */)
   {
       WALBERLA_ABORT("Not implemented");
       return real_t(0.0);
   }

   template< typename PdfField_T >
   static inline real_t get( {{class_name}} const & latticeModel,
                             const PdfField_T & /* pdf */, const cell_idx_t /* x */, const cell_idx_t /* y */, const cell_idx_t /* z */,
                             const Vector3< real_t > & /* velocity */, const real_t /* rho */ )
   {
       WALBERLA_ABORT("Not implemented");
       return real_t(0.0);
   }

   static inline real_t get( const std::vector< real_t > & /* nonEquilibrium */, const real_t /* relaxationParam */,
                             const real_t /* rho */ = real_t(1) )
   {
       WALBERLA_ABORT("Not implemented");
       return real_t(0.0);
   }
};


} // namespace {{namespace}}
} // namespace walberla



#ifdef WALBERLA_CXX_COMPILER_IS_GNU
#pragma GCC diagnostic pop
#endif

#ifdef WALBERLA_CXX_COMPILER_IS_CLANG
#pragma clang diagnostic pop
#endif
