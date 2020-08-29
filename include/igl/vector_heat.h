#ifndef IGL_VECTOR_HEAT_H
#define IGL_VECTOR_HEAT_H

#include "igl_inline.h"
#include "min_quad_with_fixed.h"
#include <Eigen/Core>
#include <vector>

namespace igl {
  // "Vector heat method"
  //
  // Inputs:
  //   V  #V by dim list of triangle mesh vertex position (e.g., dim = 3 for 3D)
  //   F  #F by 3 list of triangle indices into rows of V
  //   b  #b list of indices into V of vertices where input vectors Xb are
  //     defined
  //   Xb  #b by dim list of (extrinsic) input vectors
  //   t  timestep parameter
  // Outputs:
  //   X  #V by dim list of per-vertex output (extrinsic) vectors 
  template <
    typename DerivedV,
    typename DerivedF,
	typename Derivedb,
  	typename DerivedXb,
    typename DerivedX>
  IGL_INLINE void vector_heat_method(
    const Eigen::MatrixBase<DerivedV> & V,
    const Eigen::MatrixBase<DerivedF> & F,
	const Eigen::MatrixBase<Derivedb> & b,
    const Eigen::MatrixBase<DerivedXb> & Xb,
    const typename DerivedV::Scalar t,
    Eigen::PlainObjectBase<DerivedX> & X);
}

#ifndef IGL_STATIC_LIBRARY
#  include "vector_heat.cpp"
#endif

#endif
