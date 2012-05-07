#ifndef IGL_WRITEMESH_H
#define IGL_WRITEMESH_H
#include "igl_inline.h"

#include <string>
#include <vector>
#include <Eigen/Core>

namespace igl
{
  // save a tetrahedral volume mesh to a .mesh file
  //
  // Templates:
  //   Scalar  type for positions and vectors (will be cast as double)
  //   Index  type for indices (will be cast to int)
  // Input:
  //   mesh_file_name  path of .mesh file
  //   V  double matrix of vertex positions  #V by 3
  //   T  #T list of tet indices into vertex positions
  //   F  #F list of face indices into vertex positions
  template <typename Scalar, typename Index>
  IGL_INLINE bool writeMESH(
    const std::string mesh_file_name,
    const std::vector<std::vector<Scalar > > & V,
    const std::vector<std::vector<Index > > & T,
    const std::vector<std::vector<Index > > & F);

  // Templates:
  //   DerivedV  real-value: i.e. from MatrixXd
  //   DerivedT  integer-value: i.e. from MatrixXi
  //   DerivedF  integer-value: i.e. from MatrixXi
  // Input:
  //   mesh_file_name  path of .mesh file
  //   V  eigen double matrix #V by 3
  //   T  eigen int matrix #T by 4
  //   F  eigen int matrix #F by 3
  template <typename DerivedV, typename DerivedT, typename DerivedF>
  IGL_INLINE bool writeMESH(
    const std::string str,
    const Eigen::MatrixBase<DerivedV> & V, 
    const Eigen::MatrixBase<DerivedT> & T,
    const Eigen::MatrixBase<DerivedF> & F);
}

#ifdef IGL_HEADER_ONLY
#  include "writeMESH.cpp"
#endif

#endif