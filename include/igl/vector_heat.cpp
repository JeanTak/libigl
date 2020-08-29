#include <vector>
#include <iostream>
#include "vector_heat.h"
#include "massmatrix.h"
#include "grad.h"
#include "invert_diag.h"
#include "cotmatrix.h"
#include "internal_angles.h"
#include "matlab_format.h"
#include "adjacency_list.h"
#include "PI.h"
#include "vertex_triangle_adjacency.h"

template <
  typename DerivedF,
  typename DerivedN>
float getAngle(
  const Eigen::MatrixBase<DerivedF> & F,
  const std::vector<std::vector<size_t> > VF,
  const Eigen::MatrixBase<DerivedN> & NormalizedAngle,
  const typename DerivedF::Scalar angle0,
  const typename DerivedF::Scalar angle1,
  const typename DerivedF::Scalar angle2
){

	int VFRows = VF[angle0].size();
	int angleIdx = 0;

	for(int i = 0; i < VFRows; i++){

    int faceIdx = VF[angle0][i];
		bool found = true;

		for(int j = 0; j < 3; j++){
			if(F( faceIdx, j) != angle0 &&
				F( faceIdx, j) != angle1 &&
				F( faceIdx, j) != angle2 ){
					found = false;
					break;
				}

			if(F( faceIdx, j) == angle0){
        angleIdx = j;
      }
		}

		if(found){
      return NormalizedAngle(faceIdx, angleIdx);
    }

	}
	return 0;
}

template <
  typename DerivedV,
  typename DerivedF,
  typename Derivedb,
  typename DerivedXb,
  typename DerivedX>
IGL_INLINE void igl::vector_heat_method(
  const Eigen::MatrixBase<DerivedV> & V,
  const Eigen::MatrixBase<DerivedF> & F,
  const Eigen::MatrixBase<Derivedb> & b,
  const Eigen::MatrixBase<DerivedXb> & Xb,
  const typename DerivedV::Scalar t,
  Eigen::PlainObjectBase<DerivedX> & X)
{
  using namespace Eigen;
  using namespace std;

  const int Frows = F.rows();
  const int VRows = V.rows();

  // Preparing Connection Laplacian L∇

  // Internal Angles in each mesh, sequence of angles corresponds to edges [1,2],[2,0],[0,1]
  Matrix<
    typename DerivedV::Scalar,
    DerivedF::RowsAtCompileTime,
    DerivedF::ColsAtCompileTime> InternalAngle;

  // total interior angle at vertex i
  VectorXd Theta_i = VectorXd::Zero(V.rows(), 1);
  igl::internal_angles(V,F,InternalAngle);

  // Normalize: each angle / total interior angle
  VectorXd NormalizedAngle(InternalAngle.rows(), InternalAngle.cols());

  for(int f = 0; f < Frows; f++)
  {
    for(int corner = 0; corner < 3; corner++)
    {
      Theta_i(F(f, corner), 0) += InternalAngle(f, corner);
    }
  }

  for(int f = 0; f < Frows; f++)
  {
    for(int corner = 0; corner < 3; corner++)
    {
      NormalizedAngle(f, corner) = InternalAngle(f, corner) * igl::PI / Theta_i(F(f, corner), 0);
    }
  }

  std::vector<std::vector<size_t> > VF, VFi;
  std::vector<std::vector<int>> AdjList;

  // vertex to face mapping
  igl::vertex_triangle_adjacency(V.rows(), F, VF, VFi);

  // vetex to adjacent vertexs mapping
  igl::adjacency_list(F, AdjList, true);

  // The edge rotations from i to jk in a circular cone - à la Knöppelet al. [2013]
  MatrixXd Phi = MatrixXd::Zero(VRows, VRows);

  for (int i = 0; i < VRows; i++)
  {
    for (int j = 1; j < AdjList[i].size(); j++)
    {
      Phi(i, AdjList[i][j]) = getAngle(F, VF, NormalizedAngle, i, AdjList[i][j-1], AdjList[i][j]) + Phi(i, AdjList[i][j-1]);
    }
  }


  // Mass Matrix M
  // Laplace-Beltrami operator L
  Eigen::SparseMatrix<double> M, L;

  igl::massmatrix(V,F,igl::MASSMATRIX_TYPE_DEFAULT,M);
  igl::cotmatrix(V,F,L);

  // (M + tL)u = u0
  Eigen::SparseMatrix<double> ML = M + t*L;

  //precompute and solver
  Eigen::SparseMatrix<double> _;
  igl::min_quad_with_fixed_data<double> Neumann;

  if(!igl::min_quad_with_fixed_precompute(ML, Eigen::VectorXi(), _, false, Neumann))
  {
    std::cout<<"min_quad_with_fixed_precompute FAILED!"<<std::endl;
  }

  // DerivedX delta0 = DerivedX::Zero(VRows, 1);
  // DerivedX u0 = DerivedX::Zero(VRows, 1);

  MatrixXd delta0 = DerivedX::Zero(VRows, 1);
  MatrixXd u0 = DerivedX::Zero(VRows, 1);
  // MatrixXd _Xb_(Xb.rows(), 1);
  MatrixXd u;

  for(int b_n = 0; b_n < b.rows(); b_n++)
  {
    delta0(b(b_n)) = 1;
    u0(b(b_n)) = fabs(Xb(b(b_n)));
  }

  if(!igl::min_quad_with_fixed_solve(Neumann, u0, DerivedX(), DerivedX(), u))
  {
    std::cout<<"min_quad_with_fixed_solve FAILED!"<<std::endl;
  }

  std::cout<<u<<std::endl;

  // std::cout<<igl::matlab_format(M,"M")<<std::endl;
	// std::cout<<igl::matlab_format(L,"L")<<std::endl;
	// std::cout<<igl::matlab_format(ML,"ML")<<std::endl;


}


#ifdef IGL_STATIC_LIBRARY
// Explicit template instantiation
template void
igl::vector_heat_method<Eigen::Matrix<double, -1, -1, 0, -1, -1>, Eigen::Matrix<int, -1, -1, 0, -1, -1>, Eigen::Matrix<int, -1, -1, 0, -1, -1>, Eigen::Matrix<double, -1, -1, 0, -1, -1>, Eigen::Matrix<double, -1, -1, 0, -1, -1> >(
  Eigen::MatrixBase<Eigen::Matrix<double, -1, -1, 0, -1, -1> > const &,
  Eigen::MatrixBase<Eigen::Matrix<int, -1, -1, 0, -1, -1> > const &,
  Eigen::MatrixBase<Eigen::Matrix<int, -1, -1, 0, -1, -1> > const &,
  Eigen::MatrixBase<Eigen::Matrix<double, -1, -1, 0, -1, -1> > const &,
  Eigen::Matrix<double, -1, -1, 0, -1, -1>::Scalar,
  Eigen::PlainObjectBase<Eigen::Matrix<double, -1, -1, 0, -1, -1> > &);
#endif