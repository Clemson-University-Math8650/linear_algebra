
#ifndef __Math8650_LinearSolverInterface_HPP__
#define __Math8650_LinearSolverInterface_HPP__

#include <memory>

#include "DenseVector.hpp"

#include "DenseMatrix.hpp"
#include "SymmetricMatrix.hpp"
#include "TridiagonalMatrix.hpp"
#include "UpperTriangularMatrix.hpp"
#include "LowerTriangularMatrix.hpp"

namespace math8650
{

/**
 *
 */
class LinearSolverInterface
{
  
  static const int maxiter = 10000;
  constexpr static double tol = 1.0e-06;
  constexpr static double eps = 1.0e-15;

public:

  //! \brief method solves a symmetric linear system using the Jacobi preconditioned CG method
  static void solveSystemCG(const std::shared_ptr<const DenseMatrix>& A, 
                            const std::shared_ptr<const DenseVector>& rhs, 
                            const std::shared_ptr<DenseVector>& sol);

  //! \brief method solves a linear system using the Jacobi preconditioned BiCG method
  static void solveSystemBiCG(const std::shared_ptr<const DenseMatrix>& A, 
                              const std::shared_ptr<const DenseVector>& rhs, 
                              const std::shared_ptr<DenseVector>& sol);

  //! \brief method solves a symmetric linear system using the Jacobi preconditioned CG method
  static void solveSystemCG(const std::shared_ptr<const SymmetricMatrix>& A, 
                            const std::shared_ptr<const DenseVector>& rhs, 
                            const std::shared_ptr<DenseVector>& sol);

  //! \brief method solves a linear system using the Jacobi preconditioned BiCG method
  static void solveSystemBiCG(const std::shared_ptr<const SymmetricMatrix>& A, 
                              const std::shared_ptr<const DenseVector>& rhs, 
                              const std::shared_ptr<DenseVector>& sol);

  //! \brief method solves a symmetric linear system using the Jacobi preconditioned CG method
  static void solveSystemCG(const std::shared_ptr<const TridiagonalMatrix>& A, 
                            const std::shared_ptr<const DenseVector>& rhs, 
                            const std::shared_ptr<DenseVector>& sol);

  //! \brief method solves a linear system using the Jacobi preconditioned BiCG method
  static void solveSystemBiCG(const std::shared_ptr<const TridiagonalMatrix>& A, 
                              const std::shared_ptr<const DenseVector>& rhs, 
                              const std::shared_ptr<DenseVector>& sol);

  //! \brief method solves a symmetric linear system using the Jacobi preconditioned CG method
  static void solveSystemCG(const std::shared_ptr<const UpperTriangularMatrix>& A, 
                            const std::shared_ptr<const DenseVector>& rhs, 
                            const std::shared_ptr<DenseVector>& sol);

  //! \brief method solves a linear system using the Jacobi preconditioned BiCG method
  static void solveSystemBiCG(const std::shared_ptr<const UpperTriangularMatrix>& A, 
                              const std::shared_ptr<const DenseVector>& rhs, 
                              const std::shared_ptr<DenseVector>& sol);
  
  //! \brief method solves a symmetric linear system using the Jacobi preconditioned CG method
  static void solveSystemCG(const std::shared_ptr<const LowerTriangularMatrix>& A, 
                            const std::shared_ptr<const DenseVector>& rhs, 
                            const std::shared_ptr<DenseVector>& sol);

  //! \brief method solves a linear system using the Jacobi preconditioned BiCG method
  static void solveSystemBiCG(const std::shared_ptr<const LowerTriangularMatrix>& A, 
                              const std::shared_ptr<const DenseVector>& rhs, 
                              const std::shared_ptr<DenseVector>& sol);
  
};

}
// end namespace math8650

#endif /** __Math8650_LinearSolverInterface_HPP__ */
