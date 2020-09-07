
#include "JacobiPreconditioner.hpp"
#include "LinearSolverInterface.hpp"

using namespace math8650;

void problem1();
void problem2();
void problem3();



int main()
{
  
  problem1();
  problem2();
  

  return 0;
}

void problem1()
{
  const std::size_t n = 10;
  auto A = std::make_shared<DenseMatrix>(n,n);
  for (std::size_t i = 0; i < A->numRows(); ++i)
    (*A)(i,i) = 1.0;
  for (std::size_t i = 1; i < A->numRows(); ++i)
    (*A)(i,i-1) = -0.5;
  for (std::size_t j = 1; j < A->numCols(); ++j)
    (*A)(j-1,j) = -0.5;
  
  auto p = std::make_shared<DenseVector>(n);
  auto rhs = std::make_shared<DenseVector>(n);
  (*rhs)(0) = 0.5;
  
  LinearSolverInterface::solveSystemCG(A,rhs,p);
  std::cout << *p << std::endl;
}

void problem2()
{

  auto A = std::make_shared<DenseMatrix>(3,3);
  (*A)(0,0) = 1.0;
  (*A)(0,2) = -1.0;
  (*A)(1,0) = -0.5;
  (*A)(1,1) = 1.0;
  (*A)(1,2) = -0.25;
  (*A)(2,0) = 1.0;
  (*A)(2,1) = -0.5;
  (*A)(2,2) = 1.0;
  
  auto rhs = std::make_shared<DenseVector>(A->numRows());
  auto sol = std::make_shared<DenseVector>(A->numCols());
  
  (*rhs)(0) = 0.2;
  (*rhs)(1) = -1.425;
  (*rhs)(2) = 2.0;
  
  std::cout << *A << std::endl;
  std::cout << *rhs << std::endl;

  LinearSolverInterface::solveSystemBiCG(A,rhs,sol);
  std::cout << *sol << std::endl;
  
}

