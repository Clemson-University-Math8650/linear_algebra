
#include "JacobiPreconditioner.hpp"
#include "LinearSolverInterface.hpp"

using namespace math8650;

void problem1();
void problem2();
void problem3();



int main()
{
  
  problem1();
  

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
  
  std::cout << *A << std::endl;
  std::cout << *rhs << std::endl;
  
  LinearSolverInterface::solveSystemCG(A,rhs,p);
  
  std::cout << *p << std::endl;
  
}
