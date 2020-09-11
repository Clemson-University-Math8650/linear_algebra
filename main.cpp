
#include <vector>

#include "JacobiPreconditioner.hpp"
#include "LinearSolverInterface.hpp"

using namespace math8650;

// example 1: symmetric matrix 
void problem1();

// example 2: non-symmetric matrix
void problem2();

// example 3: ODE BVP
void problem3();



int main()
{
  
  problem1();
  problem2();
  problem3();

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


void problem3()
{
  const double a = 0.0;
  const double b = 1.0;
  const int n = 500;

  const double T_a = 0.0;
  const double T_b = 2.0;

  const double alpha = 4.0;

  const auto f = [] (const double x)->double { return -4.0*x; };
  const auto T_exact = [](const double x)->double 
  { return exp(2.0)* (exp(2.0*x)- exp(-2.0*x))/(exp(4.0)-1.0) + x;  };

  auto A = std::make_shared<DenseMatrix>(n-1,n-1);
  auto rhs = std::make_shared<DenseVector>(n-1);
  auto T = std::make_shared<DenseVector>(n-1);
  auto T_vec_exact = std::make_shared<DenseVector>(n-1);

  const double h = (b - a)/n;
  std::vector<double> x(n+1);
  for (std::size_t i (0); i < x.size(); ++i)
    x[i] = a + i*h;

  for (std::size_t i (0); i < A->numRows(); ++i)
  {
    (*A)(i,i) = - (2.0 + alpha * h*h);
    (*rhs)(i) = h*h*f(x[i+1]);
    (*T_vec_exact)(i) = T_exact(x[i+1]);
  }
  (*rhs)(0) -= T_a;
  (*rhs)(n-2) -= T_b;

  for (std::size_t i (1); i < A->numRows(); ++i)
  {
    (*A)(i,i-1) = 1.0;
    (*A)(i-1,i) = 1.0;
  }
  
  std::cout << (*A) << std::endl;
  std::cout << (*rhs) << std::endl;
  
  LinearSolverInterface::solveSystemBiCG(A,rhs,T);
  std::cout << *T << std::endl;
  
  std::cout << std::endl;
  std::cout << std::endl;
  std::cout << std::endl;
  
  // compare exact and numerical solution
  for (std::size_t i (0); i < A->numRows(); ++i)
  {
    std::cout << x[i+1] << " " << (*T)(i) << " " << (*T_vec_exact)(i) << std::endl;
  }
  
}

