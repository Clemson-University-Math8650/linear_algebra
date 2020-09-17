
#include <vector>
#include <functional>

#include "JacobiPreconditioner.hpp"
#include "LinearSolverInterface.hpp"

using namespace math8650;


// example 3: ODE BVP
void problem3();

enum class MatrixType { DenseMatrix, SymmetricMatrix, TridiagonalMatrix, UpperTriangularMatrix };

class Problem_ODE_BVP
{
public:

  //! \brief Ctor
  Problem_ODE_BVP(const double a,
                  const double b,
                  const int n,
                  const double u_a,
                  const double u_b,
                  const double alpha,
                  const std::function<double(const double)>& ode_f,
                  const MatrixType mat_type)
    : m_a(a), m_b(b), m_n(n), 
      m_u_a(u_a), m_u_b(u_b), 
      m_alpha(alpha), m_ode_f(ode_f), m_mat_type(mat_type) {}

  //! \brief Preallocate system objs
  void initialize();
  
  //! \brief Build and solve problem
  std::vector<std::shared_ptr<DenseVector>>
  run() const;
                  
private:

  double m_a;
  double m_b;
  int m_n;

  double m_u_a;
  double m_u_b;

  double m_alpha;

  std::function<double(const double)> m_ode_f;
  
  MatrixType m_mat_type;
  
  std::shared_ptr<DenseMatrix> A_d;
  std::shared_ptr<SymmetricMatrix> A_s;
  std::shared_ptr<TridiagonalMatrix> A_t;
  
  std::shared_ptr<DenseVector> m_rhs, m_U;
  
};

int main()
{
  const double a = 0.0;
  const double b = 1.0;
  const double T_a = 0.0;
  const double T_b = 2.0;
  const double alpha = 4.0;
  const auto f = [](const double x)->double
  { return -4.0*x; };
  
  // system sizes
  std::vector<int> ns = {10, 20, 40, 80, 160};
  // pointers to the problems
  std::vector<std::unique_ptr<Problem_ODE_BVP>> problems;
  // build the problems
  for (const auto it : ns)
    problems.push_back(std::make_unique<Problem_ODE_BVP>(a,b,it,T_a,T_b,alpha,f,MatrixType::DenseMatrix));
  // initialize, run, and destroy the problems
  for (auto& it : problems)
  {
    it->initialize();
    it->run();
    it.reset();
  }
  
  return 0;
}


void Problem_ODE_BVP::initialize()
{
  switch (m_mat_type) {
    
    case MatrixType::DenseMatrix :
      A_d = std::make_shared<DenseMatrix>(m_n-1,m_n-1);
      break;
    case MatrixType::SymmetricMatrix :
      A_s = std::make_shared<SymmetricMatrix>(m_n-1,m_n-1);
      break;
    case MatrixType::TridiagonalMatrix :
      A_t = std::make_shared<TridiagonalMatrix>(m_n-1,m_n-1);
      break;
    default :
      throw std::runtime_error("Error in \"Problem_ODE_BVP::initialize()\" no valid matrix type chosen");
  }
  
  m_rhs = std::make_shared<DenseVector>(m_n-1);
  m_U = std::make_shared<DenseVector>(m_n-1);
}

std::vector<std::shared_ptr<DenseVector>> 
Problem_ODE_BVP::run() const
{
  std::shared_ptr<DenseVector> x = std::make_shared<DenseVector>(m_n+1);
  const double h = (m_b - m_a)/m_n;
  for (std::size_t i (0); i < x->size(); ++i)
    (*x)[i] = m_a + i*h;

  // build system matrix
  switch (m_mat_type) {
    case MatrixType::DenseMatrix :
      {
        for (std::size_t i (0); i < A_d->numRows(); ++i)
          (*A_d)(i,i) = - (2.0 + m_alpha * h*h);
        
        for (std::size_t i (1); i < A_d->numRows(); ++i)
        {
          (*A_d)(i,i-1) = 1.0;
          (*A_d)(i-1,i) = 1.0;
        }
      }
      break;
    case MatrixType::SymmetricMatrix :
      {
        for (std::size_t i (0); i < A_s->numRows(); ++i)
          (*A_s)(i,i) = - (2.0 + m_alpha * h*h);
        for (std::size_t i (1); i < A_s->numRows(); ++i)
          (*A_d)(i,i-1) = 1.0;
      }
      break;
    case MatrixType::TridiagonalMatrix :
      {
        for (std::size_t i (0); i < A_d->numRows(); ++i)
          (*A_d)(i,i) = - (2.0 + m_alpha * h*h);
        
        for (std::size_t i (1); i < A_d->numRows(); ++i)
        {
          (*A_d)(i,i-1) = 1.0;
          (*A_d)(i-1,i) = 1.0;
        }
      }
      break;
    default :
      throw std::runtime_error("Error in \"Problem_ODE_BVP::initialize()\" no valid matrix type chosen");
  }
  
  // build RHS
  for (int i (0); i < m_n-1; ++i)
  {
    (*m_rhs)(i) = h*h*m_ode_f((*x)[i+1]);
  }
  (*m_rhs)(0) -= m_u_a;
  (*m_rhs)(m_n-2) -= m_u_b;

  switch (m_mat_type) {
    
    case MatrixType::DenseMatrix :
      LinearSolverInterface::solveSystemBiCG(A_d,m_rhs,m_U);
      break;
    case MatrixType::SymmetricMatrix :
      LinearSolverInterface::solveSystemBiCG(A_s,m_rhs,m_U);
      break;
    case MatrixType::TridiagonalMatrix :
      LinearSolverInterface::solveSystemBiCG(A_t,m_rhs,m_U);
      break;
    default :
      throw std::runtime_error("Error in \"Problem_ODE_BVP::initialize()\" no valid matrix type chosen");
  }
    
  return {m_U,x};  
}



void problem3()
{
  const double a = 0.0;
  const double b = 1.0;
  const int n = 1000;

  const double T_a = 0.0;
  const double T_b = 2.0;

  const double alpha = 4.0;

  const auto f = [] (const double x)->double { return -4.0*x; };
  const auto T_exact = [](const double x)->double 
  { return exp(2.0)* (exp(2.0*x)- exp(-2.0*x))/(exp(4.0)-1.0) + x;  };

  std::shared_ptr<SymmetricMatrix> A = std::make_shared<SymmetricMatrix>(n-1,n-1);
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
//     (*A)(i-1,i) = 1.0;
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

