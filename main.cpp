
#include <ctime>
#include <vector>
#include <functional>

#include "JacobiPreconditioner.hpp"
#include "LinearSolverInterface.hpp"

using namespace math8650;

class Chrono
{
  typedef std::clock_t clock_t;

  clock_t M_t1;
  clock_t M_t2;
  double M_diff;
  bool M_run;

public:

  Chrono();
  virtual ~Chrono() {}
  
  void reset();
  void start();
  void stop();
  double getTime();
};

enum class MatrixType { DenseMatrix, TridiagonalMatrix, SymmetricMatrix, UpperTriangularMatrix };

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
  
  mutable Chrono m_timer;
  
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

  // part 1: test run   
  //! ==========================================================================================================================
  {
    // system sizes
    std::vector<int> ns = {10, 20, 40, 80, 160, 320};
    // different matrix types
    std::vector<MatrixType> mat_types = {MatrixType::DenseMatrix, MatrixType::TridiagonalMatrix, MatrixType::SymmetricMatrix};
    // pointers to the problems
    std::vector<std::unique_ptr<Problem_ODE_BVP>> problems;
    // build the problems
    for (const auto it : ns)
      for (const auto mat_it : mat_types)
        problems.push_back(std::make_unique<Problem_ODE_BVP>(a,b,it,T_a,T_b,alpha,f,mat_it));
    // initialize, run, and destroy the problems
    for (auto& it : problems)
    {
      it->initialize();
      auto result = it->run();
      it.reset();
    }
  }
  //! ==========================================================================================================================
  
  // part 2: timing  
  //! ==========================================================================================================================
  {
    // system sizes
    std::vector<int> ns = {1000, 2000, 4000};
    // different matrix types
    std::vector<MatrixType> mat_types = {MatrixType::DenseMatrix, MatrixType::TridiagonalMatrix, MatrixType::SymmetricMatrix};
    // pointers to the problems
    std::vector<std::unique_ptr<Problem_ODE_BVP>> problems;
    // build the problems
    for (const auto it : ns)
      for (const auto mat_it : mat_types)
        problems.push_back(std::make_unique<Problem_ODE_BVP>(a,b,it,T_a,T_b,alpha,f,mat_it));
    // initialize, run, and destroy the problems
    for (auto& it : problems)
    {
      it->initialize();
      auto result = it->run();
      it.reset();
    }
  }
  //! ==========================================================================================================================
  
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
      throw std::runtime_error("Error in \"Problem_ODE_BVP::initialize()\" no valid matrix type chosen.");
  }
  
  m_rhs = std::make_shared<DenseVector>(m_n-1);
  m_U = std::make_shared<DenseVector>(m_n-1);
}

std::vector<std::shared_ptr<DenseVector>> 
Problem_ODE_BVP::run() const
{
  std::vector<std::string> mat_types_str = {"DenseMatrix", "TridiagonalMatrix", "SymmetricMatrix"};
  std::cout << "==================== Beginning of Test ========================================================" << std::endl;
  std::cout << "Running ODE BV problem with \"n           = " << m_n << "\"." << std::endl;
  std::cout << "Matrix type:                \"matrix type = " << mat_types_str.at(static_cast<int>(m_mat_type)) << "\"." << std::endl;
  std::cout << std::endl;
  std::cout << std::endl;
  std::cout << std::endl;
  
  m_timer.start();
  
  std::shared_ptr<DenseVector> x = std::make_shared<DenseVector>(m_n+1);
  const double h = (m_b - m_a)/m_n;
  for (std::size_t i (0); i < x->size(); ++i)
    (*x)[i] = m_a + i*h;

  // build system matrix
  switch (m_mat_type)
  {
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
          (*A_s)(i,i-1) = 1.0;
      }
      break;
    case MatrixType::TridiagonalMatrix :
      {
        for (std::size_t i (0); i < A_t->numRows(); ++i)
          (*A_t)(i,i) = - (2.0 + m_alpha * h*h);
        
        for (std::size_t i (1); i < A_t->numRows(); ++i)
        {
          (*A_t)(i,i-1) = 1.0;
          (*A_t)(i-1,i) = 1.0;
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
  
  m_timer.stop();
  std::cout << std::endl;
  std::cout << "Time for the simulation: " << m_timer.getTime() << " sec." << std::endl;
  std::cout << std::endl;
  std::cout << "==================== End of Test ==============================================================" << std::endl;
  std::cout << std::endl;
  std::cout << std::endl;
  std::cout << std::endl;
  return {m_U,x};  
}

Chrono::Chrono()
{
  M_t1 = clock_t(0);
  M_t2 = clock_t(0);
  M_diff = 0;
  M_run = false;
}
/// end constructor Chrono

double Chrono::getTime()
{
  return M_diff/CLOCKS_PER_SEC;
}
/// end double getTime

void Chrono::reset()
{
  M_t1 = 0;
  M_t2 = 0;
  M_diff = 0;
  M_run = false;
}
/// end void reset

void Chrono::start()
{
  M_t1 = std::clock();
  M_run = true;
}
/// end void start

void Chrono::stop()
{
  if ( M_run )
  {
    M_t2 = std::clock();
    M_diff = 1.0*(M_t2 - M_t1);
    M_run = false;
  } else {
    M_diff = 0.0;
    std::cout << "The clock was not started!!" << std::endl;
  }
}
/// end void stop

