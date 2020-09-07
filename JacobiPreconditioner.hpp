
#ifndef __Math8650_JacobiPreconditioner_HPP__
#define __Math8650_JacobiPreconditioner_HPP__

namespace math8650
{

/**
 * Simple interface for the Jacobi preconditioner
 */
template <typename Matrix, typename Vector>
class JacobiPreconditioner
{

public:
  
  // Ctor: default to Jacobi preconditioner
  JacobiPreconditioner() {}
  
  // Ctor from matrix
  JacobiPreconditioner(const Matrix& A)
  {
    A.getDiagonal(m_prec);
  }
  
  // Dtor
  virtual ~JacobiPreconditioner() { m_prec.deallocate(); }

  // generates a preconditioner of type prec_type
  void Generate(const Matrix& A)
  {
    A.getDiagonal(m_prec);
  }
  
  // invert preconditioner
  Vector solve(const Vector& rhs) const
  {
    auto sol = rhs;
    for (int i = 0; i < int(m_prec.size()); ++i) 
      sol[i] = rhs[i]/m_prec[i];
    return sol;
    
  }
  
  // invert transpose preconditioner
  Vector trans_solve(const Vector& rhs) const
  {
    auto sol = rhs;
    for (int i = 0; i < int(m_prec.size()); ++i) 
      sol[i] = rhs[i]/m_prec[i];
    return sol;
  }  

private:
  
  // diagonal JacobiPreconditioner
  Vector m_prec;

};
/// end class JacobiPreconditioner

}
// end namespace math8650


#endif /** __Math8650_JacobiPreconditioner_HPP__ */
