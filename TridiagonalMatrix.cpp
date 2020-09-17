
#include "TridiagonalMatrix.hpp"
#include "DenseVector.hpp"

namespace math8650
{

void TridiagonalMatrix::allocate(std::size_t r, std::size_t s)
{
  if (!((r > 0) && (s > 0) && (r==s)))
  {
    m_is_allocated = false;
  } else {
  
    try {
      
      m_diag = new double[r];
      m_above = new double[r-1];
      m_below = new double[r-1];
     
    } catch(std::bad_alloc& bad) {
      
      std::cerr << "bad alloc caught: " << bad.what() << std::endl;
      throw;
    }
    m_is_allocated = true;
  }
}
/// end void allocate

void TridiagonalMatrix::deallocate()
{
  if (!m_is_allocated)
    return;
  delete [] m_diag; delete [] m_above; delete [] m_below; 
  m_rows = 0; m_cols = 0;
  m_is_allocated = false;
}
/// end void deallocate

TridiagonalMatrix::TridiagonalMatrix(std::size_t p, std::size_t q, const double val) : Matrix(p,q)
{
  allocate(p,q);
  for (std::size_t i = 0; i < p; ++i)
    m_diag[i] = val;
  for (std::size_t i = 0; i < p-1; ++i)
  {
    m_above[i] = val;
    m_below[i] = val;
  }

}
/// end TridiagonalMatrix

TridiagonalMatrix::~TridiagonalMatrix()
{
  deallocate();
}
/// end Destructor ~TridiagonalMatrix

double TridiagonalMatrix::operator()(std::size_t i, std::size_t j) const
{
  if (i==j)
    return m_diag[i];
  if (i == j+1)
    return m_below[i]; 
  if (j == i+1)
    return m_above[j];
  else
    return 0.0;
}
/// end double operator()

double& TridiagonalMatrix::operator()(std::size_t i, std::size_t j)
{
  if (i==j)
    return m_diag[i];
  if (i == j+1)
    return m_below[i]; 
  if (j == i+1)
    return m_above[j];
  else
    throw std::runtime_error(" ... ");
}
/// end double& operator()

DenseVector TridiagonalMatrix::operator*(const DenseVector& vec) const
{
  assert(m_is_allocated);
  if (m_cols != vec.size())
    throw std::runtime_error("Error in \"TridiagonalMatrix::operator*\" "
                             "addition not possible, matrix sizes do not agree");
  
  DenseVector temp_vector(m_rows,0.0);
  for (std::size_t i = 0; i < m_rows; i++)
    for (std::size_t j = 0; j < m_cols; j++) temp_vector[i] += (*this)(i,j)*vec[j];
  return temp_vector;
}
/// end Vector operator*

DenseVector TridiagonalMatrix::trans_mult(const DenseVector& vec) const
{
  assert(m_is_allocated);
  if (m_cols != vec.size())
    throw std::runtime_error("Error in \"TridiagonalMatrix::operator*\" "
                             "addition not possible, matrix sizes do not agree");
  
  DenseVector temp_vector(m_rows,0.0);
  for (std::size_t i = 0; i < m_rows; i++)
    for (std::size_t j = 0; j < m_cols; j++) temp_vector[i] += (*this)(j,i)*vec[j];
  return temp_vector;
}
/// end Vector operator*

void TridiagonalMatrix::getDiagonal(DenseVector& vec) const
{
  vec.allocate(m_rows);
  for (std::size_t i = 0; i < m_rows; i++)
    vec(i) = m_diag[i];
}
/// end void Transpose

}
// end namespace math8650
