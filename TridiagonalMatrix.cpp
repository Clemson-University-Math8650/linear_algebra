
#include "TridiagonalMatrix.hpp"
#include "DenseVector.hpp"

namespace math8650
{

void TridiagonalMatrix::allocate(std::size_t nrows, std::size_t ncols)
{
  if (!((nrows > 0) && (ncols > 0) && (nrows==ncols)))
  {
    m_is_allocated = false;
  } else {
    m_rows = nrows; m_cols = ncols;
    m_diag = new double[m_rows];
    m_above = new double[m_rows-1];
    m_below = new double[m_rows-1];
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

TridiagonalMatrix::TridiagonalMatrix(std::size_t nrows, std::size_t ncols, const double val)
{
  allocate(nrows,ncols);
  for (std::size_t i = 0; i < m_rows; ++i)
    m_diag[i] = val;
  for (std::size_t i = 0; i < m_rows-1; ++i)
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
    throw std::runtime_error("Error in \"TridiagonalMatrix::operator()\" "
                             "accessing unallocated non-tridiagonal entries");
}
/// end double& operator()

DenseVector TridiagonalMatrix::operator*(const DenseVector& vec) const
{
  assert(m_is_allocated);
  if (m_cols != vec.size())
    throw std::runtime_error("Error in \"TridiagonalMatrix::operator*\" "
                             "addition not possible, matrix sizes do not agree");
  
  DenseVector temp_vector(m_rows,0.0);
  temp_vector[0] = m_diag[0] * vec[0] + m_above[0] * vec[1];
  for (std::size_t i = 1; i < m_rows-1; i++)
    temp_vector[i] = m_below[i] * vec[i-1] + m_diag[i] * vec[i] + m_above[i] * vec[i+1];
  temp_vector[m_rows-1] = m_above[m_rows-2] * vec[m_rows-2] + m_diag[m_rows-1] * vec[m_rows-1];
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
  temp_vector[0] = m_diag[0] * vec[0] + m_below[0] * vec[1];
  for (std::size_t i = 1; i < m_rows-1; i++)
    temp_vector[i] = m_above[i] * vec[i-1] + m_diag[i] * vec[i] + m_below[i] * vec[i+1];
  temp_vector[m_rows-1] = m_below[m_rows-2] * vec[m_rows-2] + m_diag[m_rows-1] * vec[m_rows-1];
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
