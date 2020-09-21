
#include "UpperTriangularMatrix.hpp"

namespace math8650
{

void UpperTriangularMatrix::allocate(std::size_t nrows, std::size_t ncols)
{
  if (!((nrows > 0) && (ncols > 0) && (nrows==ncols)))
  {
    m_is_allocated = false;
  } else {
    m_rows = nrows; m_cols = ncols;
    m_mat = new double*[m_rows];
    for (std::size_t i = 0; i < m_rows; ++i)
      m_mat[i] = new double[m_rows-i];
    m_is_allocated = true;
  }
}
/// end void allocate

void UpperTriangularMatrix::deallocate()
{
  if (!m_is_allocated)
    return;
  for (std::size_t i = 0; i < m_rows; i++) delete [] m_mat[i];
  delete [] m_mat; m_rows = 0; m_cols = 0;
  m_is_allocated = false;
}
/// end void deallocate

UpperTriangularMatrix::UpperTriangularMatrix(std::size_t p, std::size_t q, const double val)
{
  allocate(p,q);
  for (std::size_t i = 0; i < p; ++i)
    for (std::size_t j = 0; j < m_rows-i; ++j)
      m_mat[i][j] = val;
}
/// end UpperTriangularMatrix

UpperTriangularMatrix::~UpperTriangularMatrix()
{
  deallocate();
}
/// end Destructor ~UpperTriangularMatrix

double UpperTriangularMatrix::operator()(std::size_t i, std::size_t j) const
{
  if (j >= i)
    return m_mat[i][j-i];
  else
    return 0.0;
}
/// end double operator()

double& UpperTriangularMatrix::operator()(std::size_t i, std::size_t j)
{
  if (j >= i)
    return m_mat[i][j-i];
  else
    throw std::runtime_error("Error in \"UpperTriangularMatrix::operator()\" accessing unallocated lower triangular entries");
}
/// end double& operator()

DenseVector UpperTriangularMatrix::operator*(const DenseVector& vec) const
{
  assert(m_is_allocated);
  if (m_cols != vec.size())
    throw std::runtime_error("Error in \"UpperTriangularMatrix::operator*\" "
                             "addition not possible, matrix sizes do not agree");
  
  DenseVector temp_vector(m_rows,0.0);
  for (std::size_t i = 0; i < m_rows; i++)
    for (std::size_t j = i; j < m_cols; j++) temp_vector[i] += m_mat[i][j-i]*vec[j];
  return temp_vector;
}
/// end Vector operator*

DenseVector UpperTriangularMatrix::trans_mult(const DenseVector& vec) const
{
  assert(m_is_allocated);
  if (m_cols != vec.size())
    throw std::runtime_error("Error in \"UpperTriangularMatrix::trans_mult\" "
                             "addition not possible, matrix sizes do not agree");
  
  DenseVector temp_vector(m_rows,0.0);
  for (std::size_t i = 0; i < m_rows; i++)
    for (std::size_t j = 0; j < m_cols; j++)
      temp_vector[i] += (*this)(j,i)*vec[j];
  return temp_vector;
}
/// end Vector operator*

void UpperTriangularMatrix::getDiagonal(DenseVector& vec) const
{
  vec.allocate(m_rows);
  for (std::size_t i = 0; i < m_rows; i++)
    vec(i) = m_mat[i][0];
}
/// end void Transpose

}
// end namespace math8650
