
#include "SymmetricMatrix.hpp"
#include "DenseVector.hpp"

namespace math8650
{

void SymmetricMatrix::allocate(std::size_t nrows, std::size_t ncols)
{
  if (!((nrows > 0) && (ncols > 0) && (nrows==ncols)))
  {
    m_is_allocated = false;
  } else {
    m_mat = new double*[nrows];
    // create a matrix only for the lower triangular part
    for (std::size_t i = 0; i < nrows; ++i)
      m_mat[i] = new double[i+1];
    m_is_allocated = true;
    m_rows = nrows; m_cols = ncols;
  }
}
/// end void allocate

void SymmetricMatrix::deallocate()
{
  if (!m_is_allocated)
    return;
  for (std::size_t i = 0; i < m_rows; i++) delete [] m_mat[i];
  delete [] m_mat; m_rows = 0; m_cols = 0;
  m_is_allocated = false;
}
/// end void deallocate

SymmetricMatrix::SymmetricMatrix(std::size_t nrows, std::size_t ncols, const double val)
{
  allocate(nrows,ncols);
  for (std::size_t i = 0; i < m_rows; ++i)
    for (std::size_t j = 0; j < i+1; ++j)
      m_mat[i][j] = val;
}
/// end SymmetricMatrix

SymmetricMatrix::~SymmetricMatrix()
{
  deallocate();
}
/// end Destructor ~SymmetricMatrix

double SymmetricMatrix::operator()(std::size_t i, std::size_t j) const
{
  assert((i >= 0) && (i < m_rows) && (j >= 0) && (j < m_cols));
  if (i >= j)
    return m_mat[i][j];
  else
    return m_mat[j][i];
}
/// end double operator()

double& SymmetricMatrix::operator()(std::size_t i, std::size_t j)
{
  assert((i >= 0) && (i < m_rows) && (j >= 0) && (j < m_cols));
  if (i >= j)
    return m_mat[i][j];
  else 
    throw std::runtime_error("Error in \"SymmetricMatrix::operator()\" "
                             "accessing unallocated upper triangular entries");
}
/// end double& operator()

DenseVector SymmetricMatrix::operator*(const DenseVector& vec) const
{
  assert(m_is_allocated);
  if (m_cols != vec.size())
    throw std::runtime_error("Error in \"SymmetricMatrix::operator*\" "
                             "addition not possible, matrix sizes do not agree");
  
  DenseVector temp_vector(m_rows,0.0);
  for (std::size_t i = 0; i < m_rows; i++)
    for (std::size_t j = 0; j < i+1; j++) 
    {
      temp_vector[i] += (*this)(i,j)*vec[j];
      if (i != j)
        temp_vector[j] += (*this)(j,i)*vec[i];
    }
  return temp_vector;
}
/// end Vector operator*

DenseVector SymmetricMatrix::trans_mult(const DenseVector& vec) const
{
  return (*this)*(vec);
}
/// end Vector operator*

void SymmetricMatrix::getDiagonal(DenseVector& vec) const
{
  vec.allocate(m_rows);
  for (std::size_t i = 0; i < m_rows; i++)
    vec(i) = m_mat[i][i];
}
/// end void Transpose

}
// end namespace math8650
