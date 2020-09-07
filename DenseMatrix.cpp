
#include "DenseMatrix.hpp"
#include "DenseVector.hpp"

namespace math8650
{

void DenseMatrix::allocate(std::size_t r, std::size_t s)
{
  if (!((r > 0) && (s > 0)))
  {
    m_is_allocated = false;
  } else {
  
    try {
      
      m_mat = new double*[r];
      for (std::size_t i = 0; i < r; ++i)
        m_mat[i] = new double[s];
      
    } catch(std::bad_alloc& bad) {
      
      std::cerr << "bad alloc caught: " << bad.what() << std::endl;
      throw;
    }
    m_is_allocated = true;
  }
  m_rows = r; m_cols = s;
}
/// end void allocate

void DenseMatrix::deallocate()
{
  if (!m_is_allocated)
    return;
  for (std::size_t i = 0; i < m_rows; i++) delete [] m_mat[i];
  delete [] m_mat; m_rows = 0; m_cols = 0;
  m_is_allocated = false;
}
/// end void deallocate

DenseMatrix::DenseMatrix(std::size_t p, std::size_t q, const double val)
{
  allocate(p,q);
  for (std::size_t i = 0; i < p; ++i)
    for (std::size_t j = 0; j < q; ++j)
      m_mat[i][j] = val;
}
/// end DenseMatrix

DenseMatrix::DenseMatrix(const DenseMatrix& mat)
{
  deallocate();
  allocate(mat.m_rows,mat.m_cols);
  for (std::size_t i = 0; i < m_rows; i++)
    for (std::size_t j = 0; j < m_cols; j++)
      m_mat[i][j] = mat.m_mat[i][j];
}
/// end copy Constructor DenseMatrix

DenseMatrix::~DenseMatrix()
{
  deallocate();
}
/// end Destructor ~DenseMatrix

double DenseMatrix::operator()(std::size_t i, std::size_t j) const
{ return m_mat[i][j]; }
/// end double operator()

double& DenseMatrix::operator()(std::size_t i, std::size_t j)
{ return m_mat[i][j]; }
/// end double& operator()

const DenseMatrix& DenseMatrix::operator=(const DenseMatrix& mat)
{
  deallocate();
  allocate(mat.m_rows,mat.m_cols);
  for (std::size_t i = 0; i < m_rows; i++)
    for (std::size_t j = 0; j < m_cols; j++)
      m_mat[i][j] = mat.m_mat[i][j];
  return *this;
}
/// end DenseMatrix& operator=

DenseMatrix& DenseMatrix::operator+=(const DenseMatrix& mat)
{
  assert(m_is_allocated);
  if (m_cols != mat.m_cols || m_rows != mat.m_rows)
    throw std::runtime_error("Error in \"DenseMatrix::operator+=\" "
                             "matrix sizes do not agree");
  for (std::size_t i = 0; i < m_rows; i++)
    for (std::size_t j = 0; j < m_cols; j++)
      m_mat[i][j] += mat.m_mat[i][j];
  return *this;
}
/// end DenseMatrix& operator+=

DenseMatrix& DenseMatrix::operator-=(const DenseMatrix& mat)
{
  assert(m_is_allocated);
  if (m_cols != mat.m_cols || m_rows != mat.m_rows)
    throw std::runtime_error("Error in \"DenseMatrix::operator-=\" "
                             "matrix sizes do not agree");
  for (std::size_t i = 0; i < m_rows; i++)
    for (std::size_t j = 0; j < m_cols; j++)
      m_mat[i][j] += mat.m_mat[i][j];
  return *this;
}
/// end DenseMatrix& operator-=

DenseMatrix DenseMatrix::operator+(const DenseMatrix& mat) const
{
  assert(m_is_allocated);
  if (m_cols != mat.m_cols || m_rows != mat.m_rows)
    throw std::runtime_error("Error in \"DenseMatrix::operator+\" "
                             "matrix sizes do not agree");
  DenseMatrix sum(*this);
  sum += mat;
  return sum;
}
/// end DenseMatrix operator+

DenseMatrix DenseMatrix::operator-(const DenseMatrix& mat) const
{
  assert(m_is_allocated);
  if (m_cols != mat.m_cols || m_rows != mat.m_rows)
    throw std::runtime_error("Error in \"DenseMatrix::operator+\" "
                             "matrix sizes do not agree");
  DenseMatrix sum(*this);
  sum -= mat;
  return sum;
}
/// end DenseMatrix& operator-

DenseMatrix DenseMatrix::operator*(const double d)
{
  assert(m_is_allocated);
  DenseMatrix tmp(*this);
  for (std::size_t i = 0; i < m_rows; i++)
    for (std::size_t j = 0; j < m_cols; j++) tmp(i,j) *= d;
  return tmp;
}
/// end DenseMatrix& operator*

DenseVector DenseMatrix::operator*(const DenseVector& vec) const
{
  assert(m_is_allocated);
  if (m_cols != vec.size())
    throw std::runtime_error("Error in \"DenseMatrix::operator*\" "
                             "addition not possible, matrix sizes do not agree");
  
  DenseVector temp_vector(m_rows,0.0);
  for (std::size_t i = 0; i < m_rows; i++)
    for (std::size_t j = 0; j < m_cols; j++) temp_vector[i] += m_mat[i][j]*vec[j];
  return temp_vector;
}
/// end Vector operator*

DenseVector DenseMatrix::trans_mult(const DenseVector& vec) const
{
  assert(m_is_allocated);
  if (m_cols != vec.size())
    throw std::runtime_error("Error in \"DenseMatrix::operator*\" "
                             "addition not possible, matrix sizes do not agree");
  
  DenseVector temp_vector(m_rows,0.0);
  for (std::size_t i = 0; i < m_rows; i++)
    for (std::size_t j = 0; j < m_cols; j++) temp_vector[i] += m_mat[j][i]*vec[j];
  return temp_vector;
}
/// end Vector operator*

void DenseMatrix::transpose(DenseMatrix& transpose_mat)
{
  transpose_mat.allocate(m_rows,m_cols);;
  for (std::size_t i = 0; i < m_rows; i++)
    for (std::size_t j = 0; j < m_cols; j++)
      transpose_mat(j,i) = m_mat[i][j];
}
/// end void Transpose

void DenseMatrix::getDiagonal(DenseVector& vec) const
{
  vec.allocate(m_rows);
  for (std::size_t i = 0; i < m_rows; i++)
    vec(i) = m_mat[i][i];
}
/// end void Transpose

}
// end namespace math8650
