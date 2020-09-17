
#include "SymmetricMatrix.hpp"
#include "DenseVector.hpp"

namespace math8650
{

void SymmetricMatrix::allocate(std::size_t r, std::size_t s)
{
  if (!((r > 0) && (s > 0) && (r==s)))
  {
    m_is_allocated = false;
  } else {
  
    try {
      
      m_mat = new double*[r];
      // create a matrix only for the lower triangular part
      for (std::size_t i = 0; i < r; ++i)
        m_mat[i] = new double[i+1];
      
    } catch(std::bad_alloc& bad) {
      
      std::cerr << "bad alloc caught: " << bad.what() << std::endl;
      throw;
    }
    m_is_allocated = true;
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

SymmetricMatrix::SymmetricMatrix(std::size_t p, std::size_t q, const double val) : Matrix(p,q)
{
  allocate(p,q);
  for (std::size_t i = 0; i < p; ++i)
    for (std::size_t j = 0; j < i+1; ++j)
      m_mat[i][j] = val;
}
/// end SymmetricMatrix

SymmetricMatrix::SymmetricMatrix(const SymmetricMatrix& mat) : Matrix(mat.m_rows,mat.m_cols)
{
  deallocate();
  allocate(mat.m_rows,mat.m_cols);
  for (std::size_t i = 0; i < m_rows; i++)
    for (std::size_t j = 0; j < m_cols; j++)
      m_mat[i][j] = mat.m_mat[i][j];
}
/// end copy Constructor SymmetricMatrix

SymmetricMatrix::~SymmetricMatrix()
{
  deallocate();
}
/// end Destructor ~SymmetricMatrix

double SymmetricMatrix::operator()(std::size_t i, std::size_t j) const
{
  if (i >= j)	
    return m_mat[i][j];
  else
    return m_mat[j][i];
}
/// end double operator()

double& SymmetricMatrix::operator()(std::size_t i, std::size_t j)
{
  if (i >= j)
    return m_mat[i][j];
  else 
    throw std::runtime_error("Error in \"SymmetricMatrix::operator()\" accessing unallocated upper triangular entries");
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
    for (std::size_t j = 0; j < m_cols; j++) temp_vector[i] += (*this)(i,j)*vec[j];
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
