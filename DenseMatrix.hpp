#ifndef __Math8650_DenseMatrix_HPP__
#define __Math8650_DenseMatrix_HPP__

#include <cassert>
#include <iostream>

namespace math8650
{

class DenseVector;  
  
class DenseMatrix
{
public:

  DenseMatrix()
  {
    m_rows = 0;
    m_cols = 0;
    m_is_allocated = false;
  }

  DenseMatrix(std::size_t p, std::size_t q, const double val = 0.0);
  
  DenseMatrix(const DenseMatrix& w);

  const DenseMatrix& operator=(const DenseMatrix& mat);
  
  virtual ~DenseMatrix();
  
  void allocate(std::size_t r, std::size_t s);
  
  void deallocate();
  
  double& operator()(std::size_t i, std::size_t j);
  
  double operator()(std::size_t i, std::size_t j) const;

  inline std::size_t numRows() const { return m_rows; }
  
  inline std::size_t numCols() const { return m_cols; }
  
  DenseMatrix& operator+=(const DenseMatrix& mat);
  DenseMatrix& operator-=(const DenseMatrix& B);
  DenseMatrix  operator+(const DenseMatrix& B) const;
  DenseMatrix  operator-(const DenseMatrix& B) const;
  DenseMatrix  operator*(const double num);
  DenseVector  operator*(const DenseVector& vec) const;
  DenseVector  trans_mult(const DenseVector& vec) const;

  void getDiagonal(DenseVector& vec) const;
  void transpose(DenseMatrix& trans_mat);
  
  friend std::ostream& operator<<(std::ostream& out, const DenseMatrix& mat)
  {
    
    for (std::size_t i = 0; i < mat.m_rows; ++i)
    {
      for (std::size_t j = 0; j < mat.m_cols; ++j)
        out << mat.m_mat[i][j] << " ";
      out << std::endl;
    }
    return out;
  }
  
private:
  
  // matrix entries
  double** m_mat;
  // number of rows
  std::size_t m_rows;
  // number of columns
  std::size_t m_cols;
  // allocated flag
  bool m_is_allocated;

};

}
// end namespace math8650

#endif /** __Math8650_DenseMatrix_HPP__ */
