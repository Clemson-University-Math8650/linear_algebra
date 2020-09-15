#ifndef __Math8650_Matrix_HPP__
#define __Math8650_Matrix_HPP__

#include <iostream>

namespace math8650
{

class DenseVector;  
  
/**
 * \brief Abstract base class for matrices
 *
 */
class Matrix
{
public:

  Matrix()
  {
    m_rows = 0;
    m_cols = 0;
    m_is_allocated = false;
  }

  Matrix(std::size_t p, std::size_t q) : m_rows(p), m_cols(q) {}

  virtual ~Matrix() {}
  
  virtual void allocate(std::size_t r, std::size_t s) = 0;
  
  virtual void deallocate() = 0;
  
  virtual double& operator()(std::size_t i, std::size_t j) = 0;
  
  virtual double operator()(std::size_t i, std::size_t j) const = 0;

  inline std::size_t numRows() const { return m_rows; }
  
  inline std::size_t numCols() const { return m_cols; }
  
/*  virtual Matrix& operator+=(const Matrix& mat) = 0;
  virtual Matrix& operator-=(const Matrix& B) = 0;
  virtual Matrix  operator+(const Matrix& B) const = 0;
  virtual Matrix  operator-(const Matrix& B) const = 0;
  virtual Matrix  operator*(const double num) = 0;*/
  virtual DenseVector  operator*(const DenseVector& vec) const = 0;
  virtual DenseVector  trans_mult(const DenseVector& vec) const = 0;

  virtual void getDiagonal(DenseVector& vec) const = 0;
  
  friend std::ostream& operator<<(std::ostream& out, const Matrix& mat)
  {
    
    for (std::size_t i = 0; i < mat.m_rows; ++i)
    {
      for (std::size_t j = 0; j < mat.m_cols; ++j)
        out << mat(i,j) << " ";
      out << std::endl;
    }
    return out;
  }
  
protected:
  
  // number of rows
  std::size_t m_rows;
  // number of columns
  std::size_t m_cols;
  // allocated flag
  bool m_is_allocated;

};

}
// end namespace math8650

#endif /** __Math8650_Matrix_HPP__ */
