#ifndef __Math8650_LowerTriangularMatrix_HPP__
#define __Math8650_LowerTriangularMatrix_HPP__

#include <cassert>
#include <iostream>

#include "DenseVector.hpp"


namespace math8650
{

class LowerTriangularMatrix
{
public:

  LowerTriangularMatrix()
  {
    m_rows = 0;
    m_cols = 0;
    m_is_allocated = false;
  }

  LowerTriangularMatrix(std::size_t nrows, std::size_t ncols, const double val = 0.0);
  virtual ~LowerTriangularMatrix();
  void allocate(std::size_t nrows, std::size_t ncols);
  void deallocate();
  double& operator()(std::size_t i, std::size_t j);
  double operator()(std::size_t i, std::size_t j) const;
  DenseVector operator*(const DenseVector& vec) const;
  DenseVector trans_mult(const DenseVector& vec) const;
  void getDiagonal(DenseVector& vec) const;
  inline std::size_t numRows() const { return m_rows; }
  inline std::size_t numCols() const { return m_cols; }

  friend std::ostream& operator<<(std::ostream& out, const LowerTriangularMatrix& mat)
  {
    for (std::size_t i = 0; i < mat.m_rows; ++i)
    {
      for (std::size_t j = 0; j < mat.m_cols; ++j)
        out << mat(i,j) << " ";
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

#endif /** __Math8650_LowerTriangularMatrix_HPP__ */
