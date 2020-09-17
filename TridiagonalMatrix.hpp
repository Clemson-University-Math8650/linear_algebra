#ifndef __Math8650_TridiagonalMatrix_HPP__
#define __Math8650_TridiagonalMatrix_HPP__

#include <cassert>
#include <iostream>

#include "DenseVector.hpp"


namespace math8650
{

class TridiagonalMatrix
{
public:

  TridiagonalMatrix()
  {
    m_rows = 0;
    m_cols = 0;
    m_is_allocated = false;
  }

  TridiagonalMatrix(std::size_t nrows, std::size_t ncols, const double val = 0.0);
  virtual ~TridiagonalMatrix();
  void allocate(std::size_t nrows, std::size_t ncols);
  void deallocate();
  double& operator()(std::size_t i, std::size_t j);
  double operator()(std::size_t i, std::size_t j) const;
  DenseVector operator*(const DenseVector& vec) const;
  DenseVector trans_mult(const DenseVector& vec) const;
  void getDiagonal(DenseVector& vec) const;
  inline std::size_t numRows() const { return m_rows; }
  inline std::size_t numCols() const { return m_cols; }

  friend std::ostream& operator<<(std::ostream& out, const TridiagonalMatrix& mat)
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
  double* m_diag;
  double* m_above;
  double* m_below;
  // number of rows
  std::size_t m_rows;
  // number of columns
  std::size_t m_cols;
  // allocated flag
  bool m_is_allocated;  
  
};

}
// end namespace math8650

#endif /** __Math8650_TridiagonalMatrix_HPP__ */
