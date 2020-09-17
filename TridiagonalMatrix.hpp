#ifndef __Math8650_TridiagonalMatrix_HPP__
#define __Math8650_TridiagonalMatrix_HPP__

#include <cassert>
#include <iostream>

#include "Matrix.hpp"


namespace math8650
{

class TridiagonalMatrix : public Matrix
{
public:

  TridiagonalMatrix() : Matrix()
  {
  }

  TridiagonalMatrix(std::size_t p, std::size_t q, const double val = 0.0);
  
 
  virtual ~TridiagonalMatrix();
  
  virtual void allocate(std::size_t r, std::size_t s) override;
  
  virtual void deallocate() override;
  
  virtual double& operator()(std::size_t i, std::size_t j) override;
  
  virtual double operator()(std::size_t i, std::size_t j) const override;

  virtual DenseVector  operator*(const DenseVector& vec) const override;
  virtual DenseVector  trans_mult(const DenseVector& vec) const override;

  virtual void getDiagonal(DenseVector& vec) const override;
  
private:
  
  // matrix entries
  double* m_diag;
  double* m_above;
  double* m_below;

};

}
// end namespace math8650

#endif /** __Math8650_TridiagonalMatrix_HPP__ */
