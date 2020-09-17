#ifndef __Math8650_SymmetricMatrix_HPP__
#define __Math8650_SymmetricMatrix_HPP__

#include <cassert>
#include <iostream>

#include "Matrix.hpp"


namespace math8650
{

class SymmetricMatrix : public Matrix
{
public:

  SymmetricMatrix() : Matrix()
  {
  }

  SymmetricMatrix(std::size_t p, std::size_t q, const double val = 0.0);
  
  SymmetricMatrix(const SymmetricMatrix& w);
  
  virtual ~SymmetricMatrix();
  
  virtual void allocate(std::size_t r, std::size_t s) override;
  
  virtual void deallocate() override;
  
  virtual double& operator()(std::size_t i, std::size_t j) override;
  
  virtual double operator()(std::size_t i, std::size_t j) const override;

  virtual DenseVector  operator*(const DenseVector& vec) const override;
  virtual DenseVector  trans_mult(const DenseVector& vec) const override;

  virtual void getDiagonal(DenseVector& vec) const override;
  
private:
  
  // matrix entries
  double** m_mat;

};

}
// end namespace math8650

#endif /** __Math8650_SymmetricMatrix_HPP__ */
