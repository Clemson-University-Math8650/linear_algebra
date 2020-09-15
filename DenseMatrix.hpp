#ifndef __Math8650_DenseMatrix_HPP__
#define __Math8650_DenseMatrix_HPP__

#include <cassert>
#include <iostream>

#include "Matrix.hpp"


namespace math8650
{

class DenseMatrix : public Matrix
{
public:

  DenseMatrix() : Matrix()
  {
  }

  DenseMatrix(std::size_t p, std::size_t q, const double val = 0.0);
  
  DenseMatrix(const DenseMatrix& w);

//  const DenseMatrix& operator=(const DenseMatrix& mat);
  
  virtual ~DenseMatrix();
  
  virtual void allocate(std::size_t r, std::size_t s) override;
  
  virtual void deallocate() override;
  
  virtual double& operator()(std::size_t i, std::size_t j) override;
  
  virtual double operator()(std::size_t i, std::size_t j) const override;

/*  virtual DenseMatrix& operator+=(const DenseMatrix& mat) override;
  virtual DenseMatrix& operator-=(const DenseMatrix& B) override;
  virtual DenseMatrix  operator+(const DenseMatrix& B) const override;
  virtual DenseMatrix  operator-(const DenseMatrix& B) const override;
  virtual DenseMatrix  operator*(const double num)  override;*/
  virtual DenseVector  operator*(const DenseVector& vec) const override;
  virtual DenseVector  trans_mult(const DenseVector& vec) const override;

  virtual void getDiagonal(DenseVector& vec) const override;
  
private:
  
  // matrix entries
  double** m_mat;

};

}
// end namespace math8650

#endif /** __Math8650_DenseMatrix_HPP__ */
