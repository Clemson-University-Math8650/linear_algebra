
#ifndef __Math8650_DenseVector_HPP__
#define __Math8650_DenseVector_HPP__

#include <cmath>
#include <cassert>
#include <iostream>
#include <algorithm>

namespace math8650 
{

class DenseVector
{
public:
  
  DenseVector() : m_size(0) {}
  DenseVector(std::size_t n, const double val = 0.0);
  DenseVector(const DenseVector& old_vector);
  virtual ~DenseVector();
  
  void allocate(std::size_t r, const double val = 0.0);
  
  void deallocate();

  DenseVector& operator=(const DenseVector& vec);
  DenseVector& operator+=(const DenseVector& vec);
  DenseVector& operator-=(const DenseVector& vec);
  DenseVector  operator+(const DenseVector& vec) const;
  DenseVector  operator-(const DenseVector& vec) const;
  DenseVector& operator*=(const double num);

  double operator[](std::size_t index) const;
  double& operator[](std::size_t index);
  double operator()(std::size_t index) const;
  double& operator()(std::size_t index);

  void clear();
  std::size_t size() const;
  double norm() const;
  double max() const;
  double min() const;
  
  double operator*(const DenseVector& vec) const;
  
  friend DenseVector operator*(const double num, const DenseVector& vec)
  {
    DenseVector temp(vec);
    temp *= num;
    return temp;
  }
  /// end friend DenseVector operator*
  
  friend std::ostream& operator<<(std::ostream& out, const DenseVector& vec)
  {
    for (std::size_t i = 0; i < vec.size(); ++i)
      out << vec[i] << std::endl;
    return out;
  }  

private:
  
  std::size_t m_size;
  double* m_vec;
  bool m_is_allocated = false;
    
};
/// end class DenseVector

double norm(const DenseVector& r);
double dot(const DenseVector& a, const DenseVector& b);

}
// end namespace math8650

#endif /** __Math8650_DenseVector_HPP__ */
