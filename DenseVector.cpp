
#include "DenseVector.hpp"

namespace math8650
{

void DenseVector::allocate(std::size_t n, const double val)
{
  if (n == 0)
  {
    m_is_allocated = false;
  } else {
  
    try {
      
      m_vec = new double[n];
      
    } catch(std::bad_alloc& bad) {
      
      std::cerr << "bad alloc caught: " << bad.what() << std::endl;
      throw;
      
    }
    m_is_allocated = true;
    
    for (std::size_t i = 0; i < n; ++i) m_vec[i] = val;
    
  }
  m_size = n;
}
/// end void allocate

void DenseVector::deallocate()
{
  if (!m_is_allocated)
    return;
  delete [] m_vec; m_size = 0;
  m_is_allocated = false;
}
/// end void deallocate

DenseVector::DenseVector(std::size_t n, const double val)
{
  allocate(n,val);
}
/// end Constructor DenseVector

DenseVector::DenseVector(const DenseVector& vec)
{
  allocate(vec.size());
  for (std::size_t i = 0; i < vec.size(); ++i)
    m_vec[i] = vec(i);
}
/// end Constructor DenseVector

DenseVector::~DenseVector()
{
  if (m_is_allocated)
    deallocate();
}
/// end Destructor ~DenseVector

DenseVector& DenseVector::operator=(const DenseVector& vec)
{
//   if (m_is_allocated)
//     deallocate();
  allocate(vec.size());
  for (std::size_t i = 0; i < vec.size(); ++i)
    m_vec[i] = vec(i);
  return *this;
}
/// end DenseVector& operator=

double DenseVector::operator[](std::size_t index) const
{ return m_vec[index]; }
/// end double operator[]

double& DenseVector::operator[](std::size_t index)
{ return m_vec[index]; }
/// end double& operator[]

double DenseVector::operator()(std::size_t index) const
{ return m_vec[index]; }
/// end double operator()

double& DenseVector::operator()(std::size_t index)
{ return m_vec[index]; }
/// end double& operator()

std::size_t DenseVector::size() const
{ return m_size; }
/// end int getsize

void DenseVector::clear()
{ deallocate(); }
/// end void clear

double DenseVector::operator*(const DenseVector& vec) const
{
  assert(m_is_allocated);
  if (m_size != vec.m_size )
    throw std::runtime_error("Error in \"DenseVector::operator*\" vector sizes not equal");
  double sum = 0.0;
  for (std::size_t i = 0; i < m_size; ++i)
  {
    sum += m_vec[i] * vec.m_vec[i];
  }
  return sum;
}
/// end double operator*

DenseVector& DenseVector::operator*=(const double num)
{
  assert(m_is_allocated);
  for (std::size_t i = 0; i < m_size; ++i) m_vec[i] *= num;
  return *this;
}
/// end DenseVector& operator*=

DenseVector DenseVector::operator+(const DenseVector& vec) const
{
  assert(m_is_allocated);
  if (m_size != vec.m_size)
    throw std::runtime_error("Error in \"DenseVector::operator+\" vector sizes not equal");
  DenseVector temp (*this); temp += vec;
  return temp;
}
/// end DenseVector operator+

DenseVector DenseVector::operator-(const DenseVector& vec) const
{
  assert(m_is_allocated);
  if (m_size != vec.m_size)
    throw std::runtime_error("Error in \"DenseVector::operator-\" vector sizes not equal");
  DenseVector temp (*this); temp -= vec;
  return temp;
}
/// end DenseVector operator-

DenseVector& DenseVector::operator+=(const DenseVector& vec)
{
  assert(m_is_allocated);
  if (m_size != vec.m_size)
    throw std::runtime_error("Error in \"DenseVector::operator+=\" vector sizes not equal");
  for (std::size_t i = 0; i < m_size; ++ i ) m_vec[i] += vec[i];
  return *this;
}
/// end DenseVector& operator+=

DenseVector& DenseVector::operator-=(const DenseVector& vec)
{
  assert(m_is_allocated);
  if (m_size != vec.m_size)
    throw std::runtime_error("Error in \"DenseVector::operator-=\" vector sizes not equal");
  for (std::size_t i = 0; i < m_size; ++i) m_vec[i] -= vec[i];
  return *this;
}
/// end DenseVector& operator-=

double DenseVector::norm() const
{
  double s (0.0);
  for (std::size_t i(0); i < m_size; ++i) 
    s += m_vec[i]*m_vec[i];
  return std::sqrt(s);
}
/// end double norm

double DenseVector::max() const 
{
  return * std::max_element(&(m_vec[0]),&(m_vec[m_size-1]));
}
/// end double max() const

double DenseVector::min() const 
{
  return * std::min_element(&(m_vec[0]),&(m_vec[m_size-1]));
}
/// end double min() const

double norm(const DenseVector& r)
{ return r.norm(); }
/// end double norm

double dot(const DenseVector& a, const DenseVector& b)
{ return a*b; }
/// end double dot

}
// end namespace math8650
