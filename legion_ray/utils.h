#ifndef UTILS_H
#define UTILS_H

typedef unsigned char uchar;
typedef unsigned int uint;
namespace utils
{
template<typename T>
T min(T a, T b)
{
  return a < b?a : b;
}

template <typename T>
T max(T a, T b)
{
  return a > b? a:b;
}
}
#endif
