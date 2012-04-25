#include "math.hpp"

real_t half_tan(real_t cosine)
{
    real_t sine = std::sqrt(1-cosine*cosine);
    real_t halftangent = sine/(1+cosine);
    return halftangent;
}
Vec3 operator*(real_t c, const Vec3 & v)
{
  return v*c;
}
