#ifndef MAT_HPP
#define MAT_HPP
#include "math.hpp"

class Mat3
{
 public:
  /**
   * The identity matrix.
   */
  static const Mat3 Identity;

  /**
   * The zero matrix.
   */
  static const Mat3 Zero;

  static const int DIM = 3;
  static const int SIZE = DIM*DIM;
  /**@brief COLUMN Major format !!!!!*/
  explicit Mat3(real_t r[SIZE]);
  Mat3() {}
  Mat3(const Mat3 & in){for(int ii=0;ii<9;ii++){m[ii]=in.m[ii];}}
  /**@brief ROW Major format!!!!*/
  Mat3(real_t m00, real_t m10, real_t m20,
       real_t m01, real_t m11, real_t m21,
       real_t m02, real_t m12, real_t m22);
  /**@brief each ROW is passed in as a vector*/
  Mat3(const Vec3 & v1, const Vec3&v2,const Vec3 & v3);
  void inverse(Mat3& rv) const;
  //original matrix is destroyed
  void gauss_pivot(Mat3 & rv);
  union{
    real_t m[SIZE];
    real_t _m[DIM][DIM];
  };
  void transpose();
  void mult(const Mat3 & rv, Mat3& result);
  /**@TODO: not working yet, need to figure out
  how to use permutation matrix p and q*/
  void invert_full(Mat3 & rv);
private:
  void h_solve(real_t *x, int *p, int *q);
  void h_pivot_decomp(int *p, int *q);
};

#endif
