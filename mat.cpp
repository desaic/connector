#include "mat.hpp"
#include "string.h"
#include <stdlib.h>
void Mat3::h_pivot_decomp(int *p, int *q) {
  int i,j,k;
  int n=DIM;
  int pi,pj,tmp;
  real_t max;
  real_t ftmp;
  for (k=0; k<n; k++) {
    pi=-1,pj=-1,max=0.0;
    //find pivot in submatrix a(k:n,k:n)
    for (i=k; i<n; i++) {
      for (j=k; j<n; j++) {
        if (std::abs(_m[i][j])>max) {
          max = std::abs(_m[i][j]);
          pi=i;
          pj=j;
        }
      }
    }
    //Swap Row
    tmp=p[k];
    p[k]=p[pi];
    p[pi]=tmp;
    for (j=0; j<n; j++) {
      ftmp=_m[k][j];
      _m[k][j]=_m[pi][j];
      _m[pi][j]=ftmp;
    }
    //Swap Col
    tmp=q[k];
    q[k]=q[pj];
    q[pj]=tmp;
    for (i=0; i<n; i++) {
      ftmp=_m[i][k];
      _m[i][k]=_m[i][pj];
      _m[i][pj]=ftmp;
    }
    //END PIVOT

    //check pivot size and decompose
    if ((std::abs(_m[k][k])>1e-20)) {
      for (i=k+1; i<n; i++) {
        //Column normalisation
        _m[i][k]/=_m[k][k];
        ftmp=_m[i][k];
        for (j=k+1; j<n; j++) {
          //a(ik)*a(kj) subtracted from lower right submatrix elements
          _m[i][j]-=(ftmp*_m[k][j]);
        }
      }
    }
    //END DECOMPOSE
  }
}


void Mat3::h_solve(real_t *x, int *p, int *q) {
  //forward substitution; see  Golub, Van Loan 96
  //And see http://www.cs.rutgers.edu/~richter/cs510/completePivoting.pdf
  int i,ii=0,j;
  real_t ftmp;
  real_t xtmp[DIM];
  //Swap rows (x=Px)
  for (i=0; i<DIM; i++) {
    xtmp[i]=x[p[i]]; //value that should be here
  }
  //Lx=x
  for (i=0; i<DIM; i++) {
    ftmp=xtmp[i];
    if (ii != 0)
      for (j=ii-1; j<i; j++)
        ftmp-=_m[i][j]*xtmp[j];
    else if (ftmp!=0.0)
      ii=i+1;
    xtmp[i]=ftmp;
  }
  //backward substitution
  //partially taken from Sourcebook on Parallel Computing p577
  //solves Uy=z
  xtmp[DIM-1]/=_m[DIM-1][DIM-1];
  for (i=DIM-2; i>=0; i--) {
    ftmp=xtmp[i];
    for (j=i+1; j<DIM; j++) {
      ftmp-=_m[i][j]*xtmp[j];
    }
    xtmp[i]=(ftmp)/_m[i][i];
  }

  for (i=0; i<DIM; i++) {
    x[i]=xtmp[q[i]];
  }
}

void Mat3::invert_full(Mat3 & rv)
{
  int p[DIM];
  int q[DIM];
  real_t b[DIM];
  for(int ii=0; ii<DIM; ii++) {
    p[ii]=ii;
    q[ii]=ii;
    b[ii]=0;
  }
  h_pivot_decomp(p,q);
  for(int ii=0; ii<3; ii++) {
    b[ii]=1;
    h_solve(b,p,q);
    for(int jj=0; jj<3; jj++) {
      rv._m[jj][ii]=b[q[jj]];
      b[q[jj]]=0;
    }
  }
}

void Mat3::transpose()
{
  for(int ii=0; ii<DIM; ii++) {
    for(int jj=0; jj<DIM; jj++) {
      real_t tmp = _m[ii][jj];
      _m[ii][jj]=_m[jj][ii];
      _m[jj][ii]=tmp;
    }
  }
}

void Mat3::gauss_pivot(Mat3 & rv)
{
  int n=3;
  real_t * A=m;
  real_t * AInverse=rv.m;
  int i, j, iPass, imx, icol, irow;
  real_t det, temp, pivot, factor=0;
  real_t* ac = (real_t*)calloc(n*n, sizeof(real_t));
  det = 1;

  for (i = 0; i < n; i++)
  {
    for (j = 0; j < n; j++)
    {
      AInverse[n*i+j] = 0;
      ac[n*i+j] = A[n*i+j];
    }
    AInverse[n*i+i] = 1;
  }

  // The current pivot row is iPass.
  // For each pass, first find the maximum element in the pivot column.
  for (iPass = 0; iPass < n; iPass++)
  {
    imx = iPass;
    for (irow = iPass; irow < n; irow++)
    {
      if (std::abs(A[n*irow+iPass]) > std::abs(A[n*imx+iPass])) imx = irow;
    }
    // Interchange the elements of row iPass and row imx in both A and AInverse.
    if (imx != iPass)
    {
      for (icol = 0; icol < n; icol++)
      {
        temp = AInverse[n*iPass+icol];
        AInverse[n*iPass+icol] = AInverse[n*imx+icol];
        AInverse[n*imx+icol] = temp;

        if (icol >= iPass)
        {
          temp = A[n*iPass+icol];
          A[n*iPass+icol] = A[n*imx+icol];
          A[n*imx+icol] = temp;
        }
      }
    }

    // The current pivot is now A[iPass][iPass].
    // The determinant is the product of the pivot elements.
    pivot = A[n*iPass+iPass];
    det = det * pivot;
    if (det == 0)
    {
      free(ac);
      return;
    }

    for (icol = 0; icol < n; icol++)
    {
      // Normalize the pivot row by dividing by the pivot element.
      AInverse[n*iPass+icol] = AInverse[n*iPass+icol] / pivot;
      if (icol >= iPass) A[n*iPass+icol] = A[n*iPass+icol] / pivot;
    }

    for (irow = 0; irow < n; irow++)
      // Add a multiple of the pivot row to each row.  The multiple factor
      // is chosen so that the element of A on the pivot column is 0.
    {
      if (irow != iPass) factor = A[n*irow+iPass];
      for (icol = 0; icol < n; icol++)
      {
        if (irow != iPass)
        {
          AInverse[n*irow+icol] -= factor * AInverse[n*iPass+icol];
          A[n*irow+icol] -= factor * A[n*iPass+icol];
        }
      }
    }
  }
  free(ac);
}
void Mat3::mult(const Mat3 & rv, Mat3& result)
{
  for(int ii=0; ii<3; ii++) {
    for(int jj=0; jj<3; jj++) {
      result._m[ii][jj]=0;
      for(int kk=0; kk<3; kk++) {
        result._m[ii][jj]+=_m[ii][kk]*rv._m[kk][jj];
      }
    }
  }

}
void Mat3::inverse(Mat3& rv) const
{
  rv._m[0][0] = _m[1][1]*_m[2][2] - _m[1][2]*_m[2][1];
  rv._m[0][1] = _m[0][2]*_m[2][1] - _m[0][1]*_m[2][2];
  rv._m[0][2] = _m[0][1]*_m[1][2] - _m[0][2]*_m[1][1];
  rv._m[1][0] = _m[1][2]*_m[2][0] - _m[1][0]*_m[2][2];
  rv._m[1][1] = _m[0][0]*_m[2][2] - _m[0][2]*_m[2][0];
  rv._m[1][2] = _m[0][2]*_m[1][0] - _m[0][0]*_m[1][2];
  rv._m[2][0] = _m[1][0]*_m[2][1] - _m[1][1]*_m[2][0];
  rv._m[2][1] = _m[0][1]*_m[2][0] - _m[0][0]*_m[2][1];
  rv._m[2][2] = _m[0][0]*_m[1][1] - _m[0][1]*_m[1][0];

  real_t det = _m[0][0]*rv._m[0][0] +
               _m[0][1]*rv._m[1][0] +
               _m[0][2]*rv._m[2][0];

  real_t invdet = 1.0 / det;
  for (int i = 0; i < SIZE; i++)
    rv.m[i] *= invdet;
}

const Mat3 Mat3::Identity = Mat3( 1, 0, 0,
                                  0, 1, 0,
                                  0, 0, 1 );

const Mat3 Mat3::Zero = Mat3( 0, 0, 0,
                              0, 0, 0,
                              0, 0, 0 );


Mat3::Mat3(real_t r[SIZE])
{
  memcpy(m ,r, sizeof r);
}

Mat3::Mat3(real_t m00, real_t m10, real_t m20,
           real_t m01, real_t m11, real_t m21,
           real_t m02, real_t m12, real_t m22)
{
  _m[0][0] = m00;
  _m[1][0] = m10;
  _m[2][0] = m20;
  _m[0][1] = m01;
  _m[1][1] = m11;
  _m[2][1] = m21;
  _m[0][2] = m02;
  _m[1][2] = m12;
  _m[2][2] = m22;
}

Mat3::Mat3(const Vec3& v1,const Vec3& v2,const Vec3& v3)
{
  for(int ii=0; ii<3; ii++) {
    _m[ii][0]=v1.get(ii);
    _m[ii][1]=v2.get(ii);
    _m[ii][2]=v3.get(ii);
  }
}
