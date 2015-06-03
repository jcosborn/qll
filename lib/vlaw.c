typedef struct {
  vec r,i;
} Cmplx;
typedef struct {
  Cmplx c0,c1,c2;
} Vec3;
typedef struct {
  Vec3 s0,s1;
} HVec3;
typedef struct {
  Vec3 s0,s1,s2,s3;
} DVec3;
typedef struct { // stored in column major format
  Vec3 v0,v1,v2;
} Mat3;

#include "vlaw0.c"

void
P(H_eq_spinproj_D)(real *rr, real *vv, int dir, int sign, int i0, int i1)
{
  HVec3 *restrict r = (HVec3 *)rr;
  DVec3 *restrict v = (DVec3 *)vv;
  for(int i=i0; i<i1; i++) {
    HVec3 *restrict ri = &r[i];
    DVec3 *restrict vi = &v[i];
    sprojD3(ri, vi, dir, sign);
  }
}

void
P(D_eq_spinrecon_H)(real *rr, real *vv, int dir, int sign, int i0, int i1)
{
  DVec3 *restrict r = (DVec3 *)rr;
  HVec3 *restrict v = (HVec3 *)vv;
  for(int i=i0; i<i1; i++) {
    DVec3 *restrict ri = &r[i];
    HVec3 *restrict vi = &v[i];
    sreconD3(ri, vi, dir, sign);
  }
}

#if 0
void
P(D_peq_spinproj_M_times_xDp)(real *rr, real *mm, int *idx, int perm,
			      real *vv, int dir, int sign, int i0, int i1)
{
  DVec3 *restrict r = (DVec3 *)rr;
  Mat3 *restrict m = (Mat3 *)mm;
  DVec3 *restrict v = (DVec3 *)vv;
  for(int i=i0; i<i1; i++) {
    int k1 = idx[i];
    if(k1>=0) {
      DVec3 *restrict ri = &r[i];
      Mat3 *restrict mi = &m[i];
      DVec3 *restrict vi = &v[k1];
      PeqSprojM3timesD3(ri, mi, vi, dir, sign);
    } else if(k1+2<=0) {
      int k2 = -(k1+2);
      DVec3 *restrict ri = &r[i];
      Mat3 *restrict mi = &m[i];
      DVec3 *restrict vi = &v[k2];
      PeqSprojM3timesD3p(ri, mi, vi, dir, sign, perm);
    }
  }
}
#endif

void
P(D_meq_spinproj_M_times_xDp)(real *rr, real *mm, int *idx, int perm,
			      real *vv, int dir, int sign, int i0, int i1)
{
  DVec3 *restrict r = (DVec3 *)rr;
  Mat3 *restrict m = (Mat3 *)mm;
  DVec3 *restrict v = (DVec3 *)vv;
  for(int i=i0; i<i1; i++) {
    int k1 = idx[i];
    if(k1>=0) {
      DVec3 *restrict ri = &r[i];
      Mat3 *restrict mi = &m[i];
      DVec3 *restrict vi = &v[k1];
      MeqSprojM3timesD3(ri, mi, vi, dir, sign);
    } else if(k1+2<=0) {
      int k2 = -(k1+2);
      DVec3 *restrict ri = &r[i];
      Mat3 *restrict mi = &m[i];
      DVec3 *restrict vi = &v[k2];
      MeqSprojM3timesD3p(ri, mi, vi, dir, sign, perm);
    }
  }
}
