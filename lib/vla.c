#include <stdio.h>

typedef struct {
  vec r,i;
} Cmplx;
typedef struct {
  Cmplx c0,c1,c2;
} Vec3;
typedef struct { // stored in column major format
  Vec3 v0,v1,v2;
} Mat3;

static inline void
Copy(vec *x, vec *y, int n)
{
  for(int i=0; i<n; i++) {
    STNR(&(x[i]), LD(&(y[i])));
  }
}
static inline void
CopyC(Cmplx *x, Cmplx *y)
{
  STNR(&(x->r), LD(&(y->r)));
  STNR(&(x->i), LD(&(y->i)));
}
static inline void
CopyV3(Vec3 *x, Vec3 *y)
{
  CopyC(&(x->c0), &(y->c0));
  CopyC(&(x->c1), &(y->c1));
  CopyC(&(x->c2), &(y->c2));
}
#if VLEN > 1
static inline void
Perm1C(Cmplx *x, Cmplx *y)
{
  vecr yr = LD(&(y->r));
  vecr xr = perm1(yr);
  STNR(&(x->r), xr);
  vecr yi = LD(&(y->i));
  vecr xi = perm1(yi);
  STNR(&(x->i), xi);
}
static inline void
Perm1V3(Vec3 *x, Vec3 *y)
{
  Perm1C(&(x->c0), &(y->c0));
  Perm1C(&(x->c1), &(y->c1));
  Perm1C(&(x->c2), &(y->c2));
}
#endif
#if VLEN > 2
static inline void
Perm2C(Cmplx *x, Cmplx *y)
{
  vecr yr = LD(&(y->r));
  vecr xr = perm2(yr);
  STNR(&(x->r), xr);
  vecr yi = LD(&(y->i));
  vecr xi = perm2(yi);
  STNR(&(x->i), xi);
}
static inline void
Perm2V3(Vec3 *x, Vec3 *y)
{
  Perm2C(&(x->c0), &(y->c0));
  Perm2C(&(x->c1), &(y->c1));
  Perm2C(&(x->c2), &(y->c2));
}
#endif
#if VLEN > 4
static inline void
Perm4C(Cmplx *x, Cmplx *y)
{
  vecr yr = LD(&(y->r));
  vecr xr = perm4(yr);
  STNR(&(x->r), xr);
  vecr yi = LD(&(y->i));
  vecr xi = perm4(yi);
  STNR(&(x->i), xi);
}
static inline void
Perm4V3(Vec3 *x, Vec3 *y)
{
  Perm4C(&(x->c0), &(y->c0));
  Perm4C(&(x->c1), &(y->c1));
  Perm4C(&(x->c2), &(y->c2));
}
#endif
#if VLEN > 8
static inline void
Perm8C(Cmplx *x, Cmplx *y)
{
  vecr yr = LD(&(y->r));
  vecr xr = perm8(yr);
  STNR(&(x->r), xr);
  vecr yi = LD(&(y->i));
  vecr xi = perm8(yi);
  STNR(&(x->i), xi);
}
static inline void
Perm8V3(Vec3 *x, Vec3 *y)
{
  Perm8C(&(x->c0), &(y->c0));
  Perm8C(&(x->c1), &(y->c1));
  Perm8C(&(x->c2), &(y->c2));
}
#endif

#include "vla0.c"

void
P(X_veq_perm_xX)(real *rr, int perm, int *idx, real *vv,
		 int nelem, int off, int len)
{
  vec *restrict r = (vec *restrict) rr;
  vec *restrict v = (vec *restrict) vv;
#define LOOP(pf) for(int i=0; i<len; i++) {		\
    int k = idx[off+i];					\
    vec *restrict ri = &r[nelem*(off+i)];		\
    vec *restrict vi = &v[nelem*k];			\
    pf(ri, vi, nelem);					\
  }
  switch(perm) {
  case 0: LOOP(Copy); break;
#if VLEN > 1
  case 1: LOOP(Perm1); break;
#endif
#if VLEN > 2
  case 2: LOOP(Perm2); break;
#endif
#if VLEN > 4
  case 4: LOOP(Perm4); break;
#endif
#if VLEN > 8
  case 8: LOOP(Perm8); break;
#endif
  }
#undef LOOP
}

void
P(X_eq_pack_xX)(real *rr, int pack, int *idx, real *vv,
		int nelem, int i0, int i1)
{
  if(pack==0) {
    vec *restrict r = (vec *restrict) rr;
    vec *restrict v = (vec *restrict) vv;
    for(int i=i0; i<i1; i++) {
      int k = idx[i];
      vec *restrict ri = &r[nelem*i];
      vec *restrict vi = &v[nelem*k];
      Copy(ri, vi, nelem);
    }
  } else {
#if VLEN > 1
    real *restrict r = (real *restrict) rr;
    vec *restrict v = (vec *restrict) vv;
    int nv2 = nelem*VLEN/2;
#define LOOP(pf) for(int i=i0; i<i1; i++) {		\
      int k = idx[i];				\
      real *restrict ri = &r[nv2*i];		\
      vec *restrict vi = &v[nelem*k];			\
      pf(ri, vi, nelem);				\
    }
    switch(pack) {
    case  1: LOOP(Packp1); break;
    case -1: LOOP(Packm1); break;
#if VLEN > 2
    case  2: LOOP(Packp2); break;
    case -2: LOOP(Packm2); break;
#if VLEN > 4
    case  4: LOOP(Packp4); break;
    case -4: LOOP(Packm4); break;
#if VLEN > 8
    case  8: LOOP(Packp8); break;
    case -8: LOOP(Packm8); break;
#endif
#endif
#endif
    }
#undef LOOP
#endif
  }
}

void
P(X_veq_pack_xX)(real *rr, int pack, int *idx, real *vv,
		 int nelem, int off, int len)
{
  if(pack==0) {
    vec *restrict r = (vec *restrict) rr;
    vec *restrict v = (vec *restrict) vv;
    for(int i=0; i<len; i++) {
      int k = idx[off+i];
      vec *restrict ri = &r[nelem*(off+i)];
      vec *restrict vi = &v[nelem*k];
      Copy(ri, vi, nelem);
    }
  } else {
#if VLEN > 1
    real *restrict r = (real *restrict) rr;
    vec *restrict v = (vec *restrict) vv;
    int nv2 = nelem*VLEN/2;
#define LOOP(pf) for(int i=0; i<len; i++) {		\
      int k = idx[off+i];				\
      real *restrict ri = &r[nv2*(off+i)];		\
      vec *restrict vi = &v[nelem*k];			\
      pf(ri, vi, nelem);				\
    }
    switch(pack) {
    case  1: LOOP(Packp1); break;
    case -1: LOOP(Packm1); break;
#if VLEN > 2
    case  2: LOOP(Packp2); break;
    case -2: LOOP(Packm2); break;
#if VLEN > 4
    case  4: LOOP(Packp4); break;
    case -4: LOOP(Packm4); break;
#if VLEN > 8
    case  8: LOOP(Packp8); break;
    case -8: LOOP(Packm8); break;
#endif
#endif
#endif
    }
#undef LOOP
#endif
  }
}

void
P(X_veq_xXc)(real *rr, int *idx, real *vv, int nelem, int off, int len)
{
  vec *restrict r = (vec *)rr;
  vec *restrict v = (vec *)vv;
  for(int i=0; i<len; i++) {
    int k = idx[off+i];
    if(k>=0) {
      vec *restrict ri = &r[nelem*(off+i)];
      vec *restrict vi = &v[nelem*k];
      for(int j=0; j<nelem; j++) {
	STNR(&(ri[j]), LD(&(vi[j])));
      }
    }
  }
}

void
P(X_veq_xXn)(real *rr, int *idx, real *vv, int nelem, int off, int len)
{
  vec *restrict r = (vec *)rr;
  vec *restrict v = (vec *)vv;
  for(int i=0; i<len; i++) {
    int kk = idx[off+i];
    if(kk<0) {
      int k = -(kk+1);
      vec *restrict ri = &r[nelem*(off+i)];
      vec *restrict vi = &v[nelem*k];
      for(int j=0; j<nelem; j++) {
	STNR(&(ri[j]), LD(&(vi[j])));
      }
    }
  }
}

void
P(X_veq_xXp)(real *rr, int *idx1, int perm, int *idx2, real *vv,
	     int nelem, int off, int len)
{
  vec *restrict r = (vec *)rr;
  vec *restrict v = (vec *)vv;
  for(int i=0; i<len; i++) {
    int k1 = idx1[off+i];
    if(k1>=0) {
      vec *restrict ri = &r[nelem*(off+i)];
      vec *restrict vi = &v[nelem*k1];
      for(int j=0; j<nelem; j++) {
	STNR(&(ri[j]), LD(&(vi[j])));
      }
    } else {
      int k2 = idx2[-(k1+1)];
      vec *restrict ri = &r[nelem*(off+i)];
      vec *restrict vi = &v[nelem*k2];
      for(int j=0; j<nelem; j++) {
	vecr t = LD(&(vi[j]));
	vecr t2 = permN(t, perm);
	STNR(&(ri[j]), t2);
      }
    }
  }
}

void
P(X_eq_xXp)(real *rr, int *idx1, int perm, real *vv,
	    int nelem, int i0, int i1)
{
  vec *restrict r = (vec *)rr;
  vec *restrict v = (vec *)vv;
  for(int i=i0; i<i1; i++) {
    int k1 = idx1[i];
    if(k1>=0) {
      vec *restrict ri = &r[nelem*i];
      vec *restrict vi = &v[nelem*k1];
      for(int j=0; j<nelem; j++) {
	STNR(&(ri[j]), LD(&(vi[j])));
      }
    } else if(k1+2<=0) {
      int k2 = -(k1+2);
      vec *restrict ri = &r[nelem*i];
      vec *restrict vi = &v[nelem*k2];
      for(int j=0; j<nelem; j++) {
	vecr t = LD(&(vi[j]));
	vecr t2 = permN(t, perm);
	STNR(&(ri[j]), t2);
      }
    }
  }
}

void
P(X_veq_xXp2)(real *rr, int *idx1, int perm, real *vv,
	      int nelem, int off, int len)
{
  vec *restrict r = (vec *)rr;
  vec *restrict v = (vec *)vv;
  for(int i=0; i<len; i++) {
    int k1 = idx1[off+i];
    if(k1>=0) {
      vec *restrict ri = &r[nelem*(off+i)];
      vec *restrict vi = &v[nelem*k1];
      for(int j=0; j<nelem; j++) {
	STNR(&(ri[j]), LD(&(vi[j])));
      }
    } else if(k1+2<=0) {
      int k2 = -(k1+2);
      vec *restrict ri = &r[nelem*(off+i)];
      vec *restrict vi = &v[nelem*k2];
      for(int j=0; j<nelem; j++) {
	vecr t = LD(&(vi[j]));
	vecr t2 = permN(t, perm);
	STNR(&(ri[j]), t2);
      }
    }
  }
}

void
P(xX_veq_X)(int *idx1, real *rr, real *vv,
	    int nelem, int off, int len)
{
  vec *restrict r = (vec *)rr;
  vec *restrict v = (vec *)vv;
  for(int i=0; i<len; i++) {
    int k1 = idx1[off+i];
    vec *restrict ri = &r[nelem*k1];
    vec *restrict vi = &v[nelem*(off+i)];
    for(int j=0; j<nelem; j++) {
      STNR(&(ri[j]), LD(&(vi[j])));
    }
  }
}

void
P(xX_veq_yX)(int *idx1, real *rr, int *idx2, real *vv,
	     int nelem, int off, int len)
{
  vec *restrict r = (vec *)rr;
  vec *restrict v = (vec *)vv;
  for(int i=0; i<len; i++) {
    int k1 = idx1[off+i];
    int k2 = idx2[off+i];
    vec *restrict ri = &r[nelem*k1];
    vec *restrict vi = &v[nelem*k2];
    for(int j=0; j<nelem; j++) {
      STNR(&(ri[j]), LD(&(vi[j])));
    }
  }
}

void
P(xX_veq_blend_xXX)(int *idx1, real *rr, int blend, int *idx2,
		    real *vv1, real *vv2, int nelem, int off, int len)
{
#if VLEN > 1
  int nv2 = nelem*VLEN/2;
  vec *restrict r = (vec *restrict) rr;
  vec *restrict v1 = (vec *restrict) vv1;
  real *restrict v2 = (real *restrict) vv2;
#define LOOP(pf) for(int i=0; i<len; i++) {			\
    int k1 = idx1[off+i];					\
    int k2 = idx2[off+i];					\
    vec *restrict ri = &r[nelem*k1];				\
    vec *restrict v1i = &v1[nelem*k2];				\
    real *restrict v2i = &v2[nv2*(off+i)];			\
    pf(ri, v1i, v2i, nelem);					\
  }
  switch(blend) {
    //case  0: LOOP(Blend0); break;
  case  1: LOOP(Blendp1); break;
  case -1: LOOP(Blendm1); break;
#if VLEN > 2
  case  2: LOOP(Blendp2); break;
  case -2: LOOP(Blendm2); break;
#if VLEN > 4
  case  4: LOOP(Blendp4); break;
  case -4: LOOP(Blendm4); break;
#if VLEN > 8
  case  8: LOOP(Blendp8); break;
  case -8: LOOP(Blendm8); break;
#endif
#endif
#endif
  }
#undef LOOP
#endif
}

void
P(xX_veq_blend_xXxX)(int *idx0, real *rr, int blend, int *idx1, real *vv1,
		     int *idx2, real *vv2, int nelem, int off, int len)
{
  if(blend==0) {
    vec *restrict r = (vec *restrict) rr;
    vec *restrict v2 = (vec *restrict) vv2;
    for(int i=0; i<len; i++) {
      int k0 = idx0[off+i];
      int k2 = idx2[off+i];
      vec *restrict ri = &r[nelem*k0];
      vec *restrict v2i = &v2[nelem*k2];
      Copy(ri, v2i, nelem);
    }
  } else {
#if VLEN > 1
    int nv2 = nelem*VLEN/2;
    vec *restrict r = (vec *restrict) rr;
    vec *restrict v1 = (vec *restrict) vv1;
    real *restrict v2 = (real *restrict) vv2;
#define LOOP(pf) for(int i=0; i<len; i++) {			\
      int k0 = idx0[off+i];					\
      int k1 = idx1[off+i];					\
      int k2 = idx2[off+i];					\
      vec *restrict ri = &r[nelem*k0];				\
      vec *restrict v1i = &v1[nelem*k1];				\
      real *restrict v2i = &v2[nv2*k2];			\
      pf(ri, v1i, v2i, nelem);					\
    }
    switch(blend) {
    case  1: LOOP(Blendp1); break;
    case -1: LOOP(Blendm1); break;
#if VLEN > 2
    case  2: LOOP(Blendp2); break;
    case -2: LOOP(Blendm2); break;
#if VLEN > 4
    case  4: LOOP(Blendp4); break;
    case -4: LOOP(Blendm4); break;
#if VLEN > 8
    case  8: LOOP(Blendp8); break;
    case -8: LOOP(Blendm8); break;
#endif
#endif
#endif
    }
#undef LOOP
#endif
  }
}

void
P(xX_eq_blend_xXxX)(int *idx0, real *rr, int blend, int *idx1, real *vv1,
		    int *idx2, real *vv2, int nelem, int i0, int i1)
{
  if(blend==0) {
    vec *restrict r = (vec *restrict) rr;
    vec *restrict v2 = (vec *restrict) vv2;
    for(int i=i0; i<i1; i++) {
      int k0 = idx0[i];
      int k2 = idx2[i];
      vec *restrict ri = &r[nelem*k0];
      vec *restrict v2i = &v2[nelem*k2];
      Copy(ri, v2i, nelem);
    }
  } else {
#if VLEN > 1
    int nv2 = nelem*VLEN/2;
    vec *restrict r = (vec *restrict) rr;
    vec *restrict v1 = (vec *restrict) vv1;
    real *restrict v2 = (real *restrict) vv2;
#define LOOP(pf) for(int i=i0; i<i1; i++) {		\
      int k0 = idx0[i];					\
      int k1 = idx1[i];					\
      int k2 = idx2[i];						\
      vec *restrict ri = &r[nelem*k0];					\
      vec *restrict v1i = &v1[nelem*k1];				\
      real *restrict v2i = &v2[nv2*k2];					\
      pf(ri, v1i, v2i, nelem);						\
    }
    switch(blend) {
    case  1: LOOP(Blendp1); break;
    case -1: LOOP(Blendm1); break;
#if VLEN > 2
    case  2: LOOP(Blendp2); break;
    case -2: LOOP(Blendm2); break;
#if VLEN > 4
    case  4: LOOP(Blendp4); break;
    case -4: LOOP(Blendm4); break;
#if VLEN > 8
    case  8: LOOP(Blendp8); break;
    case -8: LOOP(Blendm8); break;
#endif
#endif
#endif
    }
#undef LOOP
#endif
  }
}

void
P(X_veq_blend_xXxXn)(real *rr, int blend, int *idx1, real *vv1,
		     int *idx2, real *vv2, int nelem, int off, int len)
{
  if(blend==0) {
    vec *restrict r = (vec *restrict) rr;
    //vec *restrict v1 = (vec *restrict) vv1;
    vec *restrict v2 = (vec *restrict) vv2;
    for(int i=0; i<len; i++) {
      int kk2 = idx2[off+i];
      if(kk2<0) {
	int k2 = -(kk2+1);
	//int k1 = idx1[k2];
	vec *restrict ri = &r[nelem*(off+i)];
	//vec *restrict v1i = &v1[nelem*k1];
	vec *restrict v2i = &v2[nelem*k2];
	Copy(ri, v2i, nelem);
      }
    }
  } else {
#if VLEN > 1
    int nv2 = nelem*VLEN/2;
    vec *restrict r = (vec *restrict) rr;
    vec *restrict v1 = (vec *restrict) vv1;
    real *restrict v2 = (real *restrict) vv2;
#define LOOP(pf) for(int i=0; i<len; i++) {			\
      int kk2 = idx2[off+i];					\
      if(kk2<0) {						\
	int k2 = -(kk2+1);					\
	int k1 = idx1[k2];					\
	vec *restrict ri = &r[nelem*(off+i)];			\
	vec *restrict v1i = &v1[nelem*k1];			\
	real *restrict v2i = &v2[nv2*k2];			\
	pf(ri, v1i, v2i, nelem);				\
      }								\
    }
    switch(blend) {
    case  1: LOOP(Blendp1); break;
    case -1: LOOP(Blendm1); break;
#if VLEN > 2
    case  2: LOOP(Blendp2); break;
    case -2: LOOP(Blendm2); break;
#if VLEN > 4
    case  4: LOOP(Blendp4); break;
    case -4: LOOP(Blendm4); break;
#if VLEN > 8
    case  8: LOOP(Blendp8); break;
    case -8: LOOP(Blendm8); break;
#endif
#endif
#endif
    }
#undef LOOP
#endif
  }
}

void
P(V_veq_perm_xV)(real *rr, int perm, int *idx,
		 real *vv, int off, int len)
{
  Vec3 *restrict r = off + (Vec3 *)rr;
  Vec3 *restrict v = (Vec3 *)vv;
#define LOOP(pf) for(int i=0; i<len; i++) {	\
    int k = idx[off+i];				\
    Vec3 *restrict ri = &r[i];			\
    Vec3 *restrict vi = &v[k];			\
    pf(ri, vi);					\
  }
  switch(perm) {
  case 0: LOOP(CopyV3); break;
#if VLEN > 1
  case 1: LOOP(Perm1V3); break;
#endif
#if VLEN > 2
  case 2: LOOP(Perm2V3); break;
#endif
#if VLEN > 4
  case 4: LOOP(Perm4V3); break;
#endif
#if VLEN > 8
  case 8: LOOP(Perm8V3); break;
#endif
  }
#undef LOOP
}

#if 0
void
P(r_eq_sumElem_X)(double *r, real *vv, int nelem, int i0, int i1)
{
  real zero = 0;
  vecr rv = splat(zero);
  vec *restrict v = (vec *)vv;
  int ii0 = nelem*i0;
  int ii1 = nelem*i1;
  for(int i=ii0; i<ii1; i++) {
    vec *restrict vi = &v[i];
    vecr x = LD(vi);
    rv = ADD(x,rv);
  }
  *r = vreduce(rv);
}
void
r_veq_norm2_X(double *r, real *vv, int nelem, int i0, int i1)
{
  real zero = 0;
  vecr rv = splat(zero);
  vec *restrict v = (vec *)vv;
  int ii0 = nelem*i0;
  int ii1 = nelem*i1;
  for(int i=ii0; i<ii1; i++) {
    vec *restrict vi = &v[i];
    vecr x = LD(vi);
    rv = MADD(x,x,rv);
  }
  *r = vreduce(rv);
}
void
r_veq_re_X_dot_X(double *r, real *vv, real *ww, int nelem, int i0, int i1)
{
  real zero = 0;
  vecr rv = splat(zero);
  vec *restrict v = (vec *)vv;
  vec *restrict w = (vec *)ww;
  int ii0 = nelem*i0;
  int ii1 = nelem*i1;
  for(int i=ii0; i<ii1; i++) {
    vec *restrict vi = &v[i];
    vec *restrict wi = &w[i];
    vecr x = LD(vi);
    vecr y = LD(wi);
    rv = MADD(x,y,rv);
  }
  *r = vreduce(rv);
}
#else
void
P(r_eq_sumElem_X)(double *r, real *vv, int nelem, int i0, int i1)
{
  double t = 0;
  int ii0 = nelem*i0*VLEN;
  int ii1 = nelem*i1*VLEN;
  for(int i=ii0; i<ii1; i++) {
    double x = vv[i];
    t += x;
  }
  *r = t;
}
void
r_veq_norm2_X(double *r, real *vv, int nelem, int i0, int i1)
{
  double t = 0;
  int ii0 = nelem*i0*VLEN;
  int ii1 = nelem*i1*VLEN;
  for(int i=ii0; i<ii1; i++) {
    double x = vv[i];
    t += x*x;
  }
  *r = t;
}
void
r_veq_re_X_dot_X(double *r, real *vv, real *ww, int nelem, int i0, int i1)
{
  double t = 0;
  int ii0 = nelem*i0*VLEN;
  int ii1 = nelem*i1*VLEN;
  for(int i=ii0; i<ii1; i++) {
    double x = vv[i];
    double y = ww[i];
    t += x*y;
  }
  *r = t;
}
#endif

double
sumElem_X(real *vv, int nelem, int i0, int i1)
{
  double r;
  P(r_eq_sumElem_X)(&r, vv, nelem, i0, i1);
  return r;
}

double
norm2_X(real *vv, int nelem, int i0, int i1)
{
  double r;
  r_veq_norm2_X(&r, vv, nelem, i0, i1);
  return r;
}

double
re_X_dot_X(real *xx, real *yy, int nelem, int i0, int i1)
{
  double r;
  r_veq_re_X_dot_X(&r, xx, yy, nelem, i0, i1);
  return r;
}

void
X_veq_r(real *rr, real aa, int nelem, int i0, int i1)
{
  vecr av = splat(aa);
  vec *restrict r = (vec *)rr;
  int ii0 = nelem*i0;
  int ii1 = nelem*i1;
  for(int i=ii0; i<ii1; i++) {
    vec *restrict ri = &r[i];
    STNR(ri, av);
  }
}

void
X_veq_r_times_X(real *rr, real aa, real *vv, int nelem, int i0, int i1)
{
  vecr av = splat(aa);
  vec *restrict r = (vec *)rr;
  vec *restrict v = (vec *)vv;
  int ii0 = nelem*i0;
  int ii1 = nelem*i1;
  for(int i=ii0; i<ii1; i++) {
    vec *restrict ri = &r[i];
    vec *restrict vi = &v[i];
    vecr x = LD(vi);
    vecr y = MUL(av, x);
    ST(ri, y);
  }
}

void
X_vpeq_r_times_X(real *rr, real aa, real *vv, int nelem, int i0, int i1)
{
  vecr av = splat(aa);
  vec *restrict r = (vec *)rr;
  vec *restrict v = (vec *)vv;
  int ii0 = nelem*i0;
  int ii1 = nelem*i1;
  for(int i=ii0; i<ii1; i++) {
    vec *restrict ri = &r[i];
    vec *restrict vi = &v[i];
    vecr x = LD(vi);
    vecr s = LD(ri);
    vecr y = MADD(av, x, s);
    ST(ri, y);
  }
}

void
X_veq_X_plus_r_times_X(real *rr, real *vv, real aa, real *ww,
			  int nelem, int i0, int i1)
{
  vecr av = splat(aa);
  vec *restrict r = (vec *)rr;
  vec *restrict v = (vec *)vv;
  vec *restrict w = (vec *)ww;
  int ii0 = nelem*i0;
  int ii1 = nelem*i1;
  for(int i=ii0; i<ii1; i++) {
    vec *restrict ri = &r[i];
    vec *restrict vi = &v[i];
    vec *restrict wi = &w[i];
    vecr x = LD(vi);
    vecr y = LD(wi);
    vecr z = MADD(av,y,x);
    ST(ri, z);
  }
}

void
X_eq_scale_plus_X(real *rr, real aa, real *vv, int nelem, int i0, int i1)
{
  vecr av = splat(aa);
  vec *restrict r = (vec *)rr;
  vec *restrict v = (vec *)vv;
  int ii0 = nelem*i0;
  int ii1 = nelem*i1;
  for(int i=ii0; i<ii1; i++) {
    vec *restrict ri = &r[i];
    vec *restrict vi = &v[i];
    vecr x = LD(ri);
    vecr y = LD(vi);
    vecr z = MADD(av,x,y);
    ST(ri, z);
  }
}

void
X_eq_scale_plus_r_times_X(real *rr, real aa, real bb, real *vv,
			  int nelem, int i0, int i1)
{
  vecr av = splat(aa);
  vecr bv = splat(bb);
  vec *restrict r = (vec *)rr;
  vec *restrict v = (vec *)vv;
  int ii0 = nelem*i0;
  int ii1 = nelem*i1;
  for(int i=ii0; i<ii1; i++) {
    vec *restrict ri = &r[i];
    vec *restrict vi = &v[i];
    vecr x = LD(ri);
    vecr y = LD(vi);
    vecr t = MUL(bv,y);
    vecr z = MADD(av,x,t);
    ST(ri, z);
  }
}

void
Xadj(real *rr, int nc, int i0, int i1)
{
  int nelem = nc*nc;
  cvec *restrict r = (cvec *)rr;
  int ii0 = nelem*i0;
  int ii1 = nelem*i1;
  for(int i=ii0; i<ii1; i+=nelem) {
    cvec *restrict ri = &r[i];
    for(int a=0; a<nc; a++) {
      for(int b=0; b<a; b++) {
	int k1 = a*nc+b;
	int k2 = b*nc+a;
	vecr xr = LD(&ri[k1].r);
	vecr xi = LD(&ri[k1].i);
	vecr yr = LD(&ri[k2].r);
	vecr yi = LD(&ri[k2].i);
	ST(&ri[k2].r, xr);
	ST(&ri[k2].i, NEG(xi));
	ST(&ri[k1].r, yr);
	ST(&ri[k1].i, NEG(yi));
      }
      int k = a*(nc+1);
      vecr xi = LD(&ri[k].i);
      ST(&ri[k].i, NEG(xi));
    }
  }
}

void
X_eq_Xa(real *rr, real *vv, int nr, int nc, int i0, int i1)
{
  int nelem = nr*nc;
  cvec *restrict r = (cvec *)rr;
  cvec *restrict v = (cvec *)vv;
  int ii0 = nelem*i0;
  int ii1 = nelem*i1;
  for(int i=ii0; i<ii1; i+=nelem) {
    cvec *restrict ri = &r[i];
    cvec *restrict vi = &v[i];
    for(int a=0; a<nc; a++) {
      for(int b=0; b<nr; b++) {
	vecr x = LD(&vi[b*nc+a].r);
	vecr y = LD(&vi[b*nc+a].i);
	ST(&ri[a*nr+b].r, x);
	ST(&ri[a*nr+b].i, NEG(y));
      }
    }
  }
}

void
P(V_veq_r_times_V)(real *rr, real aa, real *vv, int off, int len)
{
  if(aa) {
    vecr av = splat(aa);
    Vec3 *restrict r = off + (Vec3 *)rr;
    Vec3 *restrict v = off + (Vec3 *)vv;
    for(int i=0; i<len; i++) {
      Vec3 *restrict ri = &r[i];
      Vec3 *restrict vi = &v[i];
      RtimesV3(ri, av, vi);
    }
  } else { // set to zero
    X_veq_r(rr, aa, 6, off, off+len);
  }
}

void
P(V_vpeq_M_times_xV)(real *rr, real *mm, int *idx,
		     real *vv, int off, int len)
{
  Vec3 *restrict r = off + (Vec3 *)rr;
  Mat3 *restrict m = off + (Mat3 *)mm;
  Vec3 *restrict v = (Vec3 *)vv;
  for(int i=0; i<len; i++) {
    int k = idx[off+i];
    Vec3 *restrict ri = &r[i];
    Mat3 *restrict mi = &m[i];
    Vec3 *restrict vi = &v[k];
    PeqM3AtimesV3o(ri, mi, vi);
  }
}

void
P(V_vpeq_M_times_xVc)(real *rr, real *mm, int *idx,
		      real *vv, int i0, int i1)
{
  Vec3 *restrict r = (Vec3 *)rr;
  Mat3 *restrict m = (Mat3 *)mm;
  Vec3 *restrict v = (Vec3 *)vv;
  for(int i=i0; i<i1; i++) {
    int k = idx[i];
    if(k>=0) {
      Vec3 *restrict ri = &r[i];
      Mat3 *restrict mi = &m[i];
      Vec3 *restrict vi = &v[k];
      PeqM3AtimesV3o(ri, mi, vi);
    }
  }
}

// r[i,j,k] = m[j,l] v[i,l,k]  (i:ldv,j:nc,k:tdv,l:nc)
static inline void
PeqMtimesX(vec *r, vec *m, vec *v, int ldv, int nc, int tdv)
{
  for(int i=0; i<ldv; i++) {
    int i0 = i*nc;
    for(int k=0; k<tdv; k+=2) {
      for(int l=0; l<nc; l++) {
	int il = (i0+l)*tdv;
	vecr vilkr = LD(&v[il+k]);
	vecr vilki = LD(&v[il+k+1]);
	for(int j=0; j<nc; j++) {
	  int ij = (i0+j)*tdv;
	  int jl = 2*(j+l*nc);
	  vecr mjlr = LD(&m[jl]);
	  vecr rijkr = LD(&r[ij+k]);
	  vecr rijki = LD(&r[ij+k+1]);
	  rijkr = MADD(mjlr, vilkr, rijkr);
	  rijki = MADD(mjlr, vilki, rijki);
	  ST(&r[ij+k], rijkr);
	  ST(&r[ij+k+1], rijki);
	}
	for(int j=0; j<nc; j++) {
	  int ij = (i0+j)*tdv;
	  int jl = 2*(j+l*nc);
	  vecr mjli = LD(&m[jl+1]);
	  vecr rijkr = LD(&r[ij+k]);
	  vecr rijki = LD(&r[ij+k+1]);
	  rijkr = NMADD(mjli, vilki, rijkr);
	  rijki = MADD(mjli, vilkr, rijki);
	  ST(&r[ij+k], rijkr);
	  ST(&r[ij+k+1], rijki);
	}
      }
    }
  }
}

// r[i,j,k] = m[j,l] v[i,l,k]*  (i:ldv,j:nc,k:tdv,l:nc)
static inline void
PeqMtimesXc(vec *r, vec *m, vec *v, int ldv, int nc, int tdv)
{
  for(int i=0; i<ldv; i++) {
    int i0 = i*nc;
    for(int k=0; k<tdv; k+=2) {
      for(int l=0; l<nc; l++) {
	int il = (i0+l)*tdv;
	vecr vilkr = LD(&v[il+k]);
	vecr vilki = LD(&v[il+k+1]);
	for(int j=0; j<nc; j++) {
	  int ij = (i0+j)*tdv;
	  int jl = 2*(j+l*nc);
	  vecr mjlr = LD(&m[jl]);
	  vecr rijkr = LD(&r[ij+k]);
	  vecr rijki = LD(&r[ij+k+1]);
	  rijkr = MADD(mjlr, vilkr, rijkr);
	  rijki = NMADD(mjlr, vilki, rijki);
	  ST(&r[ij+k], rijkr);
	  ST(&r[ij+k+1], rijki);
	}
	for(int j=0; j<nc; j++) {
	  int ij = (i0+j)*tdv;
	  int jl = 2*(j+l*nc);
	  vecr mjli = LD(&m[jl+1]);
	  vecr rijkr = LD(&r[ij+k]);
	  vecr rijki = LD(&r[ij+k+1]);
	  rijkr = MADD(mjli, vilki, rijkr);
	  rijki = MADD(mjli, vilkr, rijki);
	  ST(&r[ij+k], rijkr);
	  ST(&r[ij+k+1], rijki);
	}
      }
    }
  }
}

// r[i,j,k] = m[j,l] v[k,l,i]*  (i:ldv,j:nc,k:tdv/2,l:nc)
static inline void
PeqMtimesXa(vec *r, vec *m, vec *v, int ldv, int nc, int tdv)
{
  for(int i=0; i<ldv; i++) {
    int i0 = i*nc;
    for(int k=0; k<tdv; k+=2) {
      for(int l=0; l<nc; l++) {
	int kl = ((k>>1)*nc+l)*ldv*2;
	vecr vilkr = LD(&v[kl+2*i]);
	vecr vilki = LD(&v[kl+2*i+1]);
	for(int j=0; j<nc; j++) {
	  int ij = (i0+j)*tdv;
	  int jl = 2*(j+l*nc);
	  vecr mjlr = LD(&m[jl]);
	  vecr rijkr = LD(&r[ij+k]);
	  vecr rijki = LD(&r[ij+k+1]);
	  rijkr = MADD(mjlr, vilkr, rijkr);
	  rijki = NMADD(mjlr, vilki, rijki);
	  ST(&r[ij+k], rijkr);
	  ST(&r[ij+k+1], rijki);
	}
	for(int j=0; j<nc; j++) {
	  int ij = (i0+j)*tdv;
	  int jl = 2*(j+l*nc);
	  vecr mjli = LD(&m[jl+1]);
	  vecr rijkr = LD(&r[ij+k]);
	  vecr rijki = LD(&r[ij+k+1]);
	  rijkr = MADD(mjli, vilki, rijkr);
	  rijki = MADD(mjli, vilkr, rijki);
	  ST(&r[ij+k], rijkr);
	  ST(&r[ij+k+1], rijki);
	}
      }
    }
  }
}

// r[i,j,k] = m[l,j]* v[i,l,k]  (i:ldv,j:nc,k:tdv,l:nc)
static inline void
PeqMatimesX(vec *r, vec *m, vec *v, int ldv, int nc, int tdv)
{
  for(int i=0; i<ldv; i++) {
    int i0 = i*nc;
    for(int k=0; k<tdv; k+=2) {
      for(int l=0; l<nc; l++) {
	int il = (i0+l)*tdv;
	vecr vilkr = LD(&v[il+k]);
	vecr vilki = LD(&v[il+k+1]);
	for(int j=0; j<nc; j++) {
	  int ij = (i0+j)*tdv;
	  //int jl = 2*(j+l*nc);
	  int jl = 2*(j*nc+l);
	  vecr mjlr = LD(&m[jl]);
	  vecr rijkr = LD(&r[ij+k]);
	  vecr rijki = LD(&r[ij+k+1]);
	  rijkr = MADD(mjlr, vilkr, rijkr);
	  rijki = MADD(mjlr, vilki, rijki);
	  ST(&r[ij+k], rijkr);
	  ST(&r[ij+k+1], rijki);
	}
	for(int j=0; j<nc; j++) {
	  int ij = (i0+j)*tdv;
	  //int jl = 2*(j+l*nc);
	  int jl = 2*(j*nc+l);
	  vecr mjli = LD(&m[jl+1]);
	  vecr rijkr = LD(&r[ij+k]);
	  vecr rijki = LD(&r[ij+k+1]);
	  rijkr = MADD(mjli, vilki, rijkr);
	  rijki = NMADD(mjli, vilkr, rijki);
	  ST(&r[ij+k], rijkr);
	  ST(&r[ij+k+1], rijki);
	}
      }
    }
  }
}

// r[i,j,k] = m[j,l] v[i,l,k]  (i:ldv,j:nc,k:tdv,l:nc)
static inline void
PeqMtimesXp(vec *r, vec *m, vec *v, int ldv, int nc, int tdv, int perm)
{
  for(int i=0; i<ldv; i++) {
    int i0 = i*nc;
    for(int k=0; k<tdv; k+=2) {
      for(int l=0; l<nc; l++) {
	int il = (i0+l)*tdv;
	vecr vilkrp = LD(&v[il+k]);
	vecr vilkip = LD(&v[il+k+1]);
	vecr vilkr = permN(vilkrp,perm);
	vecr vilki = permN(vilkip,perm);
	for(int j=0; j<nc; j++) {
	  int ij = (i0+j)*tdv;
	  int jl = 2*(j+l*nc);
	  vecr mjlr = LD(&m[jl]);
	  vecr rijkr = LD(&r[ij+k]);
	  vecr rijki = LD(&r[ij+k+1]);
	  rijkr = MADD(mjlr, vilkr, rijkr);
	  rijki = MADD(mjlr, vilki, rijki);
	  ST(&r[ij+k], rijkr);
	  ST(&r[ij+k+1], rijki);
	}
	for(int j=0; j<nc; j++) {
	  int ij = (i0+j)*tdv;
	  int jl = 2*(j+l*nc);
	  vecr mjli = LD(&m[jl+1]);
	  vecr rijkr = LD(&r[ij+k]);
	  vecr rijki = LD(&r[ij+k+1]);
	  rijkr = NMADD(mjli, vilki, rijkr);
	  rijki = MADD(mjli, vilkr, rijki);
	  ST(&r[ij+k], rijkr);
	  ST(&r[ij+k+1], rijki);
	}
      }
    }
  }
}

// r[i,j,k] = m[j,l] v[i,l,k]  (i:ldv,j:nc,k:tdv,l:nc)
static inline void
PeqMtimesXb(vec *r, vec *m, vec *v1, real *v2,
	    int ldv, int nc, int tdv, int blend)
{
  int vl2 = VLEN/2;
  for(int i=0; i<ldv; i++) {
    int i0 = i*nc;
    for(int k=0; k<tdv; k+=2) {
      for(int l=0; l<nc; l++) {
	int il = (i0+l)*tdv;
	vecr vilkrb = LD(&v1[il+k]);
	vecr vilkib = LD(&v1[il+k+1]);
	vecr vilkr = blendN(vilkrb,&v2[(il+k)*vl2],blend);
	vecr vilki = blendN(vilkib,&v2[(il+k+1)*vl2],blend);
	for(int j=0; j<nc; j++) {
	  int ij = (i0+j)*tdv;
	  int jl = 2*(j+l*nc);
	  vecr mjlr = LD(&m[jl]);
	  vecr rijkr = LD(&r[ij+k]);
	  vecr rijki = LD(&r[ij+k+1]);
	  rijkr = MADD(mjlr, vilkr, rijkr);
	  rijki = MADD(mjlr, vilki, rijki);
	  ST(&r[ij+k], rijkr);
	  ST(&r[ij+k+1], rijki);
	}
	for(int j=0; j<nc; j++) {
	  int ij = (i0+j)*tdv;
	  int jl = 2*(j+l*nc);
	  vecr mjli = LD(&m[jl+1]);
	  vecr rijkr = LD(&r[ij+k]);
	  vecr rijki = LD(&r[ij+k+1]);
	  rijkr = NMADD(mjli, vilki, rijkr);
	  rijki = MADD(mjli, vilkr, rijki);
	  ST(&r[ij+k], rijkr);
	  ST(&r[ij+k+1], rijki);
	}
      }
    }
  }
}

#if 0
// r[i,j,k] = m[j,l] v[i,l,k]  (i:ldv,j:nc,k:tdv,l:nc)
static inline void
PeqMtimesXp(vec *r, vec *m, vec *v, int ldv, int nc, int tdv, int perm)
{
  for(int i=0; i<ldv; i++) {
    int i0 = i*nc;
    for(int k=0; k<tdv; k+=2) {
      for(int l=0; l<nc; l++) {
	int il = (i0+l)*tdv;
	vecr vilkr = permN(LD(&v[il+k]),perm);
	vecr vilki = permN(LD(&v[il+k+1]),perm);
	for(int j=0; j<nc; j++) {
	  int ij = (i0+j)*tdv;
	  int jl = 2*(j+l*nc);
	  vecr mjlr = LD(&m[jl]);
	  //vecr mjli = LD(&m[jl+1]);
	  vecr rijkr = LD(&r[ij+k]);
	  vecr rijki = LD(&r[ij+k+1]);
	  rijkr = MADD(mjlr, vilkr, rijkr);
	  rijki = MADD(mjlr, vilki, rijki);
	  ST(&r[ij+k], rijkr);
	  ST(&r[ij+k+1], rijki);
	}
	for(int j=0; j<nc; j++) {
	  int ij = (i0+j)*tdv;
	  int jl = 2*(j+l*nc);
	  //vecr mjlr = LD(&m[jl]);
	  vecr mjli = LD(&m[jl+1]);
	  vecr rijkr = LD(&r[ij+k]);
	  vecr rijki = LD(&r[ij+k+1]);
	  rijkr = NMADD(mjli, vilki, rijkr);
	  rijki = MADD(mjli, vilkr, rijki);
	  ST(&r[ij+k], rijkr);
	  ST(&r[ij+k+1], rijki);
	}
      }
    }
  }
}
#endif

// r[a,b,c] = m[b,d] v[a,d,c]  (a:ldv,b:nc,c:tdv)
void
P(X_peq_M_times_X)(real *rr, real *mm, real *vv, int nelem, int nc, int tdv,
		   int i0, int i1)
{
  int ldv = nelem/(nc*tdv);
  int nelemM = 2*nc*nc;
  vec *restrict r = (vec *)rr;
  vec *restrict m = (vec *)mm;
  vec *restrict v = (vec *)vv;
  for(int i=i0; i<i1; i++) {
    vec *restrict ri = &r[i*nelem];
    vec *restrict mi = &m[i*nelemM];
    vec *restrict vi = &v[i*nelem];
    PeqMtimesX(ri, mi, vi, ldv, nc, tdv);
  }
}

// r[a,b,c] = m[b,d] v[a,d,c]*  (a:ldv,b:nc,c:tdv)
void
P(X_peq_M_times_Xc)(real *rr, real *mm, real *vv, int nelem, int nc, int tdv,
		    int i0, int i1)
{
  int ldv = nelem/(nc*tdv);
  int nelemM = 2*nc*nc;
  vec *restrict r = (vec *)rr;
  vec *restrict m = (vec *)mm;
  vec *restrict v = (vec *)vv;
  for(int i=i0; i<i1; i++) {
    vec *restrict ri = &r[i*nelem];
    vec *restrict mi = &m[i*nelemM];
    vec *restrict vi = &v[i*nelem];
    PeqMtimesXc(ri, mi, vi, ldv, nc, tdv);
  }
}

// r[a,b,c] = m[b,d] v[c,d,a]*  (a:ldv,b:nc,c:tdv/2)
void
P(X_peq_M_times_Xa)(real *rr, real *mm, real *vv, int nelem, int nc, int tdv,
		    int i0, int i1)
{
  int ldv = nelem/(nc*tdv);
  int nelemM = 2*nc*nc;
  vec *restrict r = (vec *)rr;
  vec *restrict m = (vec *)mm;
  vec *restrict v = (vec *)vv;
  for(int i=i0; i<i1; i++) {
    vec *restrict ri = &r[i*nelem];
    vec *restrict mi = &m[i*nelemM];
    vec *restrict vi = &v[i*nelem];
    PeqMtimesXa(ri, mi, vi, ldv, nc, tdv);
  }
}

void
P(X_peq_Ma_times_X)(real *rr, real *mm, real *vv, int nelem, int nc, int tdv,
		    int i0, int i1)
{
  int ldv = nelem/(nc*tdv);
  int nelemM = 2*nc*nc;
  vec *restrict r = (vec *)rr;
  vec *restrict m = (vec *)mm;
  vec *restrict v = (vec *)vv;
  for(int i=i0; i<i1; i++) {
    vec *restrict ri = &r[i*nelem];
    vec *restrict mi = &m[i*nelemM];
    vec *restrict vi = &v[i*nelem];
    PeqMatimesX(ri, mi, vi, ldv, nc, tdv);
  }
}

// r[a,b,c] = m[b,d] v[a,d,c]  (a:ldv,b:nc,c:tdv)
void
P(X_peq_M_times_xXp)(real *rr, real *mm, int *idx, int perm,
		     real *vv, int nelem, int nc, int tdv, int i0, int i1)
{
  int ldv = nelem/(nc*tdv);
  int nelemM = 2*nc*nc;
  vec *restrict r = (vec *)rr;
  vec *restrict m = (vec *)mm;
  vec *restrict v = (vec *)vv;
  for(int i=i0; i<i1; i++) {
    int k1 = idx[i];
    if(k1>=0) {
      vec *restrict ri = &r[i*nelem];
      vec *restrict mi = &m[i*nelemM];
      vec *restrict vi = &v[nelem*k1];
      PeqMtimesX(ri, mi, vi, ldv, nc, tdv);
    } else if(k1+2<=0) {
      int k2 = -(k1+2);
      vec *restrict ri = &r[i*nelem];
      vec *restrict mi = &m[i*nelemM];
      vec *restrict vi = &v[nelem*k2];
      PeqMtimesXp(ri, mi, vi, ldv, nc, tdv, perm);
    }
  }
}

// r[a,b,c] = m[b,d] v[a,d,c]  (a:ldv,b:nc,c:tdv)
void
P(xX_peq_xM_times_yXzXb)(int *idx0, real *rr, real *mm, int blend, int *idx1,
			 real *vv1, int *idx2, real *vv2, int nelem, int nc,
			 int tdv, int i0, int i1)
{
  int ldv = nelem/(nc*tdv);
  int nelemM = 2*nc*nc;
  if(blend==0) {
    vec *restrict r = (vec *)rr;
    vec *restrict m = (vec *)mm;
    vec *restrict v2 = (vec *)vv2;
    for(int i=i0; i<i1; i++) {
      int k0 = idx0[i];
      int k2 = idx2[i];
      vec *restrict ri = &r[k0*nelem];
      vec *restrict mi = &m[k0*nelemM];
      vec *restrict v2i = &v2[k2*nelem];
      PeqMtimesX(ri, mi, v2i, ldv, nc, tdv);
    }
  } else {
    int nv2 = nelem*(VLEN/2);
    vec *restrict r = (vec *)rr;
    vec *restrict m = (vec *)mm;
    vec *restrict v1 = (vec *)vv1;
    real *restrict v2 = (real *)vv2;
    for(int i=i0; i<i1; i++) {
      int k0 = idx0[i];
      int k1 = idx1[i];
      int k2 = idx2[i];
      vec *restrict ri = &r[k0*nelem];
      vec *restrict mi = &m[k0*nelemM];
      vec *restrict v1i = &v1[k1*nelem];
      real *restrict v2i = &v2[k2*nv2];
      PeqMtimesXb(ri, mi, v1i, v2i, ldv, nc, tdv, blend);
    }
  }
}

void
P(V_vpeq_M_times_xVp)(real *rr, real *mm, int *idx,
		      int perm, real *vv, int i0, int i1)
{
  Vec3 *restrict r = (Vec3 *)rr;
  Mat3 *restrict m = (Mat3 *)mm;
  Vec3 *restrict v = (Vec3 *)vv;
  for(int i=i0; i<i1; i++) {
    int k1 = idx[i];
    if(k1>=0) {
      Vec3 *restrict ri = &r[i];
      Mat3 *restrict mi = &m[i];
      Vec3 *restrict vi = &v[k1];
      PeqM3AtimesV3o(ri, mi, vi);
    } else if(k1+2<=0) {
      int k2 = -(k1+2);
      Vec3 *restrict ri = &r[i];
      Mat3 *restrict mi = &m[i];
      Vec3 *restrict vi = &v[k2];
      PeqM3AtimesV3op(ri, mi, vi, perm);
    }
  }
}

void
P(V_vmeq_xMap_times_xVp)(real *rr, real *mm, int *idx, int perm, real *vv,
			 int i0, int i1)
{
  Vec3 *restrict r = (Vec3 *)rr;
  Mat3 *restrict m = (Mat3 *)mm;
  Vec3 *restrict v = (Vec3 *)vv;
  for(int i=i0; i<i1; i++) {
    int k1 = idx[i];
    if(k1>=0) {
      Vec3 *restrict ri = &r[i];
      Mat3 *restrict mi = &m[k1];
      Vec3 *restrict vi = &v[k1];
      //PeqM3AtimesV3o(ri, mi, vi);
      MeqM3AtimesV3(ri, mi, vi);
    } else if(k1+2<=0) {
      int k2 = -(k1+2);
      Vec3 *restrict ri = &r[i];
      Mat3 *restrict mi = &m[k2];
      Vec3 *restrict vi = &v[k2];
      //PeqM3AtimesV3o(ri, mi, vi);
      //MeqM3AtimesV3(ri, mi, vi);
      MeqM3AtimesV3p(ri, mi, vi, perm);
    }
  }
}

void
P(V_vpeq_M_times_xVn)(real *rr, real *mm, int *idx, real *vv, int i0, int i1)
{
  Vec3 *restrict r = (Vec3 *)rr;
  Mat3 *restrict m = (Mat3 *)mm;
  Vec3 *restrict v = (Vec3 *)vv;
  for(int i=i0; i<i1; i++) {
    int kk = idx[i];
    if(kk<0) {
      int k = -(kk+1);
      Vec3 *restrict ri = &r[i];
      Mat3 *restrict mi = &m[i];
      Vec3 *restrict vi = &v[k];
      PeqM3AtimesV3o(ri, mi, vi);
    }
  }
}

#if 0
void
V_vpeq_M_times_xVb(real *rr, real *mm, int blend, int *idx1, real *vv1,
		   int *idx2, real *vv2, int i0, int i1)
{
  if(blend==0) {
    Vec3 *restrict r = (Vec3 *)rr;
    Mat3 *restrict m = (Mat3 *)mm;
    Vec3 *restrict v2 = (Vec3 *)vv2;
    for(int i=i0; i<i1; i++) {
      int kk2 = idx2[i];
      if(kk2<0) {
	int k2 = -(kk2+1);
	Vec3 *restrict ri = &r[i];
	Mat3 *restrict mi = &m[i];
	Vec3 *restrict v2i = &v2[k2];
	PeqM3AtimesV3o(ri, mi, v2i);
      }
    }
  } if(blend==1) {
    int nv2 = 3*VLEN;
    Vec3 *restrict r = (Vec3 *)rr;
    Mat3 *restrict m = (Mat3 *)mm;
    Vec3 *restrict v1 = (Vec3 *)vv1;
    real *restrict v2 = (real *)vv2;
    for(int i=i0; i<i1; i++) {
      int kk2 = idx2[i];
      if(kk2<0) {
	int k2 = -(kk2+1);
	int k1 = idx1[k2];
	Vec3 *restrict ri = &r[i];
	Mat3 *restrict mi = &m[i];
	Vec3 *restrict v1i = &v1[k1];
	real *restrict v2i = &v2[nv2*k2];
	PeqM3AtimesV3obp1(ri, mi, v1i, v2i);
      }
    }
  } else {
    int nv2 = 3*VLEN;
    Vec3 *restrict r = (Vec3 *)rr;
    Mat3 *restrict m = (Mat3 *)mm;
    Vec3 *restrict v1 = (Vec3 *)vv1;
    real *restrict v2 = (real *)vv2;
    for(int i=i0; i<i1; i++) {
      int kk2 = idx2[i];
      if(kk2<0) {
	int k2 = -(kk2+1);
	int k1 = idx1[k2];
	Vec3 *restrict ri = &r[i];
	Mat3 *restrict mi = &m[i];
	Vec3 *restrict v1i = &v1[k1];
	real *restrict v2i = &v2[nv2*k2];
	PeqM3AtimesV3ob(ri, mi, v1i, v2i, blend);
      }
    }
  }
}

void
xV_vpeq_xM_times_yVzVb(int *idx0, real *rr, real *mm, int blend, int *idx1,
		       real *vv1, int *idx2, real *vv2, int i0, int i1)
{
  if(blend==0) {
    Vec3 *restrict r = (Vec3 *)rr;
    Mat3 *restrict m = (Mat3 *)mm;
    Vec3 *restrict v2 = (Vec3 *)vv2;
    for(int i=i0; i<i1; i++) {
      int k0 = idx0[i];
      int k2 = idx2[i];
      Vec3 *restrict ri = &r[k0];
      Mat3 *restrict mi = &m[k0];
      Vec3 *restrict v2i = &v2[k2];
      PeqM3AtimesV3o(ri, mi, v2i);
    }
  } else {
    int nv2 = 3*VLEN;
    Vec3 *restrict r = (Vec3 *)rr;
    Mat3 *restrict m = (Mat3 *)mm;
    Vec3 *restrict v1 = (Vec3 *)vv1;
    real *restrict v2 = (real *)vv2;
    for(int i=i0; i<i1; i++) {
      int k0 = idx0[i];
      int k1 = idx1[i];
      int k2 = idx2[i];
      Vec3 *restrict ri = &r[k0];
      Mat3 *restrict mi = &m[k0];
      Vec3 *restrict v1i = &v1[k1];
      real *restrict v2i = &v2[nv2*k2];
      PeqM3AtimesV3ob(ri, mi, v1i, v2i, blend);
    }
  }
}
#endif

void
P(xV_vpeq_M_times_V)(int *idx, real *rr, real *mm, real *vv, int off, int len)
{
  Vec3 *restrict r = (Vec3 *)rr;
  Mat3 *restrict m = off + (Mat3 *)mm;
  Vec3 *restrict v = off + (Vec3 *)vv;
  for(int i=0; i<len; i++) {
    int k = idx[off+i];
    Vec3 *restrict ri = &r[k];
    Mat3 *restrict mi = &m[i];
    Vec3 *restrict vi = &v[i];
    PeqM3AtimesV3o(ri, mi, vi);
  }
}

void
P(xV_vpeq_xM_times_V)(int *idx, real *rr, real *mm, real *vv, int off, int len)
{
  Vec3 *restrict r = (Vec3 *)rr;
  Mat3 *restrict m = (Mat3 *)mm;
  Vec3 *restrict v = off + (Vec3 *)vv;
  for(int i=0; i<len; i++) {
    int k = idx[off+i];
    Vec3 *restrict ri = &r[k];
    Mat3 *restrict mi = &m[k];
    Vec3 *restrict vi = &v[i];
    PeqM3AtimesV3o(ri, mi, vi);
  }
}

void
P(xV_vpeq_xM_times_yV)(int *idx1, real *rr, real *mm, int *idx2,
		       real *vv, int off, int len)
{
  Vec3 *restrict r = (Vec3 *)rr;
  Mat3 *restrict m = (Mat3 *)mm;
  Vec3 *restrict v = (Vec3 *)vv;
  for(int i=0; i<len; i++) {
    int k1 = idx1[off+i];
    int k2 = idx2[off+i];
    Vec3 *restrict ri = &r[k1];
    Mat3 *restrict mi = &m[k1];
    Vec3 *restrict vi = &v[k2];
    PeqM3AtimesV3o(ri, mi, vi);
  }
}

void
P(xV_vpeq_yM_times_yV)(int *idx1, real *rr, int *idx2,
		       real *mm, real *vv, int begin, int end)
{
  Vec3 *restrict r = (Vec3 *)rr;
  Mat3 *restrict m = (Mat3 *)mm;
  Vec3 *restrict v = (Vec3 *)vv;
  for(int i=begin; i<end; i++) {
    int k1 = idx1[i];
    int k2 = idx2[i];
    Vec3 *restrict ri = &r[k1];
    Mat3 *restrict mi = &m[k2];
    Vec3 *restrict vi = &v[k2];
    PeqM3AtimesV3o(ri, mi, vi);
  }
}

static inline void
cvscale(vec *v, vecr s, int n)
{
  for(int i=0; i<2*n; i++) {
    vecr vi = LD(&v[i]);
    vi = MUL(s,vi);
    ST(&v[i], vi);
  }
}

static inline vecr
cvnorm2(vec *v, int n)
{
  vecr v0 = LD(&v[0]);
  vecr s = MUL(v0,v0);
  for(int i=1; i<2*n; i++) {
    vecr vi = LD(&v[i]);
    s = MADD(vi,vi,s);
  }
  return s;
}

static inline cvecr
cvdot(vec *x, vec *y, int n)
{
  cvecr s;
  vecr x0 = LD(&x[0]);
  vecr x1 = LD(&x[1]);
  vecr y0 = LD(&y[0]);
  vecr y1 = LD(&y[1]);
  s.r = MUL(x0,y0);
  s.i = MUL(x0,y1);
  s.r = MADD(x1,y1,s.r);
  s.i = NMADD(x1,y0,s.i);
  for(int i=2; i<2*n; i+=2) {
    vecr xi0 = LD(&x[i]);
    vecr xi1 = LD(&x[i+1]);
    vecr yi0 = LD(&y[i]);
    vecr yi1 = LD(&y[i+1]);
    s.r = MADD(xi0,yi0,s.r);
    s.i = MADD(xi0,yi1,s.i);
    s.r = MADD(xi1,yi1,s.r);
    s.i = NMADD(xi1,yi0,s.i);
  }
  return s;
}

static inline void
cv_meq_c_times_v(vec *x, cvecr c, vec *y, int n)
{
  for(int i=0; i<2*n; i+=2) {
    vecr xi0 = LD(&x[i]);
    vecr xi1 = LD(&x[i+1]);
    vecr yi0 = LD(&y[i]);
    vecr yi1 = LD(&y[i+1]);
    xi0 = NMADD(c.r, yi0, xi0);
    xi1 = NMADD(c.r, yi1, xi1);
    xi0 =  MADD(c.i, yi1, xi0);
    xi1 = NMADD(c.i, yi0, xi1);
    ST(&x[i], xi0);
    ST(&x[i+1], xi1);
  }
}

static inline void
makeUsite(vec *m, int nc)
{
  for(int i=0; i<nc; i++) {
    vec *mi = &m[i*2*nc];
    vecr n = cvnorm2(mi, nc);
    vecr ni = INVSQRT(n);
    cvscale(mi, ni, nc);
    for(int j=0; j<i; j++) {
      vec *mj = &m[j*2*nc];
      cvecr d = cvdot(mj, mi, nc);
      cv_meq_c_times_v(mi, d, mj, nc);
      //d = cvdot(mi, mj, nc);
      //printf("%g\n", *(double*)&d.r);
      n = cvnorm2(mi, nc);
      ni = INVSQRT(n);
      cvscale(mi, ni, nc);
    }
  }
}

void
P(makeU)(real *mm, int nc, int i0, int i1)
{
  int nelem = 2*nc*nc;
  vec *restrict m = (vec *)mm;
  for(int i=i0; i<i1; i++) {
    vec *restrict mi = &m[i*nelem];
    makeUsite(mi, nc);
  }
}

#if 0
void
XfmaXX(real *xx, real *yy, real*zz, int *c, int n, int nelem, int i0, int i1)
{
  const int nbits = 10;
  const int mask = (1<<nbits)-1;
  const int sign = 1<<(3*nbits);
  vec *restrict x = (vec *)xx;
  vec *restrict y = (vec *)yy;
  vec *restrict z = (vec *)zz;
  for(int i=i0; i<i1; i++) {
    vec *restrict xi = &x[i*nelem];
    vec *restrict yi = &y[i*nelem];
    vec *restrict zi = &z[i*nelem];
    for(int j=0; j<n; j++) {
      int cj = c[j];
      int kz = cj & mask;
      int ky = (cj>>nbits) & mask;
      int kx = (jc>>(2*nbits)) & mask;
      vecr zk = LD(&zi[kz]);
      vecr yk = LD(&yi[ky]);
      vecr xk = LD(&xi[kx]);
      if(cj&sign) {
	NMADD(yk,zk,xk);
      } else {
	MADD(yk,zk,xk);
      }
      ST(&xi[kx], xk);
    }
  }
}
#endif
