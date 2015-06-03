#ifdef C1
#define VLEN 1
#endif

#ifdef C2
#define VLEN 2
#define perm1xD(r,s) (r)[0] = (s)[1]; (r)[1] = (s)[0]
//#define perm1D(r,s) (r) = __builtin_shuffle(s, (int __attribute__((vector_size(8)))){1,0})
#define packp1xD(r,s) (r)[0] = (s)[1]
#define packm1xD(r,s) (r)[0] = (s)[0]
#define blendp1xD(r,s,t) (r)[0] = (t)[0]; (r)[1] = (s)[0]
#define blendm1xD(r,s,t) (r)[0] = (s)[1]; (r)[1] = (t)[0]
#define perm1D(r,s) perm1xD((double*)&(r),(double*)&(s))
#define packp1D(r,s) packp1xD(r,(double*)&(s))
#define packm1D(r,s) packm1xD(r,(double*)&(s))
#define blendp1D(r,s,t) blendp1xD((double*)&(r),(double*)&(s),t)
#define blendm1D(r,s,t) blendm1xD((double*)&(r),(double*)&(s),t)
#endif

#ifdef C4
#define VLEN 4
#define perm1xD(r,s) (r)[0]=(s)[1];(r)[1]=(s)[0];(r)[2]=(s)[3];(r)[3]=(s)[2]
#define perm2xD(r,s) (r)[0]=(s)[2];(r)[2]=(s)[0];(r)[1]=(s)[3];(r)[3]=(s)[1]
//#define perm1D(r,s) (r) = __builtin_shuffle(s, (int __attribute__((vector_size(16)))){1,0,3,2})
//#define perm2(r,s) (r) = __builtin_shuffle(s, (int __attribute__((vector_size(16)))){2,3,0,1})
#define packp1xD(r,s) (r)[0] = (s)[1]; (r)[1] = (s)[3]
#define packm1xD(r,s) (r)[0] = (s)[0]; (r)[1] = (s)[2]
#define blendp1xD(r,s,t) (r)[0]=(t)[0];(r)[1]=(s)[0];(r)[2]=(t)[1];(r)[3]=(s)[2]
#define blendm1xD(r,s,t) (r)[0]=(s)[1];(r)[1]=(t)[0];(r)[2]=(s)[3];(r)[3]=(t)[1]
#define packp2xD(r,s) (r)[0] = (s)[2]; (r)[1] = (s)[3]
#define packm2xD(r,s) (r)[0] = (s)[0]; (r)[1] = (s)[1]
#define blendp2xD(r,s,t) (r)[0]=(t)[0];(r)[1]=(t)[1];(r)[2]=(s)[0];(r)[3]=(s)[1]
#define blendm2xD(r,s,t) (r)[0]=(s)[2];(r)[1]=(s)[3];(r)[2]=(t)[0];(r)[3]=(t)[1]
#define perm1D(r,s) perm1xD((double*)&(r),(double*)&(s))
#define packp1D(r,s) packp1xD(r,(double*)&(s))
#define packm1D(r,s) packm1xD(r,(double*)&(s))
#define blendp1D(r,s,t) blendp1xD((double*)&(r),(double*)&(s),t)
#define blendm1D(r,s,t) blendm1xD((double*)&(r),(double*)&(s),t)
#define perm2D(r,s) perm2xD((double*)&(r),(double*)&(s))
#define packp2D(r,s) packp2xD(r,(double*)&(s))
#define packm2D(r,s) packm2xD(r,(double*)&(s))
#define blendp2D(r,s,t) blendp2xD((double*)&(r),(double*)&(s),t)
#define blendm2D(r,s,t) blendm2xD((double*)&(r),(double*)&(s),t)
#endif

#ifdef SSE
//#include <immintrin.h>
#define VLEN 4
#define perm1D(r,s) (r) = __builtin_shuffle(s, (long __attribute__((vector_size(32)))){1,0,3,2})
#define perm2D(r,s) (r) = __builtin_shuffle(s, (long __attribute__((vector_size(32)))){2,3,0,1})
#define packp1xD(r,s) (r)[0] = (s)[1]; (r)[1] = (s)[3]
#define packm1xD(r,s) (r)[0] = (s)[0]; (r)[1] = (s)[2]
#define blendp1xD(r,s,t) (r)[0]=(t)[0];(r)[1]=(s)[0];(r)[2]=(t)[1];(r)[3]=(s)[2]
#define blendm1xD(r,s,t) (r)[0]=(s)[1];(r)[1]=(t)[0];(r)[2]=(s)[3];(r)[3]=(t)[1]
#define packp2xD(r,s) (r)[0] = (s)[2]; (r)[1] = (s)[3]
#define packm2xD(r,s) (r)[0] = (s)[0]; (r)[1] = (s)[1]
#define blendp2xD(r,s,t) (r)[0]=(t)[0];(r)[1]=(t)[1];(r)[2]=(s)[0];(r)[3]=(s)[1]
#define blendm2xD(r,s,t) (r)[0]=(s)[2];(r)[1]=(s)[3];(r)[2]=(t)[0];(r)[3]=(t)[1]
#define packp1D(r,s) packp1xD(r,(double*)&(s))
#define packm1D(r,s) packm1xD(r,(double*)&(s))
#define blendp1D(r,s,t) blendp1xD((double*)&(r),(double*)&(s),t)
#define blendm1D(r,s,t) blendm1xD((double*)&(r),(double*)&(s),t)
#define packp2D(r,s) packp2xD(r,(double*)&(s))
#define packm2D(r,s) packm2xD(r,(double*)&(s))
#define blendp2D(r,s,t) blendp2xD((double*)&(r),(double*)&(s),t)
#define blendm2D(r,s,t) blendm2xD((double*)&(r),(double*)&(s),t)
#endif // SSE

#ifdef QPX
#define VLEN 4
typedef __vector4double vdr;
#define perm1L(r,s) do { vdr _s = vec_ld(0,(s).v); vdr _r; vdr _perm = vec_gpci(01032); _r = vec_perm(_s,_s,_perm); vec_st(_r,0,(r).v); } while(0)
#define perm2L(r,s) do { vdr _s = vec_ld(0,(s).v); vdr _r; vdr _perm = vec_gpci(02301); _r = vec_perm(_s,_s,_perm); vec_st(_r,0,(r).v); } while(0)
#define perm4L(r,s) vec_st(vec_ld(0,(s).v),0,(r).v)
#define perm8L(r,s) vec_st(vec_ld(0,(s).v),0,(r).v)
#define perm1D(r,s) (r) = vec_perm(s,s,vec_gpci(01032))
#define perm2D(r,s) (r) = vec_perm(s,s,vec_gpci(02301))
#define perm4D(r,s) (r) = (s)
#define perm8D(r,s) (r) = (s)
//#define packp1D(r,s) (r)[0] = (s)[1]; (r)[1] = (s)[3]
//#define packm1D(r,s) (r)[0] = (s)[0]; (r)[1] = (s)[2]
#define packp1D(r,s) do { vdr _t=vec_perm(s,s,vec_gpci(01302)); vec_st2(_t,0,r); } while(0)
#define packm1D(r,s) do { vdr _t=vec_perm(s,s,vec_gpci(00213)); vec_st2(_t,0,r); } while(0)
//#define blendp1(r,s,t) (r)[0]=(t)[0];(r)[1]=(s)[0];(r)[2]=(t)[1];(r)[3]=(s)[2]
//#define blendm1(r,s,t) (r)[0]=(s)[1];(r)[1]=(t)[0];(r)[2]=(s)[3];(r)[3]=(t)[1]
#define blendp1D(r,s,t) do { vdr _t=vec_ld2(0,t); (r) = vec_perm(s,_t,vec_gpci(04052)); } while(0)
#define blendm1D(r,s,t) do { vdr _t=vec_ld2(0,t); (r) = vec_perm(s,_t,vec_gpci(01435)); } while(0)
//#define packp2(r,s) (r)[0] = (s)[2]; (r)[1] = (s)[3]
//#define packm2(r,s) (r)[0] = (s)[0]; (r)[1] = (s)[1]
#define packp2D(r,s) do { vdr _t=vec_perm(s,s,vec_gpci(02301)); vec_st2(_t,0,r); } while(0)
#define packm2D(r,s) do { vec_st2(s,0,r); } while(0)
//#define blendp2(r,s,t) (r)[0]=(t)[0];(r)[1]=(t)[1];(r)[2]=(s)[0];(r)[3]=(s)[1]
//#define blendm2(r,s,t) (r)[0]=(s)[2];(r)[1]=(s)[3];(r)[2]=(t)[0];(r)[3]=(t)[1]
#define blendp2D(r,s,t) do { vdr _t=vec_ld2(0,t); (r) = vec_perm(s,_t,vec_gpci(04501)); } while(0)
#define blendm2D(r,s,t) do { vdr _t=vec_ld2(0,t); (r) = vec_perm(s,_t,vec_gpci(02345)); } while(0)
#define vload(x) vec_ld(0,(double*)(x))
#define vstore(x,y) vec_st(y,0,(x)->v)
#define NEG(x) vec_neg(x)
#define MUL(x,y) vec_mul(x,y)
#define MADD(x,y,z) vec_madd(x,y,z)
#define NMADD(x,y,z) vec_nmsub(x,y,z)
#endif

#ifdef AVX
#include <immintrin.h>
#define VLEN 8
#ifdef INTRINSICS
#define splat(x) _mm256_broadcast_ss(&(x))
#define perm1D(r,s) (r) = _mm256_shuffle_ps(s, s, BASE4(2,3,0,1))
#define perm2D(r,s) (r) = _mm256_shuffle_ps(s, s, BASE4(1,0,3,2))
//#define perm3(r,s) (r) = _mm256_shuffle_ps(s, s, BASE4(0,1,2,3))
#define perm4D(r,s) (r) = _mm256_permute2f128_ps(s, s, 1)
#else
#define perm1D(r,s) (r) = __builtin_shuffle(s,(long __attribute__((vector_size(64)))){1,0,3,2,5,4,7,6})
#define perm2D(r,s) (r) = __builtin_shuffle(s,(long __attribute__((vector_size(64)))){2,3,0,1,6,7,4,5})
#define perm4D(r,s) (r) = __builtin_shuffle(s,(long __attribute__((vector_size(64)))){4,5,6,7,0,1,2,3})
#endif
#define packp1xD(r,s) (r)[0]=(s)[1];(r)[1]=(s)[3];(r)[2]=(s)[5];(r)[3]=(s)[7]
#define packm1xD(r,s) (r)[0]=(s)[0];(r)[1]=(s)[2];(r)[2]=(s)[4];(r)[3]=(s)[6]
#define blendp1xD(r,s,t) (r)[0]=(t)[0];(r)[1]=(s)[0];(r)[2]=(t)[1];(r)[3]=(s)[2];(r)[4]=(t)[2];(r)[5]=(s)[4];(r)[6]=(t)[3];(r)[7]=(s)[6]
#define blendm1xD(r,s,t) (r)[0]=(s)[1];(r)[1]=(t)[0];(r)[2]=(s)[3];(r)[3]=(t)[1];(r)[4]=(s)[5];(r)[5]=(t)[2];(r)[6]=(s)[7];(r)[7]=(t)[3]
#define packp2xD(r,s) (r)[0]=(s)[2];(r)[1]=(s)[3];(r)[2]=(s)[6];(r)[3]=(s)[7]
#define packm2xD(r,s) (r)[0]=(s)[0];(r)[1]=(s)[1];(r)[2]=(s)[4];(r)[3]=(s)[5]
#define blendp2xD(r,s,t) (r)[0]=(t)[0];(r)[1]=(t)[1];(r)[2]=(s)[0];(r)[3]=(s)[1];(r)[4]=(t)[2];(r)[5]=(t)[3];(r)[6]=(s)[4];(r)[7]=(s)[5]
#define blendm2xD(r,s,t) (r)[0]=(s)[2];(r)[1]=(s)[3];(r)[2]=(t)[0];(r)[3]=(t)[1];(r)[4]=(s)[6];(r)[5]=(s)[7];(r)[6]=(t)[2];(r)[7]=(t)[3]
#define packp4xD(r,s) (r)[0]=(s)[4];(r)[1]=(s)[5];(r)[2]=(s)[6];(r)[3]=(s)[7]
#define packm4xD(r,s) (r)[0]=(s)[0];(r)[1]=(s)[1];(r)[2]=(s)[2];(r)[3]=(s)[3]
#define blendp4xD(r,s,t) (r)[0]=(t)[0];(r)[1]=(t)[1];(r)[2]=(t)[2];(r)[3]=(t)[3];(r)[4]=(s)[0];(r)[5]=(s)[1];(r)[6]=(s)[2];(r)[7]=(s)[3]
#define blendm4xD(r,s,t) (r)[0]=(s)[4];(r)[1]=(s)[5];(r)[2]=(s)[6];(r)[3]=(s)[7];(r)[4]=(t)[0];(r)[5]=(t)[1];(r)[6]=(t)[2];(r)[7]=(t)[3]
#define packp1D(r,s) packp1xD(r,(double*)&(s))
#define packm1D(r,s) packm1xD(r,(double*)&(s))
#define blendp1D(r,s,t) blendp1xD((double*)&(r),(double*)&(s),t)
#define blendm1D(r,s,t) blendm1xD((double*)&(r),(double*)&(s),t)
#define packp2D(r,s) packp2xD(r,(double*)&(s))
#define packm2D(r,s) packm2xD(r,(double*)&(s))
#define blendp2D(r,s,t) blendp2xD((double*)&(r),(double*)&(s),t)
#define blendm2D(r,s,t) blendm2xD((double*)&(r),(double*)&(s),t)
#define packp4D(r,s) packp4xD(r,(double*)&(s))
#define packm4D(r,s) packm4xD(r,(double*)&(s))
#define blendp4D(r,s,t) blendp4xD((double*)&(r),(double*)&(s),t)
#define blendm4D(r,s,t) blendm4xD((double*)&(r),(double*)&(s),t)
#endif // AVX

#ifdef MIC
#include <immintrin.h>
#define VLEN 16
#define splat(x) _mm512_extload_ps(&(x), _MM_UPCONV_PS_NONE, _MM_BROADCAST_1X16, 0)
#define perm1D(r,s) (r) = _mm512_swizzle_ps(s, _MM_SWIZ_REG_CDAB)
#define perm2(r,s) (r) = _mm512_swizzle_ps(s, _MM_SWIZ_REG_BADC)
#define perm4(r,s) (r) = _mm512_permute4f128_ps(s, BASE4(2,3,0,1))
#define perm8(r,s) (r) = _mm512_permute4f128_ps(s, BASE4(1,0,3,2))
#endif // MIC

#define BASE4(a,b,c,d) (64*(a)+16*(b)+4*(c)+(d))

#ifdef C1
typedef double vd;
typedef vd vdr;
#else
#ifdef QPX
//typedef double vd[4];
typedef struct { double v[4]; } vd;
#else
typedef double vd __attribute__ ((vector_size(8*VLEN))) __attribute__ ((aligned(8*VLEN)));
typedef vd vdr;
#endif
#endif

#define CAT(a,b) CAT1(a,b)
#define CAT1(a,b) a ## b
#define REP1(x) x
#define REP2(x) x,x
#define REP4(x) x,x,x,x
#define REP8(x) x,x,x,x,x,x,x,x
#define REP16(x) REP8(x),REP8(x)
#define VINDEX(x,n) (((double*)&(x))[n])
#define VRED1(x) VINDEX(x,0)
#define VRED2(x) (VRED1(x)+VRED1(VINDEX(x,1)))
#define VRED4(x) (VRED2(x)+VRED2(VINDEX(x,2)))
#define VRED8(x) (VRED4(x)+VRED4(VINDEX(x,4)))

#ifndef splat
//#define splat(x) (vd)(x)
#define splat(x) (vdr){CAT(REP,VLEN)(x)}
#endif
#ifndef vreduce
#define vreduce(r,x) (r) = CAT(VRED,VLEN)(x)
#endif
#ifndef perm1D
#define perm1D(r,s) (r) = (s)
#endif
#ifndef perm2D
#define perm2D(r,s) (r) = (s)
#endif
#ifndef perm4D
#define perm4D(r,s) (r) = (s)
#endif
#ifndef perm8D
#define perm8D(r,s) (r) = (s)
#endif
#ifndef packp1D
#define packp1D(r,s)
#define packm1D(r,s)
#define blendp1D(r,s,t)
#define blendm1D(r,s,t)
#endif
#ifndef packp2D
#define packp2D(r,s)
#define packm2D(r,s)
#define blendp2D(r,s,t)
#define blendm2D(r,s,t)
#endif
#ifndef packp4D
#define packp4D(r,s)
#define packm4D(r,s)
#define blendp4D(r,s,t)
#define blendm4D(r,s,t)
#endif
#ifndef packp8D
#define packp8D(r,s)
#define packm8D(r,s)
#define blendp8D(r,s,t)
#define blendm8D(r,s,t)
#endif

#ifndef perm1L
#define perm1L(r,s) perm1D(r,s)
#endif
#ifndef perm2L
#define perm2L(r,s) perm2D(r,s)
#endif
#ifndef perm4L
#define perm4L(r,s) perm4D(r,s)
#endif
#ifndef perm8L
#define perm8L(r,s) perm8D(r,s)
#endif

#ifndef vload
#define vload(x) (*(vdr*)(x))
#define vstore(x,y) (*(vdr*)(x)) = (y)
#endif
#define LD(x) vload(x)
#define ST(x,y) vstore((x),(y))
//#define PEQ(x,y) (x) += (y)
//#define MEQ(x,y) (x) -= (y)
#ifndef MUL
#define NEG(x) (-(x))
#define MUL(x,y) ((x)*(y))
#define MADD(x,y,z) (((x)*(y))+(z))
#define NMADD(x,y,z) ((z)-((x)*(y)))
#endif

extern int vlen;

void X_veq_perm_xX(double *rr, int perm, int *idx, double *vv, int nelem, int off, int len);
void xX_veq_yX(int *idx1, double *rr, int *idx2, double *vv, int nelem, int off, int len);
void xX_veq_X(int *idx1, double *rr, double *vv, int nelem, int off, int len);
void X_veq_xXc(double *rr, int *idx, double *vv, int nelem, int off, int len);
void X_veq_xXn(double *rr, int *idx, double *vv, int nelem, int off, int len);
void X_veq_xXp(double *rr, int *idx1, int perm, int *idx2, double *vv,
	       int nelem, int off, int len);
void X_veq_xXp2(double *rr, int *idx1, int perm, double *vv,
		int nelem, int off, int len);


void V_veq_perm_xV(double *restrict rr, int perm, int *restrict idx, double *restrict vv, int off, int len);
void V_veq_r_times_V(double *restrict rr, double aa, double *restrict vv, int off, int len);
void
V_vpeq_M_times_xV(double *restrict rr, double *restrict mm, int *restrict idx, double *restrict vv,
		  int off, int len);
void
V_vpeq_M_times_xVc(double *restrict rr, double *restrict mm, int *restrict idx,
		   double *restrict vv, int off, int len);
void
V_vpeq_M_times_xVp(double *restrict rr, double *restrict mm, int *restrict idx,
		   int perm, double *restrict vv, int i0, int i1);
void
V_vpeq_M_times_xVn(double *restrict rr, double *restrict mm, int *restrict idx,
		   double *restrict vv, int off, int len);
void
V_vpeq_M_times_xVb(double *rr, double *mm, int blend, int *idx1, double *vv1,
		   int *idx2, double *vv2, int i0, int i1);
void
xV_vpeq_M_times_V(int *restrict idx1, double *restrict rr, double *restrict mm, double *restrict vv, int off, int len);
void
xV_vpeq_xM_times_V(int *restrict idx, double *restrict rr, double *restrict mm, double *restrict vv, int off, int len);
void
xV_vpeq_xM_times_yV(int *restrict idx1, double *restrict rr, double *restrict mm, int *restrict idx2, double *restrict vv, int off, int len);
void
xV_vpeq_yM_times_yV(int *restrict idx1, double *restrict rr, int *restrict idx2, double *restrict mm, double *restrict vv, int begin, int end);
void
xV_vpeq_xM_times_yVzVb(int *idx0, double *rr, double *mm, int blend, int *idx1, double *vv1,
		       int *idx2, double *vv2, int i0, int i1);

void
X_veq_pack_xX(double *rr, int pack, int *idx, double *vv,
	      int nelem, int off, int len);
void
xX_veq_blend_xXX(int *idx1, double *rr, int blend, int *idx2,
		 double *vv1, double *vv2, int nelem, int off, int len);
void
xX_veq_blend_xXxX(int *idx0, double *rr, int blend, int *idx1, double *vv1,
		  int *idx2, double *vv2, int nelem, int off, int len);
void
X_veq_blend_xXxXn(double *rr, int blend, int *idx1, double *vv1,
		  int *idx2, double *vv2, int nelem, int off, int len);

void r_veq_norm2_X(double *r, double *vv, int nelem, int i0, int i1);
void r_veq_re_X_dot_X(double *r, double *vv, double *ww, int nelem, int i0, int i1);
void X_veq_r(double *rr, double aa, int nelem, int i0, int i1);
void X_veq_r_times_X(double *rr, double aa, double *vv,
		     int nelem, int i0, int i1);
void X_veq_X_plus_r_times_X(double *rr, double *vv, double aa, double *ww,
			    int nelem, int i0, int i1);
