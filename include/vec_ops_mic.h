#include <immintrin.h>
#define BASE4(a,b,c,d) (64*(a)+16*(b)+4*(c)+(d))

#ifndef v16fperm1
#define v16fperm1(s) _mm512_swizzle_ps(s, _MM_SWIZ_REG_CDAB)
#endif
#ifndef v16fperm2
#define v16fperm2(s) _mm512_swizzle_ps(s, _MM_SWIZ_REG_BADC)
#endif
#ifndef v16fperm4
#define v16fperm4(s) _mm512_permute4f128_ps(s, BASE4(2,3,0,1))
#endif
#ifndef v16fperm8
#define v16fperm8(s) _mm512_permute4f128_ps(s, BASE4(1,0,3,2))
#endif

#ifndef v4fpackp1
#define v4fpackp1(r,s) vec_st2(vec_perm(s,s,vec_gpci(01302)),0,r)
#endif
#ifndef v4fpackm1
#define v4fpackm1(r,s) vec_st2(vec_perm(s,s,vec_gpci(00213)),0,r)
#endif
#ifndef v4fpackp2
#define v4fpackp2(r,s) vec_st2(vec_perm(s,s,vec_gpci(02301)),0,r)
#endif
#ifndef v4fpackm2
#define v4fpackm2(r,s) vec_st2(s,0,r)
#endif
#ifndef v4fblendp1
#define v4fblendp1(s,t) vec_perm(s,vec_ld2(0,t),vec_gpci(04052))
#endif
#ifndef v4fblendm1
#define v4fblendm1(s,t) vec_perm(s,vec_ld2(0,t),vec_gpci(01435))
#endif
#ifndef v4fblendp2
#define v4fblendp2(s,t) vec_perm(s,vec_ld2(0,t),vec_gpci(04501))
#endif
#ifndef v4fblendm2
#define v4fblendm2(s,t) vec_perm(s,vec_ld2(0,t),vec_gpci(02345))
#endif

#ifndef v16fld
#define v16fld(x) _mm512_load_ps(x)
#endif
#ifndef v16fst
#define v16fst(x,y) _mm512_store_ps(x,y)
#endif
#ifndef v16fstnr
#define v16fstnr(x,y) _mm512_store_ps(x,y)
#endif
#ifndef v16fsplat
#define v16fsplat(x) _mm512_extload_ps((float[]){(float)(x)}, _MM_UPCONV_PS_NONE, _MM_BROADCAST_1X16, 0)
#endif
//#ifndef v16fneg
//#define v16fneg(x) _mm512_neg_ps(x)
//#endif
#ifndef v16fadd
#define v16fadd(x,y) _mm512_add_ps(x,y)
#endif
#ifndef v16fsub
#define v16fsub(x,y) _mm512_sub_ps(x,y)
#endif
#ifndef v16fmul
#define v16fmul(x,y) _mm512_mul_ps(x,y)
#endif
#ifndef v16fmadd
#define v16fmadd(x,y,z) _mm512_fmadd_ps(x,y,z)
#endif
#ifndef v16fnmadd
#define v16fnmadd(x,y,z) _mm512_fnmadd_ps(x,y,z)
#endif

#if 0
#ifndef v4dperm1
#define v4dperm1(s) vec_perm(s,s,vec_gpci(01032))
#endif
#ifndef v4dperm2
#define v4dperm2(s) vec_perm(s,s,vec_gpci(02301))
#endif
#ifndef v4dpackp1
#define v4dpackp1(r,s) vec_st2(vec_perm(s,s,vec_gpci(01302)),0,r)
#endif
#ifndef v4dpackm1
#define v4dpackm1(r,s) vec_st2(vec_perm(s,s,vec_gpci(00213)),0,r)
#endif
#ifndef v4dpackp2
#define v4dpackp2(r,s) vec_st2(vec_perm(s,s,vec_gpci(02301)),0,r)
#endif
#ifndef v4dpackm2
#define v4dpackm2(r,s) vec_st2(s,0,r)
#endif
#ifndef v4dblendp1
#define v4dblendp1(s,t) vec_perm(s,vec_ld2(0,t),vec_gpci(04052))
#endif
#ifndef v4dblendm1
#define v4dblendm1(s,t) vec_perm(s,vec_ld2(0,t),vec_gpci(01435))
#endif
#ifndef v4dblendp2
#define v4dblendp2(s,t) vec_perm(s,vec_ld2(0,t),vec_gpci(04501))
#endif
#ifndef v4dblendm2
#define v4dblendm2(s,t) vec_perm(s,vec_ld2(0,t),vec_gpci(02345))
#endif
#endif

#ifndef v8dld
#define v8dld(x) _mm512_load_pd(x)
#endif
#ifndef v8dst
#define v8dst(x,y) _mm512_store_pd(x,y)
#endif
#ifndef v8dstnr
#define v8dstnr(x,y) _mm512_store_pd(x,y)
#endif
#ifndef v8dsplat
#define v8dsplat(x) _mm512_extload_pd((double[]){(double)(x)}, _MM_UPCONV_PD_NONE, _MM_BROADCAST_1X8, 0)
#endif
//#ifndef v4dneg
//#define v4dneg(x) vec_neg(x)
//#endif
#ifndef v8dadd
#define v8dadd(x,y) _mm512_add_pd(x,y)
#endif
#ifndef v8dmul
#define v8dmul(x,y) _mm512_mul_pd(x,y)
#endif
#ifndef v8dmadd
#define v8dmadd(x,y,z) _mm512_fmadd_pd(x,y,z)
#endif
#ifndef v8dnmadd
#define v8dnmadd(x,y,z) _mm512_fnmadd_pd(x,y,z)
#endif
