/**
 *ld,*st,*splat
 *neg,*add,*sub,*mul,*madd,*nmadd
 *red
 *perm1,2,4,8
 *packp1,*packm1,2,4,8
 *blendp1,*blendm1,2,4,8
 vf8ldvf4(a)
 vf8ldvf4x2(a,b)
 vf4stvf8(a,x)
 vf4x2stvf8(a,b,x)
**/

#ifdef QPX
#include <vec_ops_qpx.h>
#endif
#ifdef AVX
#include <vec_ops_avx.h>
#endif
#ifdef MIC
#include <vec_ops_mic.h>
#endif
#ifdef USE_GCC_VECTOR_OPS
#include <vec_ops_gcc.h>
#endif
#include <vec_ops_c.h>

#ifdef DOUBLE
#define PREC d
#else
#define PREC f
#endif

#define CAT(a,b) CATX(a,b)
#define CATX(a,b) a ## b

#define VF(x) CAT(CAT3(v,VLEN,PREC),x)

#define LD(x) VF(ld)(x)
#define ST(x,y) VF(st)(x,y)
#define STNR(x,y) VF(stnr)(x,y)

#define perm1(x) VF(perm1)(x)
#define perm2(x) VF(perm2)(x)
#define perm4(x) VF(perm4)(x)
#define perm8(x) VF(perm8)(x)

#define packp1(x,y) VF(packp1)(x,y)
#define packm1(x,y) VF(packm1)(x,y)
#define packp2(x,y) VF(packp2)(x,y)
#define packm2(x,y) VF(packm2)(x,y)
#define packp4(x,y) VF(packp4)(x,y)
#define packm4(x,y) VF(packm4)(x,y)
#define packp8(x,y) VF(packp8)(x,y)
#define packm8(x,y) VF(packm8)(x,y)

#define blendp1(x,y) VF(blendp1)(x,y)
#define blendm1(x,y) VF(blendm1)(x,y)
#define blendp2(x,y) VF(blendp2)(x,y)
#define blendm2(x,y) VF(blendm2)(x,y)
#define blendp4(x,y) VF(blendp4)(x,y)
#define blendm4(x,y) VF(blendm4)(x,y)
#define blendp8(x,y) VF(blendp8)(x,y)
#define blendm8(x,y) VF(blendm8)(x,y)

#define splat(x) VF(splat)(x)
#define vreduce(x) VF(reduce)(x)
#define NEG(x) VF(neg)(x)
#define ADD(x,y) VF(add)(x,y)
#define SUB(x,y) VF(sub)(x,y)
#define MUL(x,y) VF(mul)(x,y)
#define MADD(x,y,z) VF(madd)(x,y,z)
#define NMADD(x,y,z) VF(nmadd)(x,y,z)
#define INVSQRT(x) VF(rsqrt)(x)
