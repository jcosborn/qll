#ifdef C1
#define VLEN 1
#define MAXVLENF 1
#define MAXVLEND 1
#endif
#ifdef C2
#define VLEN 2
#define MAXVLENF 2
#define MAXVLEND 2
#endif
#ifdef C4
#define VLEN 4
#define MAXVLENF 4
#define MAXVLEND 4
#endif
#ifdef C8
#define VLEN 8
#define MAXVLENF 8
#define MAXVLEND 8
#endif
#ifdef SSE
#ifdef DOUBLE
#define VLEN 2
#else
#define VLEN 4
#endif
#define MAXVLENF 4
#define MAXVLEND 2
#endif
#ifdef QPX
#define VLEN 4
#define MAXVLENF 4
#define MAXVLEND 4
#endif
#ifdef AVX
#ifdef DOUBLE
#define VLEN 4
#else
#define VLEN 8
#endif
#define MAXVLENF 8
#define MAXVLEND 4
#endif
#ifdef MIC
#ifdef DOUBLE
#define VLEN 8
#else
#define VLEN 16
#endif
#define MAXVLENF 16
#define MAXVLEND 8
#endif

#ifndef QPX
#define HAVE_GCC_VECTOR_TYPE
#ifndef __INTEL_COMPILER
#define CAN_INDEX_VECTOR
#endif
#endif

#define TVCF(s,n) typedef struct { float v[n]; } s	\
  __attribute__ ((aligned(4*n)))
#define TVCD(s,n) typedef struct { double v[n]; } s	\
  __attribute__ ((aligned(8*n)))
#define ELEMCF(x,i) ((x).v[i])
#define ELEMCD(x,i) ((x).v[i])
#define INITC(s,n,...) ((s){{CAT(CATLIST,n)(__VA_ARGS__)}})
#ifdef CAN_INDEX_VECTOR
#define ELEMGF(x,i) ((x)[i])
#define ELEMGD(x,i) ((x)[i])
#else
#define ELEMGF(x,i) (((float*)&(x))[i])
#define ELEMGD(x,i) (((double*)&(x))[i])
#endif
#define INITG(s,n,...) ((s){CAT(CATLIST,n)(__VA_ARGS__)})

#ifdef HAVE_GCC_VECTOR_TYPE
#define TVF(s,n) typedef float s					\
  __attribute__ ((vector_size(4*n))) __attribute__ ((aligned(4*n)))
#define TVD(s,n) typedef double s					\
  __attribute__ ((vector_size(8*n))) __attribute__ ((aligned(8*n)))
#define TVDA(s,n,m) typedef struct { \
  double x[m] __attribute__ ((vector_size(8*n))) __attribute__ ((aligned(8*n))); \
  } s
#define ELEMF(x,i) ELEMGF(x,i)
#define ELEMD(x,i) ELEMGD(x,i)
#define INIT(s,n,...) INITG(s,n,__VA_ARGS__)
#else
#define TVF(s,n) TVCF(s,n)
#define TVD(s,n) TVCD(s,n)
#define ELEMF(x,i) ELEMCF(x,i)
#define ELEMD(x,i) ELEMCD(x,i)
#define INIT(s,n,...) INITC(s,n,__VA_ARGS__)
#endif
#define CATLIST1(a) (a)
#define CATLIST2(a,b) (a),(b)
#define CATLIST4(a,b,c,d) (a),(b),(c),(d)
#define CATLIST8(a,b,c,d,...) CATLIST4(a,b,c,d),CATLIST4(__VA_ARGS__)
#define CATLIST16(a,b,c,d,e,f,g,h,...) CATLIST8(a,b,c,d,e,f,g,h),CATLIST8(__VA_ARGS__)
#define CATLIST32(a,b,c,d,e,f,g,h,i,j,k,l,m,n,o,p,...) CATLIST16(a,b,c,d,e,f,g,h,i,j,k,l,m,n,o,p),CATLIST16(__VA_ARGS__)
#define CAT(a,b) CATX(a,b)
#define CATX(a,b) a ## b

#if MAXVLENF >= 1
TVF(vec1f, 1);
#define v1felem(x,i) ELEMF(x,i)
#define v1finit(...) INIT(vec1fr,1,__VA_ARGS__)
#else
TVCF(vec1f, 1);
#define v1felem(x,i) ELEMCF(x,i)
#define v1finit(...) INITC(vec1fr,1,__VA_ARGS__)
#endif
#if MAXVLEND >= 1
TVD(vec1d, 1);
#define v1delem(x,i) ELEMD(x,i)
#define v1dinit(...) INIT(vec1dr,1,__VA_ARGS__)
#else
TVCD(vec1d, 1);
#define v1delem(x,i) ELEMCD(x,i)
#define v1dinit(...) INITC(vec1dr,1,__VA_ARGS__)
#endif

#if MAXVLENF >= 2
TVF(vec2f, 2);
#define v2felem(x,i) ELEMF(x,i)
#define v2finit(...) INIT(vec2fr,2,__VA_ARGS__)
#else
TVCF(vec2f, 2);
#define v2felem(x,i) ELEMCF(x,i)
#define v2finit(...) INITC(vec2fr,2,__VA_ARGS__)
#endif
#if MAXVLEND >= 2
TVD(vec2d, 2);
#define v2delem(x,i) ELEMD(x,i)
#define v2dinit(...) INIT(vec2dr,2,__VA_ARGS__)
#else
TVCD(vec2d, 2);
#define v2delem(x,i) ELEMCD(x,i)
#define v2dinit(...) INITC(vec2dr,2,__VA_ARGS__)
#endif

#if MAXVLENF >= 4
TVF(vec4f, 4);
#ifdef QPX
#define v4felem(x,i) ELEMGF(x,i)
#define v4finit(...) INITGF(vec4fr,4,__VA_ARGS__)
#else
#define v4felem(x,i) ELEMF(x,i)
#define v4finit(...) INIT(vec4fr,4,__VA_ARGS__)
#endif
#else
TVCF(vec4f, 4);
#define v4felem(x,i) ELEMCF(x,i)
#define v4finit(...) INITC(vec4fr,4,__VA_ARGS__)
#endif
#if MAXVLEND >= 4
TVD(vec4d, 4);
#ifdef QPX
//typedef __vector4double vec4d;
#define v4delem(x,i) ELEMGD(x,i)
#define v4dinit(...) INITGD(vec4dr,4,__VA_ARGS__)
#else
#define v4delem(x,i) ELEMD(x,i)
#define v4dinit(...) INIT(vec4dr,4,__VA_ARGS__)
#endif
#else
TVCD(vec4d, 4);
#define v4delem(x,i) ELEMCD(x,i)
#define v4dinit(...) INITC(vec4dr,4,__VA_ARGS__)
#endif

#if MAXVLENF >= 8
TVF(vec8f, 8);
#define v8felem(x,i) ELEMF(x,i)
#define v8finit(...) INIT(vec8fr,8,__VA_ARGS__)
#else
TVCF(vec8f, 8);
#define v8felem(x,i) ELEMCF(x,i)
#define v8finit(...) INITC(vec8fr,8,__VA_ARGS__)
#endif
#if MAXVLEND >= 8
TVD(vec8d, 8);
#define v8delem(x,i) ELEMD(x,i)
#define v8dinit(...) INIT(vec8dr,8,__VA_ARGS__)
#else
TVCD(vec8d, 8);
#define v8delem(x,i) ELEMCD(x,i)
#define v8dinit(...) INITC(vec8dr,8,__VA_ARGS__)
#endif

#if MAXVLENF >= 16
TVF(vec16f, 16);
#define v16felem(x,i) ELEMF(x,i)
#define v16finit(...) INIT(vec16fr,16,__VA_ARGS__)
#else
TVCF(vec16f, 16);
#define v16felem(x,i) ELEMCF(x,i)
#define v16finit(...) INITC(vec16fr,16,__VA_ARGS__)
#endif
#if MAXVLEND >= 16
TVD(vec16d, 16);
#define v16delem(x,i) ELEMD(x,i)
#define v16dinit(...) INIT(vec16dr,16,__VA_ARGS__)
#else
TVCD(vec16d, 16);
#define v16delem(x,i) ELEMCD(x,i)
#define v16dinit(...) INITC(vec16dr,16,__VA_ARGS__)
#endif

#if MAXVLENF >= 32
TVF(vec32f, 32);
#define v32felem(x,i) ELEMF(x,i)
#define v32finit(...) INIT(vec32fr,32,__VA_ARGS__)
#else
TVCF(vec32f, 32);
#define v32felem(x,i) ELEMCF(x,i)
#define v32finit(...) INITC(vec32fr,32,__VA_ARGS__)
#endif
#if MAXVLEND >= 32
TVD(vec32d, 32);
#define v32delem(x,i) ELEMD(x,i)
#define v32dinit(...) INIT(vec32dr,32,__VA_ARGS__)
#else
TVCD(vec32d, 32);
#define v32delem(x,i) ELEMCD(x,i)
#define v32dinit(...) INITC(vec32dr,32,__VA_ARGS__)
#endif

#ifdef QPX
typedef __vector4double vec4fr;
typedef __vector4double vec4dr;
typedef vec1d  vec1fr;
typedef vec2d  vec2fr;
//typedef vec4d  vec4fr;
typedef vec8d  vec8fr;
typedef vec16d vec16fr;
typedef vec32d vec32fr;
#else
typedef vec1f  vec1fr;
typedef vec2f  vec2fr;
typedef vec4f  vec4fr;
typedef vec8f  vec8fr;
typedef vec16f vec16fr;
typedef vec32f vec32fr;
typedef vec4d  vec4dr;
#endif

typedef vec1d  vec1dr;
typedef vec2d  vec2dr;
typedef vec8d  vec8dr;
typedef vec16d vec16dr;
typedef vec32d vec32dr;

#define CAT3(a,b,c) CAT3X(a,b,c)
#define CAT3X(a,b,c) a ## b ## c

#ifdef DOUBLE

typedef double real;
typedef CAT3(vec,VLEN,d)  vec;
typedef CAT3(vec,VLEN,dr) vecr;

#else

typedef float real;
typedef CAT3(vec,VLEN,f)  vec;
typedef CAT3(vec,VLEN,fr) vecr;

#endif

typedef struct {
  vec r,i;
} cvec;

typedef struct {
  vecr r,i;
} cvecr;
