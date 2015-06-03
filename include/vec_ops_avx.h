#include <immintrin.h>

#define BASE2(a,b,c,d) (8*(a)+4*(b)+2*(c)+(d))
#define BASE4(a,b,c,d) (64*(a)+16*(b)+4*(c)+(d))

#define v8fsplat(x) _mm256_broadcast_ss((float[]){(float)(x)})
#define v8fperm1(s) _mm256_shuffle_ps(s, s, BASE4(2,3,0,1))
#define v8fperm2(s) _mm256_shuffle_ps(s, s, BASE4(1,0,3,2))
#define v8fperm4(s) _mm256_permute2f128_ps(s, s, 1)

#define v4dsplat(x) _mm256_broadcast_sd((double[]){(double)(x)})
#define v4dperm1(s) _mm256_shuffle_pd(s, s, BASE2(0,1,0,1))
#define v4dperm2(s) _mm256_permute2f128_pd(s, s, BASE4(0,0,0,1))
