
#ifndef v2fperm1
#define v2fperm1(s) __builtin_shuffle(s, (int __attribute__((vector_size(8)))){1,0})
#endif

#ifndef v4fperm1
#define v4fperm1(s) __builtin_shuffle(s, (int __attribute__((vector_size(16)))){1,0,3,2})
#endif
#ifndef v4fperm2
#define v4fperm2(s) __builtin_shuffle(s, (int __attribute__((vector_size(16)))){2,3,0,1})
#endif

#ifndef v8fperm1
#define v8fperm1(s) __builtin_shuffle(s,(int __attribute__((vector_size(32)))){1,0,3,2,5,4,7,6})
#endif
#ifndef v8fperm2
#define v8fperm2(s) __builtin_shuffle(s,(int __attribute__((vector_size(32)))){2,3,0,1,6,7,4,5})
#endif
#ifndef v8fperm4
#define v8fperm4(s) __builtin_shuffle(s,(int __attribute__((vector_size(32)))){4,5,6,7,0,1,2,3})
#endif

#ifndef v4dperm1
#define v4dperm1(s) __builtin_shuffle(s, (long __attribute__((vector_size(32)))){1,0,3,2})
#endif
#ifndef v4dperm2
#define v4dperm2(s) __builtin_shuffle(s, (long __attribute__((vector_size(32)))){2,3,0,1})
#endif

#ifndef v8dperm1
#define v8dperm1(s) __builtin_shuffle(s,(long __attribute__((vector_size(64)))){1,0,3,2,5,4,7,6})
#endif
#ifndef v8dperm2
#define v8dperm2(s) __builtin_shuffle(s,(long __attribute__((vector_size(64)))){2,3,0,1,6,7,4,5})
#endif
#ifndef v8dperm4
#define v8dperm4(s) __builtin_shuffle(s,(long __attribute__((vector_size(64)))){4,5,6,7,0,1,2,3})
#endif

#ifndef v8fadd
#define v8fadd(x,y) ((x)+(y))
#endif
#ifndef v8fsub
#define v8fsub(x,y) ((x)-(y))
#endif
#ifndef v8fmul
#define v8fmul(x,y) ((x)*(y))
#endif
#ifndef v8fmadd
#define v8fmadd(x,y,z) ((x)*(y)+(z))
#endif
#ifndef v8fnmadd
#define v8fnmadd(x,y,z) ((z)-(x)*(y))
#endif

#ifndef v4dadd
#define v4dadd(x,y) ((x)+(y))
#endif
#ifndef v4dsub
#define v4dsub(x,y) ((x)-(y))
#endif
#ifndef v4dmul
#define v4dmul(x,y) ((x)*(y))
#endif
#ifndef v4dmadd
#define v4dmadd(x,y,z) ((x)*(y)+(z))
#endif
#ifndef v4dnmadd
#define v4dnmadd(x,y,z) ((z)-(x)*(y))
#endif
