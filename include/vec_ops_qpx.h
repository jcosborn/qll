#ifndef v4fperm1
#define v4fperm1(s) vec_perm(s,s,vec_gpci(01032))
#endif
#ifndef v4fperm2
#define v4fperm2(s) vec_perm(s,s,vec_gpci(02301))
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
#ifndef v4fld
#define v4fld(x) vec_ld(0,(x)->v)
#endif
#ifndef v4fst
#define v4fst(x,y) vec_st(vec_rsp(y),0,(x)->v)
#endif
#ifndef v4fstnr
#define v4fstnr(x,y) vec_st(y,0,(x)->v)
#endif
#ifndef v4fsplat
#define v4fsplat(x) vec_splats((double)(x))
#endif
#ifndef v4fneg
#define v4fneg(x) vec_neg(x)
#endif
#ifndef v4fadd
#define v4fadd(x,y) vec_add(x,y)
#endif
#ifndef v4fsub
#define v4fsub(x,y) vec_sub(x,y)
#endif
#ifndef v4fmul
#define v4fmul(x,y) vec_mul(x,y)
#endif
#ifndef v4fmadd
#define v4fmadd(x,y,z) vec_madd(x,y,z)
#endif
#ifndef v4fnmadd
#define v4fnmadd(x,y,z) vec_nmsub(x,y,z)
#endif

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
#ifndef v4dld
#define v4dld(x) vec_ld(0,(x)->v)
#endif
#ifndef v4dst
#define v4dst(x,y) vec_st(y,0,(x)->v)
#endif
#ifndef v4dstnr
#define v4dstnr(x,y) vec_st(y,0,(x)->v)
#endif
#ifndef v4dsplat
#define v4dsplat(x) vec_splats((double)(x))
#endif
#ifndef v4dneg
#define v4dneg(x) vec_neg(x)
#endif
#ifndef v4dadd
#define v4dadd(x,y) vec_add(x,y)
#endif
#ifndef v4dsub
#define v4dsub(x,y) vec_sub(x,y)
#endif
#ifndef v4dmul
#define v4dmul(x,y) vec_mul(x,y)
#endif
#ifndef v4dmadd
#define v4dmadd(x,y,z) vec_madd(x,y,z)
#endif
#ifndef v4dnmadd
#define v4dnmadd(x,y,z) vec_nmsub(x,y,z)
#endif
