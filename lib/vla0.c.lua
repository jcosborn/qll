function printf(...) io.write(string.format(...)) end

function LC(v, e)
  printf("  vecr %sr = LD(&((%s).r));\n", v, e)
  printf("  vecr %si = LD(&((%s).i));\n", v, e)
end

function LCP(v, e, p)
  printf("  vecr %sr = permN(LD(&((%s).r)),%s);\n", v, e, p)
  printf("  vecr %si = permN(LD(&((%s).i)),%s);\n", v, e, p)
end

function LCB(v, e1, e2, b)
  printf("  vecr %sr = blendN(LD(&((%s).r)),(%s),%s);\n", v, e1,e2, b)
  printf("  vecr %si = blendN(LD(&((%s).i)),(%s)+(VLEN/2),%s);\n", v, e1,e2, b)
end

function LCBi(v, e1, e2, b)
  printf("  vecr %sr = blend%s(LD(&((%s).r)),(%s));\n",v,b,e1,e2)
  printf("  vecr %si = blend%s(LD(&((%s).i)),(%s)+(VLEN/2));\n",v,b,e1,e2)
end

function LV(v, e)
  for i=0,2 do
    printf("  vecr %s%ir = LD(&((%s).c%i.r));\n", v, i, e, i)
    printf("  vecr %s%ii = LD(&((%s).c%i.i));\n", v, i, e, i)
  end
end

function LVZ(v)
  printf("  real zero = 0;\n")
  for i=0,2 do
    printf("  vecr %s%ir = splat(zero);\n", v, i, e, i)
    printf("  vecr %s%ii = splat(zero);\n", v, i, e, i)
  end
end

function AVP(r, v, p)
  for i=0,2 do
    printf("  vecr %s%irp = permN(%s%ir,%s);\n", v, i, v, i, p)
    printf("  %s%ir = ADD(%s%ir,%s%irp);\n", r, i, r, i, v, i, p)
    printf("  vecr %s%iip = permN(%s%ii,%s);\n", v, i, v, i, p)
    printf("  %s%ii = ADD(%s%ii,%s%iip);\n", r, i, r, i, v, i, p)
  end
end

function LMC(v, m, e)
  for i=0,2 do
    printf("  vecr %s%ir = LD(&((%s)->v%i.c%i.r));\n", v, i, m, i, e)
    printf("  vecr %s%ii = LD(&((%s)->v%i.c%i.i));\n", v, i, m, i, e)
  end
end

function SV(v, e)
  for i=0,2 do
    printf("  ST(&((%s).c%i.r), %s%ir);\n", v, i, e, i)
    printf("  ST(&((%s).c%i.i), %s%ii);\n", v, i, e, i)
  end
end

function PEQTIMES(x, y, z)
  --printf("  PEQ(%s,MUL(%s,%s));\n", x, y, z)
  printf("  %s = MADD(%s,%s,%s);\n", x, y, z, x)
end

function MEQTIMES(x, y, z)
  --printf("  MEQ(%s,MUL(%s,%s));\n", x, y, z)
  printf("  %s = NMADD(%s,%s,%s);\n", x, y, z, x)
end

------ start output ------

for i=0,3 do
  local k = 2^i
  printf([[
#if VLEN > %i
static inline void
Perm%i(vec *x, vec *y, int n)
{
  for(int i=0; i<n; i++) {
    vecr yy = LD(&y[i]);
    vecr xx = perm%i(yy);
    STNR(&x[i], xx);
  }
}
]],k,k,k)
printf([[
static inline void
Packp%i(real *x, vec *y, int n)
{
  int nv2 = VLEN/2;
  for(int i=0; i<n; i++) {
    vecr yy = LD(&y[i]);
    packp%i(&x[nv2*i], yy);
  }
}
static inline void
Packm%i(real *x, vec *y, int n)
{
  int nv2 = VLEN/2;
  for(int i=0; i<n; i++) {
    vecr yy = LD(&y[i]);
    packm%i(&x[nv2*i], yy);
  }
}
]],k,k,k,k,k,k,k)
printf([[
static inline void
Blendp%i(vec *x, vec *y, real *z, int n)
{
  int nv2 = VLEN/2;
  for(int i=0; i<n; i++) {
    vecr yy = LD(&y[i]);
    vecr xx = blendp%i(yy, &z[nv2*i]);
    STNR(&x[i], xx);
  }
}
static inline void
Blendm%i(vec *x, vec *y, real *z, int n)
{
  int nv2 = VLEN/2;
  for(int i=0; i<n; i++) {
    vecr yy = LD(&y[i]);
    vecr xx = blendm%i(yy, &z[nv2*i]);
    STNR(&x[i], xx);
  }
}
#endif
]],k,k,k,k,k,k)
end

printf([[
static inline vecr
permN(vecr x, int p)
{
  //if(p==0) return x;
  //vecr y = x;
  vecr y;
  switch(p) {
]])
for i=0,3 do
  local k = 2^i
printf([[
#if VLEN > %i
  case %i: y = perm%i(x); break;
#endif
]],k,k,k)
end
printf([[
  default: y = x;
  }
  return y;
}
]])

printf([[
static inline vecr
blendN(vecr y, real *z, int b)
{
  vecr x;
  switch(b) {
]])
for i=0,3 do
  local k = 2^i
printf([[
#if VLEN > %i
  case  %i: x = blendp%i(y,z); break;
  case -%i: x = blendm%i(y,z); break;
#endif
]],k,k,k,k,k)
end
printf([[
    default: x = LD((vec*)z);
  }
  return x;
}
]])

printf([[
static inline void
RtimesV3(Vec3 *x, vecr c, Vec3 *y)
{
]])
LV("y","*y")
for i=0,2 do
  printf("  vecr x%ir = MUL(c,y%ir);\n", i, i)
  printf("  ST(&(x->c%i.r), x%ir);\n", i, i)
  printf("  vecr x%ii = MUL(c,y%ii);\n", i, i)
  printf("  ST(&(x->c%i.i), x%ii);\n", i, i)
end
printf("}\n")

printf([[
static inline void
PeqM3AtimesV3o(Vec3 *r, Mat3 *m, Vec3 *v)
{
]])
LV("r","*r")
for i=0,2 do
  LC("v"..i,"v->c"..i)
  LV("m"..i,"m->v"..i)
  for j=0,2 do
    PEQTIMES("r"..j.."r", "v"..i.."r", "m"..i..j.."r")
    PEQTIMES("r"..j.."i", "v"..i.."r", "m"..i..j.."i")
  end
  for j=0,2 do
    MEQTIMES("r"..j.."r", "v"..i.."i", "m"..i..j.."i")
    PEQTIMES("r"..j.."i", "v"..i.."i", "m"..i..j.."r")
  end
end
SV("*r","r")
printf("}\n")

printf([[
static inline void
PeqM3AtimesV3op(Vec3 *r, Mat3 *m, Vec3 *v, int perm)
{
]])
LV("r","*r")
for i=0,2 do
  LCP("v"..i,"v->c"..i,"perm")
  LV("m"..i,"m->v"..i)
  for j=0,2 do
    PEQTIMES("r"..j.."r", "v"..i.."r", "m"..i..j.."r")
    PEQTIMES("r"..j.."i", "v"..i.."r", "m"..i..j.."i")
  end
  for j=0,2 do
    MEQTIMES("r"..j.."r", "v"..i.."i", "m"..i..j.."i")
    PEQTIMES("r"..j.."i", "v"..i.."i", "m"..i..j.."r")
  end
end
SV("*r","r")
printf("}\n")

printf([[
static inline void
PeqM3AtimesV3ob(Vec3 *r, Mat3 *m, Vec3 *vv1, real *vv2, int blend)
{
]])
LV("r","*r")
for i=0,2 do
  LCB("v"..i,"vv1->c"..i,"&vv2["..i.."*VLEN]","blend")
  LV("m"..i,"m->v"..i)
  for j=0,2 do
    PEQTIMES("r"..j.."r", "v"..i.."r", "m"..i..j.."r")
    PEQTIMES("r"..j.."i", "v"..i.."r", "m"..i..j.."i")
  end
  for j=0,2 do
    MEQTIMES("r"..j.."r", "v"..i.."i", "m"..i..j.."i")
    PEQTIMES("r"..j.."i", "v"..i.."i", "m"..i..j.."r")
  end
end
SV("*r","r")
printf("}\n")

for i=0,3 do
  for s=1,-1,-2 do
    local p = 2^i
    b = (s==1) and ("p"..p) or ("m"..p)
    if(s==1) then printf("#if VLEN > %i\n", p) end
    printf([[
static inline void
PeqM3AtimesV3ob%s(Vec3 *r, Mat3 *m, Vec3 *vv1, real *vv2)
{
]],b)
    LV("r","*r")
    for i=0,2 do
      LCBi("v"..i,"vv1->c"..i,"&vv2["..i.."*VLEN]",b)
      LV("m"..i,"m->v"..i)
      for j=0,2 do
	PEQTIMES("r"..j.."r", "v"..i.."r", "m"..i..j.."r")
	PEQTIMES("r"..j.."i", "v"..i.."r", "m"..i..j.."i")
      end
      for j=0,2 do
	MEQTIMES("r"..j.."r", "v"..i.."i", "m"..i..j.."i")
	PEQTIMES("r"..j.."i", "v"..i.."i", "m"..i..j.."r")
      end
    end
    SV("*r","r")
    printf("}\n")
    if(s==-1) then printf("#endif\n", i) end
  end
end

printf([[
void
P(V_vpeq_M_times_xVb)(real *rr, real *mm, int blend, int *idx1, real *vv1,
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
]])
for i=0,3 do
  for s=1,-1,-2 do
    local p = 2^i
    v = s*p
    b = (s==1) and ("p"..p) or ("m"..p)
    if(s==1) then printf("#if VLEN > %i\n", p) end
    printf([[
  } else if(blend==%s) {
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
	PeqM3AtimesV3ob%s(ri, mi, v1i, v2i);
      }
    }
]],v,b)
    if(s==-1) then printf("#endif\n", i) end
  end
end
printf([[
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
]])

printf([[
void
P(xV_vpeq_xM_times_yVzVb)(int *idx0, real *rr, real *mm, int blend, int *idx1,
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
]])
for i=0,3 do
  for s=1,-1,-2 do
    local p = 2^i
    v = s*p
    b = (s==1) and ("p"..p) or ("m"..p)
    if(s==1) then printf("#if VLEN > %i\n", p) end
    printf([[
  } else if(blend==%i) {
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
      PeqM3AtimesV3ob%s(ri, mi, v1i, v2i);
    }
]],v,b)
    if(s==-1) then printf("#endif\n", i) end
  end
end
printf([[
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
]])

printf([[
static inline void
MeqM3AtimesV3(Vec3 *r, Mat3 *m, Vec3 *v)
{
]])
LV("r","*r")
for i=0,2 do
  LC("v"..i,"v->c"..i)
  --LV("m"..i,"m->v"..i)
  LMC("m"..i,"m",i)
  for j=0,2 do
    MEQTIMES("r"..j.."r", "v"..i.."r", "m"..i..j.."r")
    PEQTIMES("r"..j.."i", "v"..i.."r", "m"..i..j.."i")
  end
  for j=0,2 do
    MEQTIMES("r"..j.."r", "v"..i.."i", "m"..i..j.."i")
    MEQTIMES("r"..j.."i", "v"..i.."i", "m"..i..j.."r")
  end
end
SV("*r","r")
printf("}\n")

printf([[
static inline void
MeqM3AtimesV3p(Vec3 *r, Mat3 *m, Vec3 *v, int perm)
{
]])
LVZ("t")
LV("r","*r")
for i=0,2 do
  --LCP("v"..i,"v->c"..i,"perm")
  LC("v"..i,"v->c"..i)
  --LV("m"..i,"m->v"..i)
  LMC("m"..i,"m",i)
  for j=0,2 do
    MEQTIMES("t"..j.."r", "v"..i.."r", "m"..i..j.."r")
    PEQTIMES("t"..j.."i", "v"..i.."r", "m"..i..j.."i")
  end
  for j=0,2 do
    MEQTIMES("t"..j.."r", "v"..i.."i", "m"..i..j.."i")
    MEQTIMES("t"..j.."i", "v"..i.."i", "m"..i..j.."r")
  end
end
AVP("r","t","perm")
SV("*r","r")
printf("}\n")

--[=[
-- r[a,b,c] = m[b,d] v[a,d,c]  (a:ldv,b:nc,c:tdv)
printf([[
static inline void
PeqMtimesX(vec *r, vec *m, vec *v, int ldv, int nc, int tdv)
{
]])
LV("r","*r")
for i=0,2 do
  LC("v"..i,"v->c"..i)
  --LV("m"..i,"m->v"..i)
  LMC("m"..i,"m",i)
  for j=0,2 do
    MEQTIMES("r"..j.."r", "v"..i.."r", "m"..i..j.."r")
    PEQTIMES("r"..j.."i", "v"..i.."r", "m"..i..j.."i")
  end
  for j=0,2 do
    MEQTIMES("r"..j.."r", "v"..i.."i", "m"..i..j.."i")
    MEQTIMES("r"..j.."i", "v"..i.."i", "m"..i..j.."r")
  end
end
SV("*r","r")
printf("}\n")

printf([[
static inline void
PeqMtimesXp(vec *r, vec *m, vec *v, int perm)
{
]])
LVZ("t")
LV("r","*r")
for i=0,2 do
  --LCP("v"..i,"v->c"..i,"perm")
  LC("v"..i,"v->c"..i)
  --LV("m"..i,"m->v"..i)
  LMC("m"..i,"m",i)
  for j=0,2 do
    MEQTIMES("t"..j.."r", "v"..i.."r", "m"..i..j.."r")
    PEQTIMES("t"..j.."i", "v"..i.."r", "m"..i..j.."i")
  end
  for j=0,2 do
    MEQTIMES("t"..j.."r", "v"..i.."i", "m"..i..j.."i")
    MEQTIMES("t"..j.."i", "v"..i.."i", "m"..i..j.."r")
  end
end
AVP("r","t","perm")
SV("*r","r")
printf("}\n")
--]=]
