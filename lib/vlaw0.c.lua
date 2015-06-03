require 'Complex'
require 'Matrix'
require 'Gamma'
function printf(...) io.write(string.format(...)) end

function DC(v)
  printf("  vecr %sr;\n", v)
  printf("  vecr %si;\n", v)
end

function LC(v, e)
  printf("  vecr %sr = LD(&((%s).r));\n", v, e)
  printf("  vecr %si = LD(&((%s).i));\n", v, e)
end

function SC(e, v)
  printf("  ST(&((%s).r), %sr);\n", e, v)
  printf("  ST(&((%s).i), %si);\n", e, v)
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
  --printf("  real zero = 0;\n")
  for i=0,2 do
    --printf("  vecr %s%ir = splat(zero);\n", v, i, e, i)
    --printf("  vecr %s%ii = splat(zero);\n", v, i, e, i)
    printf("  vecr %s%ir = splat(0);\n", v, i, e, i)
    printf("  vecr %s%ii = splat(0);\n", v, i, e, i)
  end
end

function SV(v, e)
  for i=0,2 do
    printf("  ST(&((%s).c%i.r), %s%ir);\n", v, i, e, i)
    printf("  ST(&((%s).c%i.i), %s%ii);\n", v, i, e, i)
  end
end

function DH(v, e)
  for i=0,1 do
    for j=0,2 do
      printf("  vecr %s%i%ir;\n", v, i, j)
      printf("  vecr %s%i%ii;\n", v, i, j)
    end
  end
end

function LHZ(v)
  --printf("  real zero = 0;\n")
  for i=0,1 do
    for j=0,2 do
      --printf("  vecr %s%i%ir = splat(zero);\n", v, i, j)
      --printf("  vecr %s%i%ii = splat(zero);\n", v, i, j)
      printf("  vecr %s%i%ir = splat(0);\n", v, i, j)
      printf("  vecr %s%i%ii = splat(0);\n", v, i, j)
    end
  end
end

function LH(v, e)
  for i=0,1 do
    for j=0,2 do
      printf("  vecr %s%i%ir = LD(&((%s).s%i.c%i.r));\n", v, i, j, e, i, j)
      printf("  vecr %s%i%ii = LD(&((%s).s%i.c%i.i));\n", v, i, j, e, i, j)
    end
  end
end

function LFZ(v)
  --printf("  real zero = 0;\n")
  for i=0,3 do
    for j=0,2 do
      --printf("  vecr %s%i%ir = splat(zero);\n", v, i, j)
      --printf("  vecr %s%i%ii = splat(zero);\n", v, i, j)
      printf("  vecr %s%i%ir = splat(0);\n", v, i, j)
      printf("  vecr %s%i%ii = splat(0);\n", v, i, j)
    end
  end
end

function LF(v, e)
  for i=0,3 do
    for j=0,2 do
      printf("  vecr %s%i%ir = LD(&((%s).s%i.c%i.r));\n", v, i, j, e, i, j)
      printf("  vecr %s%i%ii = LD(&((%s).s%i.c%i.i));\n", v, i, j, e, i, j)
    end
  end
end

function LFP(v, e, p)
  for i=0,3 do
    for j=0,2 do
      printf("  vecr %s%i%ir = permN(LD(&((%s).s%i.c%i.r)),%s);\n", v, i, j, e, i, j, p)
      printf("  vecr %s%i%ii = permN(LD(&((%s).s%i.c%i.i)),%s);\n", v, i, j, e, i, j, p)
    end
  end
end

function SH(v, e)
  for i=0,1 do
    for j=0,2 do
      printf("  ST(&((%s).s%i.c%i.r), %s%i%ir);\n", v, i, j, e, i, j)
      printf("  ST(&((%s).s%i.c%i.i), %s%i%ii);\n", v, i, j, e, i, j)
    end
  end
end

function SF(v, e)
  for i=0,3 do
    for j=0,2 do
      printf("  ST(&((%s).s%i.c%i.r), %s%i%ir);\n", v, i, j, e, i, j)
      printf("  ST(&((%s).s%i.c%i.i), %s%i%ii);\n", v, i, j, e, i, j)
    end
  end
end

function LM(v, e)
  for i=0,2 do
    for j=0,2 do
      printf("  vecr %s%i%ir = LD(&((%s).v%i.c%i.r));\n", v, j, i, e, i, j)
      printf("  vecr %s%i%ii = LD(&((%s).v%i.c%i.i));\n", v, j, i, e, i, j)
    end
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

function PEQTIMES(x, y, z)
  printf("  %s = MADD(%s,%s,%s);\n", x, y, z, x)
end

function MEQTIMES(x, y, z)
  printf("  %s = NMADD(%s,%s,%s);\n", x, y, z, x)
end

function CEQTIMES(x, y, z)
  printf("  %sr =  MUL(%sr,%sr);\n", x, y, z)
  printf("  %si =  MUL(%si,%sr);\n", x, y, z)
  printf("  %sr = NMADD(%si,%si,%sr);\n", x, y, z, x)
  printf("  %si =  MADD(%sr,%si,%si);\n", x, y, z, x)
end

function CPEQTIMES(x, y, z)
  printf("  %sr =  MADD(%sr,%sr,%sr);\n", x, y, z, x)
  printf("  %si =  MADD(%si,%sr,%si);\n", x, y, z, x)
  printf("  %sr = NMADD(%si,%si,%sr);\n", x, y, z, x)
  printf("  %si =  MADD(%sr,%si,%si);\n", x, y, z, x)
end

function peqCxV(r, c, v)
  for i=0,2 do
    PEQTIMES(r..i..'r', c..'r', v..i..'r')
    PEQTIMES(r..i..'i', c..'r', v..i..'i')
  end
  for i=0,2 do
    MEQTIMES(r..i..'r', c..'i', v..i..'i')
    PEQTIMES(r..i..'i', c..'i', v..i..'r')
  end
end

function HeqMtimesH(r, m, h)
  for j=0,2 do
    for i=0,2 do
      for k=0,1 do
	if j==0 then
	  CEQTIMES(r..k..i, m..i..j, h..k..j)
	else
	  CPEQTIMES(r..k..i, m..i..j, h..k..j)
	end
      end
    end
  end
end

function HpeqMtimesH(r, m, h)
  for j=0,2 do
    for i=0,2 do
      for k=0,1 do
	CPEQTIMES(r..k..i, m..i..j, h..k..j)
      end
    end
  end
end

function z4(p,z)
  if ceq(p,1) then
    return {{sign='+',val=z..'r'},{sign='+',val=z..'i'}}
  elseif ceq(p,-1) then
    return {{sign='-',val=z..'r'},{sign='-',val=z..'i'}}
  elseif ceq(p,I) then
    return {{sign='-',val=z..'i'},{sign='+',val=z..'r'}}
  end
  return {{sign='+',val=z..'i'},{sign='-',val=z..'r'}}
end

function axPhase(r, p1, x1)
  local y1 = z4(p1,x1)
  --printf("%s %s %s\n", p1, y1[1].sign, y1[1].val)
  for i=1,2 do
    local rri = r..({'r','i'})[i]
    if y1[i].sign == '+' then
      printf("  %s = %s;\n", rri, y1[i].val)
    else
      printf("  %s = NEG(%s);\n", rri, y1[i].val)
    end
  end
end

function peqaxPhase(r, p1, x1)
  local y1 = z4(p1,x1)
  --printf("%s %s %s\n", p1, y1[1].sign, y1[1].val)
  for i=1,2 do
    local rri = r..({'r','i'})[i]
    if y1[i].sign == '+' then
      printf("  %s = ADD(%s,%s);\n", rri, rri, y1[i].val)
    else
      printf("  %s = SUB(%s,%s);\n", rri, rri, y1[i].val)
    end
  end
end

function axbyPhase(r, p1, x1, p2, x2)
  local y1 = z4(p1,x1)
  local y2 = z4(p2,x2)
  --printf("%s %s %s\n", p1, y1[1].sign, y1[1].val)
  for i=1,2 do
    local rri = r..({'r','i'})[i]
    if y1[i].sign == '+' then
      if y2[i].sign == '+' then
	printf("  %s = ADD(%s,%s);\n", rri, y1[i].val, y2[i].val)
      else
	printf("  %s = SUB(%s,%s);\n", rri, y1[i].val, y2[i].val)
      end
    else
      if y2[i].sign == '+' then
	printf("  %s = SUB(%s,%s);\n", rri, y2[i].val, y1[i].val)
      else
	printf("  %s = NEG(ADD(%s,%s));\n", rri, y1[i].val, y2[i].val)
      end
    end
  end
end

function sprojDir(h,f,d,s)
  local p = sprojMat[d][s]
  for i=0,1 do
    local ks = {}
    for k=1,4 do
      if not ceq(p(i+1,k),0) then ks[#ks+1] = k end
      --printf("%i %i %s\n", i+1, k, p(i+1,k))
    end
    --printf("%i %i\n", ks[1], ks[2])
    for j=0,2 do
      if #ks == 2 then
	axbyPhase(h..i..j, p(i+1,ks[1]), f..(ks[1]-1)..j,
		  p(i+1,ks[2]), f..(ks[2]-1)..j)
      else -- ks==1
	axPhase(h..i..j, p(i+1,ks[1]), f..(ks[1]-1)..j)
      end
    end
  end
end

function LsprojDir(h,f,e,i,j,d,s)
  local p = sprojMat[d][s]
  local ks = {}
  for k=1,4 do
    if not ceq(p(i+1,k),0) then ks[#ks+1] = k end
    --printf("%i %i %s\n", i+1, k, p(i+1,k))
  end
  --printf("%i %i\n", ks[1], ks[2])
  DC(h..i..j)
  local f1 = f..(ks[1]-1)..j
  LC(f1, string.format("(%s).s%i.c%i", e, ks[1]-1, j))
  if #ks == 2 then
    local f2 = f..(ks[2]-1)..j
    LC(f2, string.format("(%s).s%i.c%i", e, ks[2]-1, j))
    axbyPhase(h..i..j, p(i+1,ks[1]), f1, p(i+1,ks[2]), f2)
  else -- ks==1
    axPhase(h..i..j, p(i+1,ks[1]), f1)
  end
end

function LsprojDirP(h,f,e,i,j,d,s,perm)
  local p = sprojMat[d][s]
  local ks = {}
  for k=1,4 do
    if not ceq(p(i+1,k),0) then ks[#ks+1] = k end
    --printf("%i %i %s\n", i+1, k, p(i+1,k))
  end
  --printf("%i %i\n", ks[1], ks[2])
  DC(h..i..j)
  local f1 = f..(ks[1]-1)..j
  LCP(f1, string.format("(%s).s%i.c%i", e, ks[1]-1, j), perm)
  if #ks == 2 then
    local f2 = f..(ks[2]-1)..j
    LCP(f2, string.format("(%s).s%i.c%i", e, ks[2]-1, j), perm)
    axbyPhase(h..i..j, p(i+1,ks[1]), f1, p(i+1,ks[2]), f2)
  else -- ks==1
    axPhase(h..i..j, p(i+1,ks[1]), f1)
  end
end

function peqsreconDir(r,h,d,s)
  local p = sreconMat[d][s]
  for i=0,3 do
    local ks = {}
    for k=1,2 do
      if not ceq(p(i+1,k),0) then ks[#ks+1] = k end
      --printf("%i %i %s\n", i+1, k, p(i+1,k))
    end
    --printf("%i %i\n", ks[1], ks[2])
    if #ks > 0 then
      for j=0,2 do
	peqaxPhase(r..i..j, p(i+1,ks[1]), h..(ks[1]-1)..j)
      end
    end
  end
end

function meqsreconDir(r,h,d,s)
  local p = sreconMat[d][s]
  for i=0,3 do
    local ks = {}
    for k=1,2 do
      if not ceq(p(i+1,k),0) then ks[#ks+1] = k end
      --printf("%i %i %s\n", i+1, k, p(i+1,k))
    end
    --printf("%i %i\n", ks[1], ks[2])
    if #ks > 0 then
      for j=0,2 do
	peqaxPhase(r..i..j, -p(i+1,ks[1]), h..(ks[1]-1)..j)
      end
    end
  end
end

function LSmeqsreconDir(r,e,h,k,d,s)
  local p = sreconMat[d][s]
  for i=0,3 do
    if not ceq(p(i+1,k+1),0) then
      for j=0,2 do
	LC(r..i..j, string.format("(%s).s%i.c%i", e, i, j))
	peqaxPhase(r..i..j, -p(i+1,k+1), h..k..j)
	SC(string.format("(%s).s%i.c%i", e, i, j), r..i..j)
      end
    end
  end
end

function SPROJ(h,f,dir,sign)
  printf("  switch(%s*(%s+1)) {\n", sign, dir)
  for d=0,4 do
    for s = 0,1 do
      printf("  case %i: {\n", (1-2*s)*(d+1))
      sprojDir(h,f,d,s)
      printf("  }; break;\n")
    end
  end
  printf("  }\n")
end

function peqSRECON(r,h,dir,sign)
  printf("  switch(%s*(%s+1)) {\n", sign, dir)
  for d=0,4 do
    for s = 0,1 do
      printf("  case %i: {\n", (1-2*s)*(d+1))
      peqsreconDir(r,h,d,s)
      printf("  }; break;\n")
    end
  end
  printf("  }\n")
end

function meqSRECON(r,h,dir,sign)
  printf("  switch(%s*(%s+1)) {\n", sign, dir)
  for d=0,4 do
    for s = 0,1 do
      printf("  case %i: {\n", (1-2*s)*(d+1))
      meqsreconDir(r,h,d,s)
      printf("  }; break;\n")
    end
  end
  printf("  }\n")
end

------ start output ------

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
static inline void
sprojD3(HVec3 *r, DVec3 *f, int dir, int sign)
{
]])
LF('f','*f')
LHZ('h')
SPROJ('h','f','dir','sign')
SH('*r','h')
printf("}\n")

printf([[
static inline void
sreconD3(DVec3 *r, HVec3 *h, int dir, int sign)
{
]])
LH('h','*h')
LFZ('f')
peqSRECON('f','h','dir','sign')
SF('*r','f')
printf("}\n")

--[=[
printf([[
static inline void
PeqSprojM3timesD3(DVec3 *r, Mat3 *m, DVec3 *f, int dir, int sign)
{
]])
LF('f','*f')
LHZ('h')
SPROJ('h','f','dir','sign')
--DH('h2')
LHZ('h2')
LM('m','*m')
--HeqMtimesH('h2','m','h')
HpeqMtimesH('h2','m','h')
LF('r','*r')
peqSRECON('r','h2','dir','sign')
SF("*r","r")
printf("}\n")

printf([[
static inline void
PeqSprojM3timesD3p(DVec3 *r, Mat3 *m, DVec3 *f, int dir, int sign, int perm)
{
]])
LFP('f','*f','perm')
LHZ('h')
SPROJ('h','f','dir','sign')
--DH('h2')
LHZ('h2')
LM('m','*m')
--HeqMtimesH('h2','m','h')
HpeqMtimesH('h2','m','h')
LF('r','*r')
peqSRECON('r','h2','dir','sign')
SF("*r","r")
printf("}\n")
--]=]

printf([[
static inline void
MeqSprojM3timesD3(DVec3 *r, Mat3 *m, DVec3 *f, int dir, int sign)
{
]])
LF('f','*f')
LHZ('h')
SPROJ('h','f','dir','sign')
--DH('h2')
LHZ('h2')
LM('m','*m')
--HeqMtimesH('h2','m','h')
HpeqMtimesH('h2','m','h')
LF('r','*r')
meqSRECON('r','h2','dir','sign')
SF("*r","r")
printf("}\n")

printf([[
static inline void
MeqSprojM3timesD3p(DVec3 *r, Mat3 *m, DVec3 *f, int dir, int sign, int perm)
{
]])
LFP('f','*f','perm')
LHZ('h')
SPROJ('h','f','dir','sign')
--DH('h2')
LHZ('h2')
LM('m','*m')
--HeqMtimesH('h2','m','h')
HpeqMtimesH('h2','m','h')
LF('r','*r')
meqSRECON('r','h2','dir','sign')
SF("*r","r")
printf("}\n")
--]=]

for dir=0,4 do
  for s=0,1 do
    local sgn = ({'p','m'})[s+1]

printf([[
static inline void
MeqSproj%i%sM3timesD3(DVec3 *r, Mat3 *m, DVec3 *f)
{
]],dir,sgn)
for k=0,1 do
  LVZ('t'..k)
  for j=0,2 do
    LsprojDir('h','f','*f',k,j,dir,s)
    if k==0 then
      LV('m'..j, 'm->v'..j)
    end
    peqCxV('t'..k, 'h'..k..j, 'm'..j)
  end
  LSmeqsreconDir('r','*r','t',k,dir,s)
end
printf("}\n")

--[=[
printf([[
static inline void
MeqSproj%i%sM3timesD3a(DVec3 *r, Mat3 *m, DVec3 *f)
{
]],dir,sgn)
LF('f','*f')
DH('h')
sprojDir('h','f',dir,s)
LM('m','*m')
DH('h2')
--LHZ('h2')
HeqMtimesH('h2','m','h')
--HpeqMtimesH('h2','m','h')
LF('r','*r')
meqsreconDir('r','h2',dir,s)
SF("*r","r")
printf("}\n")
--]=]

printf([[
static inline void
MeqSproj%i%sM3timesD3p(DVec3 *r, Mat3 *m, DVec3 *f, int perm)
{
]],dir,sgn)
for k=0,1 do
  LVZ('t'..k)
  for j=0,2 do
    LsprojDirP('h','f','*f',k,j,dir,s,'perm')
    if k==0 then
      LV('m'..j, 'm->v'..j)
    end
    peqCxV('t'..k, 'h'..k..j, 'm'..j)
  end
  LSmeqsreconDir('r','*r','t',k,dir,s)
end
printf("}\n")

--[=[
printf([[
static inline void
MeqSproj%i%sM3timesD3p(DVec3 *r, Mat3 *m, DVec3 *f, int perm)
{
]],dir,sgn)
LFP('f','*f','perm')
DH('h')
sprojDir('h','f',dir,s)
DH('h2')
--LHZ('h2')
LM('m','*m')
HeqMtimesH('h2','m','h')
--HpeqMtimesH('h2','m','h')
LF('r','*r')
meqsreconDir('r','h2',dir,s)
SF("*r","r")
printf("}\n")
--]=]

printf([[
void
P(D_meq_spinproj%i%s_M_times_xDp)(real *rr, real *mm, int *idx, int perm,
			          real *vv, int i0, int i1)
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
      MeqSproj%i%sM3timesD3(ri, mi, vi);
    } else if(k1+2<=0) {
      int k2 = -(k1+2);
      DVec3 *restrict ri = &r[i];
      Mat3 *restrict mi = &m[i];
      DVec3 *restrict vi = &v[k2];
      MeqSproj%i%sM3timesD3p(ri, mi, vi, perm);
    }
  }
}
]], dir, sgn, dir, sgn, dir, sgn)

  end
end
