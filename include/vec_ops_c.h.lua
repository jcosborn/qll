function printf(...) io.write(string.format(...)) end

function rtype(p)
  if p=='f' then return 'float' end
  return 'double'
end
function vtype(n,p)
  return "vec"..n..p
end
function vtyper(n,p)
  return "vec"..n..p.."r"
end
function vfunc(n,p,f)
  return "v"..n..p..f
end
function V(s,i,n,p)
  --return "((("..rtype(p).."*)&("..s.."))["..i.."])"
  local e = vfunc(n,p,'elem')
  return e.."("..s..","..i..")"
end

local band, bxor, srep
if bit32 then
  band = bit32.band
  bxor = bit32.bxor
else
  band = function(a,b)
    --printf("band(%g,%g)", a, b)
    local c,n = 0,1
    while a>0 or b>0 do
      local ta = (a%2)==1
      local tb = (b%2)==1
      local t = ta and tb
      if t then c = c + n end
      n = 2*n
      a = math.floor(a/2)
      b = math.floor(b/2)
    end
    --printf(" = %g\n", c)
    return c
  end
  bxor = function(a,b)
    --printf("bxor(%g,%g)", a, b)
    local c,n = 0,1
    while a>0 or b>0 do
      local ta = (a%2)==1
      local tb = (b%2)==1
      local t = ta == (not tb)
      if t then c = c + n end
      n = 2*n
      a = math.floor(a/2)
      b = math.floor(b/2)
    end
    --printf(" = %g\n", c)
    return c
  end
end
if string.rep('a',2,'b') == 'aba' then
  srep = string.rep
else
  srep = function(s,n,p)
    if n>0 then
      return s..string.rep(p..s, n-1)
    end
  end
end

--[[
 *ld,*st,*splat,*reduce
 *neg,*add,*sub,*mul,*madd,*nmadd
 *perm1,2,4,8
 *packp1,*packm1,2,4,8
 *blendp1,*blendm1,2,4,8
 v8fldv4f(a)
 v8fldv4f2(a,b)
 v4fstv8f(a,x)
 v4f2stv8f(a,b,x)
--]]

function ld(n,p)
  local tm = vtype(n,p)
  local t = vtyper(n,p)
  local f = vfunc(n,p,"ld")
  printf("#ifndef %s\n", f)
  printf("#define %s(x) (*(%s*)(x))\n", f, t)
--[[  printf("#define %s(x) %sf(x)\n", f, f)
  local function e(i) return V('*x',i,p) end
  printf("static inline %s %sf(%s *x)\n{\n  return (%s){%s", t, f, tm, t, e(0))
  for i=1,n-1 do printf(","..e(i)) end
  printf("};\n}\n")
--]]
  printf("#endif\n")
end

function st(n,p)
  local t = vtyper(n,p)
  local f = vfunc(n,p,"st")
  printf("#ifndef %s\n", f)
  printf("#define %s(x,y) (*(%s*)(x)) = (y)\n", f, t)
  printf("#endif\n")
end

function stnr(n,p)
  local t = vtyper(n,p)
  local f = vfunc(n,p,"stnr")
  printf("#ifndef %s\n", f)
  printf("#define %s(x,y) (*(%s*)(x)) = (y)\n", f, t)
  printf("#endif\n")
end

function splat(n,p)
  local t = vtyper(n,p)
  local f = vfunc(n,p,"splat")
  local init = vfunc(n,p,"init")
  --local a = string.rep("(x)", n, ",")
  local a = srep("x", n, ",")
  printf("#ifndef %s\n", f)
  --printf("#define %s(x) (%s){%s}\n", f, t, a)
  printf("#define %s(x) %s(%s)\n", f, init, a)
  printf("#endif\n")
end

function reduce(n,p)
  local t = vtyper(n,p)
  local f = vfunc(n,p,"reduce")
  local function e(i) return V('x',i,n,p) end
  printf("#ifndef %s\n", f)
  printf("#define %s(x) (%s", f, e(0))
  for i=1,n-1 do printf("+"..e(i)) end
  printf(")\n")
  printf("#endif\n")
end

function neg(n,p)
  local t = vtyper(n,p)
  local f = vfunc(n,p,"neg")
  local init = vfunc(n,p,"init")
  local function e(i) return "-"..V('x',i,n,p) end
  printf("#ifndef %s\n", f)
  printf("#define %s(x) %s(%s", f, init, e(0))
  for i=1,n-1 do printf(","..e(i)) end
  printf(")\n")
  printf("#endif\n")
end

function add(n,p)
  local t = vtyper(n,p)
  local f = vfunc(n,p,"add")
  local init = vfunc(n,p,"init")
  local function e(i) return V('x',i,n,p).."+"..V('y',i,n,p) end
  printf("#ifndef %s\n", f)
  printf("#define %s(x,y) %s(%s", f, init, e(0))
  for i=1,n-1 do printf(","..e(i)) end
  printf(")\n")
  printf("#endif\n")
end

function sub(n,p)
  local t = vtyper(n,p)
  local f = vfunc(n,p,"sub")
  local init = vfunc(n,p,"init")
  local function e(i) return V('x',i,n,p).."-"..V('y',i,n,p) end
  printf("#ifndef %s\n", f)
  printf("#define %s(x,y) %s(%s", f, init, e(0))
  for i=1,n-1 do printf(","..e(i)) end
  printf(")\n")
  printf("#endif\n")
end

function mul(n,p)
  local t = vtyper(n,p)
  local f = vfunc(n,p,"mul")
  local init = vfunc(n,p,"init")
  local function e(i) return V('x',i,n,p).."*"..V('y',i,n,p) end
  printf("#ifndef %s\n", f)
  printf("#define %s(x,y) %s(%s", f, init, e(0))
  for i=1,n-1 do printf(","..e(i)) end
  printf(")\n")
  printf("#endif\n")
end

function madd(n,p)
  local t = vtyper(n,p)
  local f = vfunc(n,p,"madd")
  local init = vfunc(n,p,"init")
  local function e(i)
    return '('..V('z',i,n,p)..'+'..V('x',i,n,p).."*"..V('y',i,n,p)..')'
  end
  printf("#ifndef %s\n", f)
--[[
  printf("#define %s(x,y,z) %s(%s", f, init, e(0))
  for i=1,n-1 do printf(","..e(i)) end
  printf(")\n")
--]]
  printf("#define %s(x,y,z) %sf(x,y,z)\n", f, f)
  printf("static inline %s %sf(%s x, %s y, %s z)\n{\n", t, f, t, t, t)
  printf("  %s r;\n", t)
  for i=0,n-1 do
    printf("  %s = %s * %s + %s;\n", V('r',i,n,p), V('x',i,n,p), V('y',i,n,p), V('z',i,n,p))
  end
  printf("  return r;\n}\n")
  printf("#endif\n")
end

function nmadd(n,p)
  local t = vtyper(n,p)
  local f = vfunc(n,p,"nmadd")
  local init = vfunc(n,p,"init")
  local function e(i)
    return '('..V('z',i,n,p)..'-'..V('x',i,n,p).."*"..V('y',i,n,p)..')'
  end
  printf("#ifndef %s\n", f)
  printf("#define %s(x,y,z) %s(%s", f, init, e(0))
  for i=1,n-1 do printf(","..e(i)) end
  printf(")\n")
  printf("#endif\n")
end

function rsqrt(n,p)
  local t = vtyper(n,p)
  local f = vfunc(n,p,"rsqrt")
  printf("#ifndef %s\n", f)
  printf("#define %s(x) %sf(x)\n", f, f)
  printf("static inline %s %sf(%s x)\n{\n", t, f, t)
  printf("  %s r;\n", t)
  for i=0,n-1 do
    printf("  %s = 1/sqrt(%s);\n", V('r',i,n,p), V('x',i,n,p))
  end
  printf("  return r;\n}\n")
  printf("#endif\n")
end

function perm(n,p,k)
  local t = vtyper(n,p)
  local f = vfunc(n,p,"perm"..k)
  local init = vfunc(n,p,"init")
  local function e(i) return V('x',bxor(i,k)%n,n,p) end
  printf("#ifndef %s\n", f)
  printf("#define %s(x) %sf(x)\n", f, f)
  printf("static inline %s %sf(%s x)\n{\n  return %s(%s", t, f, t, init, e(0))
  for i=1,n-1 do printf(","..e(i)) end
  printf(");\n}\n")
  printf("#endif\n")
end

function pack(n,p,s,k)
  local tp = rtype(p)
  local t = vtyper(n,p)
  local f = vfunc(n,p,"pack"..(({'m','p'})[s+1])..k)
  printf("#ifndef %s\n", f)
  printf("#define %s(x,y) %sf(x,y)\n", f, f)
  printf("static inline void %sf(%s *x, %s y)\n{\n", f, tp, t)
  local j = 0
  for i=0,n-1 do
    if band(i,k) == s*k then
      printf("  x[%i] = %s;\n", j, V('y',i,n,p))
      j = j + 1
    end
  end
  printf("}\n")
  printf("#endif\n")
end

function blend(n,p,s,k)
  local tp = rtype(p)
  local t = vtyper(n,p)
  local f = vfunc(n,p,"blend"..(({'m','p'})[s+1])..k)
  local init = vfunc(n,p,"init")
  local function e(i,j)
    if band(i,k) == s*k then
      return V('x',bxor(i,k)%n,n,p), j
    else
      return 'y['..j..']', j+1
    end
  end
  printf("#ifndef %s\n", f)
  printf("#define %s(x,y) %sf(x,y)\n", f, f)
  printf("static inline %s %sf(%s x, %s *y)\n{\n", t, f, t, tp)
  local x,j = e(0,0)
  printf("  return %s(%s", init, x)
  for i=1,n-1 do
    x,j = e(i,j)
    printf(",%s", x)
  end
  printf(");\n}\n")
  printf("#endif\n")
end

-- begin output

printf('#include "math.h"\n')

for k=0,5 do
  local n = 2^k
  for ip,p in ipairs({'f','d'}) do
    printf("\n// vector length: %i  precision: %s\n\n", n, p)
    ld(n,p)
    st(n,p)
    stnr(n,p)
    splat(n,p)
    reduce(n,p)
    neg(n,p)
    add(n,p)
    sub(n,p)
    mul(n,p)
    madd(n,p)
    nmadd(n,p)
    rsqrt(n,p)
    for i=0,k-1 do
      local j = 2^i
      perm(n,p,j)
    end
    for i=0,k-1 do
      local j = 2^i
      pack(n,p,1,j)
      pack(n,p,0,j)
    end
    for i=0,k-1 do
      local j = 2^i
      blend(n,p,1,j)
      blend(n,p,0,j)
    end
  end
end
