local MT = {}
MT.__index = MT

function ComplexRI(x,y)
  return setmetatable({r=x,i=y},MT)
end

function isComplex(x)
  return getmetatable(x)==MT
end

function Complex(x,y)
  if isComplex(x) then
    return ComplexRI(x.r,x.i)
  else
    return ComplexRI(x,y or 0)
  end
end

function ceq(x,y)
  local xr,xi,yr,yi
  if isComplex(x) then
    xr = x.r
    xi = x.i
  else
    xr = x
    xi = 0
  end
  if isComplex(y) then
    yr = y.r
    yi = y.i
  else
    yr = y
    yi = 0
  end
  --print(xr,xi,yr,yi)
  return (xr == yr) and (xi == yi)
end

function cabs2(x)
  if isComplex(x) then return x.r*x.r+x.i*x.i end
  return x*x
end

function MT:__tostring()
  local r = ''
  local needSign = false
  if self.r~=0 then
    r = string.format("%g", self.r)
    needSign = true
  end
  if self.i==0 then
    if not needSign then r = '0' end
  else
    if self.i>0 then
      if needSign then r = r..' + ' end
    else
      if needSign then r = r..' - '
      else r = '-' end
    end
    local a = math.abs(self.i)
    if a==1 then r = r..'i'
    else r = r..string.format("%gi", a) end
  end
  return r
end

function MT:__unm()
  return ComplexRI(-self.r, -self.i)
end

function MT:__add(x)
  if isComplex(self) then
    if isComplex(x) then
      return ComplexRI(self.r+x.r, self.i+x.i)
    else
      return ComplexRI(self.r+x, self.i)
    end
  else
    return ComplexRI(self+x.r, x.i)
  end
end

function MT:__sub(x)
  if isComplex(self) then
    if isComplex(x) then
      return ComplexRI(self.r-x.r, self.i-x.i)
    else
      return ComplexRI(self.r-x, self.i)
    end
  else
    return ComplexRI(self-x.r, -x.i)
  end
end

function MT:__mul(x)
  if isComplex(self) then
    if isComplex(x) then
      return ComplexRI(self.r*x.r-self.i*x.i, self.r*x.i+self.i*x.r)
    else
      return ComplexRI(self.r*x, self.i*x)
    end
  else
    return ComplexRI(self*x.r, self*x.i)
  end
end

function MT:conj()
  return ComplexRI(self.r, -self.i)
end

function MT:__eq(x)
  local sr = self.r or self
  local xr = x.r or x
  if sr ~= xr then return false end
  local si = self.i or 0
  local xi = x.i or 0
  return si == xi
end

--[[
print(Complex(1,0))
print(Complex(0,1))
print(Complex(-1,0))
print(Complex(0,-1))
print(Complex(1,1))
print(Complex(1,-1))
print(Complex(-1,1))
print(Complex(-1,-1))
print(Complex(2,2))
print(Complex(2,-2))
print(Complex(-2,2))
print(Complex(-2,-2))
--]]
