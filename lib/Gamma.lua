require 'Complex'
require 'Matrix'
function printf(...) io.write(string.format(...)) end

gammaMatrix = {}
I = Complex(0,1)
gammaMatrix[0] = Matrix(4,4,{1,0,0,0, 0,1,0,0, 0,0,1,0, 0,0,0,1})
gammaMatrix[1] = Matrix(4,4,{0,0,0,I, 0,0,I,0, 0,-I,0,0, -I,0,0,0})
gammaMatrix[2] = Matrix(4,4,{0,0,0,-1, 0,0,1,0, 0,1,0,0, -1,0,0,0})
gammaMatrix[4] = Matrix(4,4,{0,0,I,0, 0,0,0,-I, -I,0,0,0, 0,I,0,0})
gammaMatrix[8] = Matrix(4,4,{0,0,1,0, 0,0,0,1, 1,0,0,0, 0,1,0,0})
for g=0,15 do
  if not gammaMatrix[g] then
    local m = gammaMatrix[0]
    local h,b = g,1
    repeat
      if h%2==1 then
        --printf("%i\t%i\t%i\n", g, h, b)
        m = m * gammaMatrix[b]
      end
      h,b = math.floor(h/2),2*b
    until h == 0
    gammaMatrix[g] = m
  end
  --myprint("g[",tostring(g),"]=",gammaelem[g],"\n")
end

gammaDir = {}
gammaDir[0] = gammaMatrix[1]
gammaDir[1] = gammaMatrix[2]
gammaDir[2] = gammaMatrix[4]
gammaDir[3] = gammaMatrix[8]
gammaDir[4] = gammaMatrix[15]

function conj(x)
  if type(x) == "number" then return x end
  return x:conj()
end

function adjMat(m)
  local r = Matrix(m.nc,m.nr)
  for i=1,m.nr do
    for j=1,m.nc do
      r(j,i, conj(m(i,j)))
    end
  end
  return r
end

local gproj = Matrix(2,4,{1,0,0,0, 0,1,0,0})
local gproj2 = Matrix(2,4,{0,0,1,0, 0,0,0,1})

sprojMat = {}
sreconMat = {}
for d=0,4 do
  sprojMat[d] = {}
  sreconMat[d] = {}
  for s=0,1 do
    local sign = 1-2*s
    local t = gammaMatrix[0]+sign*gammaDir[d]
    if d==4 then t = 0.5*t end
    if d==4 and s==1 then
      sprojMat[d][s] = gproj2*t
    else
      sprojMat[d][s] = gproj*t
    end
    if d==4 then
      sreconMat[d][s] = t * adjMat(sprojMat[d][s])
    else
      sreconMat[d][s] = 0.5*t * adjMat(sprojMat[d][s])
    end
    --printf("sprojMat[%i][%i]:\n%s\n", d, s, tostring(sprojMat[d][s]))
    --printf("sreconMat[%i][%i]:\n%s\n", d, s, tostring(sreconMat[d][s]))
  end
end
