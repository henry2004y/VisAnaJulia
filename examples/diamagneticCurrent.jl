using PyCall, PyPlot#, ScatteredInterpolation

function get_diamagnetic_current(filehead::Dict, data::Array{T,2}) where
   T<:AbstractFloat

   @views X = data[1,:]
   @views Y = data[2,:]
   @views Z = data[3,:]
   @views Bx0 = data[11,:]
   @views By0 = data[12,:]
   @views Bz0 = data[13,:]
   @views P0 = data[14,:]

   points = hcat(X,Y,Z);
   #points = [X Y Z]' # for ScatteredInterpolation

   np = pyimport("numpy")
   scipy = pyimport("scipy")
   interpolate = pyimport("scipy.interpolate")
   griddata = interpolate.griddata

   #=
   xMin = minimum(X); xMax = maximum(X)
   yMin = minimum(Y); yMax = maximum(Y)
   zMin = minimum(Z); zMax = maximum(Z)
   =#

   xMin, xMax = -3.0, 3.0
   yMin, yMax = -3.0, 3.0
   zMin, zMax = -3.0, 3.0

   nX, nY, nZ = 100, 3, 100

   s = pybuiltin(:slice)
   XYZ = get(np.mgrid,
      (s(xMin,xMax,nX*1im), s(yMin,yMax,nY*1im), s(zMin,zMax,nZ*1im)))
   grid_x, grid_y, grid_z = XYZ[1,:,:,:], XYZ[2,:,:,:], XYZ[3,:,:,:]

   dx = (xMax - xMin) / nX

   #= # The major issue with ScatteredInterpolation is memory!
   XYZ = reshape(XYZ,3,nX*nY*nZ)

   # Let me try ScatteredInterpolation in Julia
   itp = interpolate(Shepard(), points, Bx0)
   Bx = evaluate(itp, XYZ)
   Bx = reshape(Bx,nX,nY,nZ)

   itp = interpolate(Shepard(), points, By0)
   By = evaluate(itp, XYZ)
   By = reshape(By,nX,nY,nZ)

   itp = interpolate(Shepard(), points, Bz0)
   Bz = evaluate(itp, XYZ)
   Bz = reshape(Bz,nX,nY,nZ)

   itp = interpolate(Shepard(), points, P0)
   P = evaluate(itp, XYZ)
   P = reshape(P,nX,nY,nZ)
   =#

   Bx = griddata(points, Bx0, (grid_x, grid_y, grid_z), method="linear")
   By = griddata(points, By0, (grid_x, grid_y, grid_z), method="linear")
   Bz = griddata(points, Bz0, (grid_x, grid_y, grid_z), method="linear")
   P  = griddata(points, P0,  (grid_x, grid_y, grid_z), method="linear")

   px, py, pz = np.gradient(P, dx) # ▽P, unit length should be provided!!!

   jDiaMag  = Array{Float64,4}(undef,3,nX,nY,nZ)
   jCurrent = Array{Float64,3}(undef,nX,nY,nZ)

   # j = B×▽p/B^2
   @inbounds for k = 1:nZ, j = 1:nY, i = 1:nX
      Bmag = hypot(Bx[i,j,k], By[i,j,k], Bz[i,j,k])
      jDiaMag[1,i,j,k] = (By[i,j,k]*pz[i,j,k] - Bz[i,j,k]*py[i,j,k]) / Bmag
      jDiaMag[2,i,j,k] = (Bz[i,j,k]*px[i,j,k] - Bx[i,j,k]*pz[i,j,k]) / Bmag
      jDiaMag[3,i,j,k] = (Bx[i,j,k]*py[i,j,k] - By[i,j,k]*px[i,j,k]) / Bmag
      jCurrent[i,j,k] = hypot(jDiaMag[:,i,j,k]...)
   end

   x = range(xMin, stop=xMax, length=nX)
   z = range(zMin, stop=zMax, length=nZ)

   contourf(x,z,jCurrent[:,2,:]',50)
   colorbar()
   xlabel("x")
   ylabel("z")
end

#get_diamagnetic_current(data.head., data)

using Glob

const q = 1.6021765e-19 #[C]
const dx = 1/32 #[Rg]
const Rg = 2634000. #[m]
const Rm = 2444000. #[m]

filename = "3d_var_region0_0_*.out"
#filename = "box*.outs"

function test(filename::String)


np = pyimport("numpy")

filedir = "/Users/hyzhou/Ganymede/MOP2018/runG8_PIC_1200s/PC";

filenames = glob(filename, filedir)

#for fname in filenames
for i = 80:100#length(filenames)
   fname = filenames[i]
   dir, name = splitdir(fname)
   data = readdata(name,dir=dir,verbose=false);
   @views x = data.x[:,:,:,1]
   @views y = data.x[:,:,:,2]
   @views z = data.x[:,:,:,3]
   w = data.w

   ne_ = findfirst(isequal("rhoS0"), data.head.wnames)
   ni_ = findfirst(isequal("rhoS1"), data.head.wnames)
   pi_ = findfirst(isequal("pS1"),   data.head.wnames)
   pe_ = findfirst(isequal("pS0"),   data.head.wnames)

   bx_ = findfirst(isequal("Bx"),    data.head.wnames)
   by_ = findfirst(isequal("By"),    data.head.wnames)
   bz_ = findfirst(isequal("Bz"),    data.head.wnames)


   #ne = w[:,:,:,ne_]
   ni = w[:,:,:,ni_] / 14 *1e6  # [/m^3]
   pi = w[:,:,:,pi_]*1e-9 # [Pa]
   Bx = w[:,:,:,bx_]*1e-9 # [T]
   By = w[:,:,:,by_]*1e-9
   Bz = w[:,:,:,bz_]*1e-9

   px, py, pz = np.gradient(pi,dx*Rg)

   nX, nY, nZ = size(w)

   vDiaMag  = Array{Float64,4}(undef,3,nX,nY,nZ)
   vMag = Array{Float64,3}(undef,nX,nY,nZ)

   # v = B×▽p/neB^2
   @inbounds for k = 1:nZ, j = 1:nY, i = 1:nX
      Bmag2 = Bx[i,j,k]^2 + Bx[i,j,k]^2 + Bz[i,j,k]^2
      vDiaMag[1,i,j,k] = (By[i,j,k]*pz[i,j,k] - Bz[i,j,k]*py[i,j,k]) /
         (ni[i,j,k]*q*Bmag2)
      vDiaMag[2,i,j,k] = (Bz[i,j,k]*px[i,j,k] - Bx[i,j,k]*pz[i,j,k]) /
         (ni[i,j,k]*q*Bmag2)
      vDiaMag[3,i,j,k] = (Bx[i,j,k]*py[i,j,k] - By[i,j,k]*px[i,j,k]) /
         (ni[i,j,k]*q*Bmag2)
      vMag[i,j,k] = hypot(vDiaMag[:,i,j,k]...)
   end

   yMid = floor(Int,nY/2)
   xrange = floor(Int,nX/4):floor(Int,nX/4*3)
   zrange = floor(Int,nZ/4):floor(Int,nZ/4*3)

   clf()


   #cf = contourf(x[:,yMid,:],z[:,yMid,:],vMag[:,yMid,:],levels=np.logspace(6,8,30))
   #cf = contourf(x[:,yMid,:],z[:,yMid,:],vMag[:,yMid,:],
   #   locator=matplotlib.ticker.LogLocator(subs=collect(1:10)))
   #cf = contourf(x[:,yMid,:],z[:,yMid,:],vMag[:,yMid,:], levels=50)
   cf = contourf(x[xrange,yMid,zrange],z[xrange,yMid,zrange],vMag[xrange,yMid,zrange], levels=50)

   plt.axis("equal")
   cbar = colorbar(cf)
   #plt.clim(0, 3e6) # This creates the same color for values above
   xlabel("x"); ylabel("z")
   title(L"$\mathbf{B}\times \nabla P/(neB^2)$"*", t=$(i)")
   plt.savefig("testoutput/diamagnetic_$(i).png")


end


end


function ganymede_plot(filename::String)

filename = "y*.out"

filedir = "/Users/hyzhou/Ganymede/MOP2018/runG8_PIC_1200s/GM";

filenames = glob(filename, filedir)

for i = 80:100#length(filenames)
   fname = filenames[i]
   dir, name = splitdir(fname)
   data = readdata(name,dir=dir,verbose=false);

   plotdata(data,"p bx;bz",plotmode="contbar streamover", plotrange=[-2.2,-1.5,-1.,1.],
   	density=1.5, plotinterval=0.05)
   plt.axis("equal")

   plt.savefig("testcut/test_$(i).png")
   close()

end


end


function test_mercury(filename::String)


np = pyimport("numpy")


filedir = "/Users/hyzhou/Ganymede/run_mercury_80s/GM"

filenames = glob(filename, filedir)

dir, name = splitdir(filenames[1])
data = readdata(name,dir=dir,verbose=false)
npict = data.list.npictinfiles

for ipict = 1:npict
   data = readdata(name,npict=ipict,dir=dir,verbose=false)
   @views x = data.x[:,:,:,1]
   @views y = data.x[:,:,:,2]
   @views z = data.x[:,:,:,3]
   w = data.w

   ni_ = findfirst(isequal("Rho"), data.head.wnames)
   pi_ = findfirst(isequal("P"),  data.head.wnames)

   bx_ = findfirst(isequal("Bx"), data.head.wnames)
   by_ = findfirst(isequal("By"), data.head.wnames)
   bz_ = findfirst(isequal("Bz"), data.head.wnames)


   #ne = w[:,:,:,ne_]
   ni = w[:,:,:,ni_] / 14 *1e6  # [/m^3]
   pi = w[:,:,:,pi_]*1e-9 # [Pa]
   Bx = w[:,:,:,bx_]*1e-9 # [T]
   By = w[:,:,:,by_]*1e-9
   Bz = w[:,:,:,bz_]*1e-9

   px, py, pz = np.gradient(pi,dx*Rg)

   nX, nY, nZ = size(w)

   vDiaMag  = Array{Float64,4}(undef,3,nX,nY,nZ)
   vMag = Array{Float64,3}(undef,nX,nY,nZ)

   # v = B×▽p/neB^2
   @inbounds for k = 1:nZ, j = 1:nY, i = 1:nX
      Bmag2 = Bx[i,j,k]^2 + Bx[i,j,k]^2 + Bz[i,j,k]^2
      vDiaMag[1,i,j,k] = (By[i,j,k]*pz[i,j,k] - Bz[i,j,k]*py[i,j,k]) /
         (ni[i,j,k]*q*Bmag2)
      vDiaMag[2,i,j,k] = (Bz[i,j,k]*px[i,j,k] - Bx[i,j,k]*pz[i,j,k]) /
         (ni[i,j,k]*q*Bmag2)
      vDiaMag[3,i,j,k] = (Bx[i,j,k]*py[i,j,k] - By[i,j,k]*px[i,j,k]) /
         (ni[i,j,k]*q*Bmag2)
      vMag[i,j,k] = hypot(vDiaMag[:,i,j,k]...)
   end

   yMid = floor(Int,nY/2)

   clf()

   cf = contourf(x[:,yMid,:],z[:,yMid,:],vMag[:,yMid,:],levels=50)
   cbar = colorbar(cf)
   plt.clim(0,8000000) # This creates the same color for values above


   xlabel("x"); ylabel("z")
   title(L"$\mathbf{B}\times \nabla P/(neB^2)$")
   plt.savefig("testoutput/test_"*string(ipict)*".png")


end


end



function test_mercuryplot(filename::String)



filedir = "/Users/hyzhou/Ganymede/run_mercury_80s/GM"

filenames = glob(filename, filedir)

dir, name = splitdir(filenames[1])
data = readdata(name,dir=dir,verbose=false)
npict = data.list.npictinfiles

for ipict = 1:npict
   data= readdata(name,npict=ipict,dir=dir,verbose=false)

   plotdata(data,"p bx;bz",plotmode="contbar streamover",cut="y")


   clf()


   plt.clim(0,8000000) # This creates the same color for values above


   xlabel("x"); ylabel("z")
   title(L"$\mathbf{B}\times \nabla P/(neB^2)$")
   plt.savefig("testoutput/test_"*string(ipict)*".png")


end


end


#=


   gradP  = Array{Float64,4}(undef,3,nX,nY,nZ)
   gradPMag = Array{Float64,3}(undef,nX,nY,nZ)
   @inbounds for k = 1:nZ, j = 1:nY, i = 1:nX
       gradP[1,i,j,k] = (px[i,j,k]) / (ni[i,j,k]*q)
       gradP[2,i,j,k] = (py[i,j,k]) / (ni[i,j,k]*q)
       gradP[3,i,j,k] = (pz[i,j,k]) / (ni[i,j,k]*q)
       gradPMag[i,j,k] = hypot(gradP[:,i,j,k]...)
   end

   figure()
   contourf(x[:,yMid,:],z[:,yMid,:],gradPMag[:,yMid,:],50)
   colorbar()

   ex_ = findfirst(isequal("Ex"),    data.head.wnames)
   ey_ = findfirst(isequal("Ey"),    data.head.wnames)
   ez_ = findfirst(isequal("Ez"),    data.head.wnames)

   Ex = @view w[:,:,:,ex_]
   Ey = @view w[:,:,:,ey_]
   Ez = @view w[:,:,:,ez_]

   E  = @. hypot(Ex, Ey, Ez)

   figure()
   contourf(x[:,yMid,:],z[:,yMid,:],E[:,yMid,:],50)
   colorbar()

=#
