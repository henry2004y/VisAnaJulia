# Script for further analysis of simulation results.
# Diamagnetic velocity, plasma beta
#
# Hongyang Zhou, hyzhou@umich.edu 10/01/2019

using Batsrus, VisAna, PyPlot, PyCall, Glob, Printf

include("constants.jl")

"Series of contour plots from box output."
function plot_contour_from_box(filename::String, dir=".")

   data = readdata(filename, dir=dir, verbose=false);
   npict = data.list.npictinfiles

   saveDir = "ycutsFromBox"
   mkdir(saveDir)

   for ipict = 1:npict
      data = readdata(filename, npict=ipict, dir=dir)

      plotdata(data, "by bx;bz",plotmode="contbar streamover",
         cut="y",cutPlaneIndex=3, density=1.5)

      fig = matplotlib.pyplot.gcf()
      fig.set_size_inches(3.5,6)
      plt.axis("scaled")

      plt.savefig(saveDir*"/y_"*string(ipict)*".png")
      close()
   end

end

"Series of contour plots from cuts."
function plot_contour_from_cuts(filename::String, dir=".")

   data = readdata(filename, dir=dir)
   npict = data.list.npictinfiles
   saveDir = "ycuts"
   mkdir(saveDir)

   for ipict = 1:npict
      data = readdata(filename,npict=ipict,dir=dir)

      plotdata(data, "by bx;bz", plotmode="contbar streamover",
         plotrange=[1.0 1.8 -2.0 2.0], density=1.5)

      fig = matplotlib.pyplot.gcf()
      fig.set_size_inches(3.5,6)
      plt.axis("scaled")

      plt.savefig(saveDir*"/y_"*string(ipict)*".png")
      close()
   end

end



"""
   get_diamagnetic_velocity(filename, filedir, nType)

Calculate diamagnetic velocity from outputs.

nType = 1 is Hall MHD, nType = 2 is PIC.
"""
function get_diamagnetic_velocity(filename::String, filedir=".", nType=1)

   saveDir = "ycuts"
   mkdir(saveDir)

   if filename[end-3:end] == "outs"
      data = readdata(filename, dir=filedir, verbose=false)
      npict = data.list.npictinfiles
      for i = 1:npict
         data = readdata(filename, dir=filedir, npict=i)
         processing(data, saveDir, nType)
      end
   else
      filenames = glob(filename, filedir)
      npict = length(filenames)
      for i = 1:npict
         data = readdata(filenames[i], dir=filedir)
         processing(data, saveDir, nType)
      end
   end

end

function processing(data::Data, saveDir::String, nType=1)

   np = pyimport("numpy")

   @views x = data.x[:,:,:,1]
   @views y = data.x[:,:,:,2]
   @views z = data.x[:,:,:,3]
   w = data.w

   if nType == 1 # Hall MHD
      ni_ = findfirst(isequal("Rho"), data.head.wnames)
      pi_ = findfirst(isequal("P"),   data.head.wnames)
   elseif nType == 2 # PIC
      ni_ = findfirst(isequal("rhos1"), data.head.wnames)
      pi_ = findfirst(isequal("ps1"),   data.head.wnames)
   end

   bx_ = findfirst(isequal("Bx"),  data.head.wnames)
   by_ = findfirst(isequal("By"),  data.head.wnames)
   bz_ = findfirst(isequal("Bz"),  data.head.wnames)

   ni = w[:,:,:,ni_] # [/cm^3]
   pi = w[:,:,:,pi_] # [nPa]
   Bx = w[:,:,:,bx_] # [nT]
   By = w[:,:,:,by_] # [nT]
   Bz = w[:,:,:,bz_] # [nT]

   # This assumes uniform grid resolution!
   px, py, pz = np.gradient(pi, x[2]-x[1])

   nX, nY, nZ = size(w)

   vDiaMag = Array{Float64,4}(undef,3,nX,nY,nZ)
   vMag    = Array{Float64,3}(undef,nX,nY,nZ)

   cUnit = 1.0e7/(14*1.602*2.634)

   # v = B×▽p/neB^2
   @inbounds for k = 1:nZ, j = 1:nY, i = 1:nX
      Bmag2 = Bx[i,j,k]^2 + Bx[i,j,k]^2 + Bz[i,j,k]^2
      vDiaMag[1,i,j,k] = (By[i,j,k]*pz[i,j,k] - Bz[i,j,k]*py[i,j,k]) /
         (ni[i,j,k]*Bmag2) * cUnit
      vDiaMag[2,i,j,k] = (Bz[i,j,k]*px[i,j,k] - Bx[i,j,k]*pz[i,j,k]) /
         (ni[i,j,k]*Bmag2) * cUnit
      vDiaMag[3,i,j,k] = (Bx[i,j,k]*py[i,j,k] - By[i,j,k]*px[i,j,k]) /
         (ni[i,j,k]*Bmag2) * cUnit
      vMag[i,j,k] = hypot(vDiaMag[:,i,j,k]...)
   end

   yMid = floor(Int,nY/2)
   xrange = floor(Int,nX/4):floor(Int,nX/4*3)
   zrange = floor(Int,nZ/4):floor(Int,nZ/4*3)

   clf()

   fig, ax = subplots()
   cont_levels = range(0.0, 10000.0, length=50)
   c = ax.contourf(x[xrange,yMid,zrange], z[xrange,yMid,zrange],
      vMag[xrange,yMid,zrange], levels=cont_levels)
   fig.colorbar(c, ax=ax, ticks=range(0.0, 10000.0, step=1000.0))
   xlabel("x"); ylabel("z")
   title(L"$\mathbf{B}\times \nabla P/(neB^2)$"*", t=$(i)")

   fig.set_size_inches(3.5,6)
   plt.axis("scaled")

   plt.savefig(saveDir*"/diamagnetic_$(i).png")
end

#=
filename = "box*outs";
dir = "/Users/hyzhou/Ganymede/run_mercury_80s_newbox/GM";

filename = "y*outs";
dir = "/Users/hyzhou/Ganymede/run_mercury_80s_newbox/GM";
=#

#filename = "cut*.dat";
#dir = "GM/IO2";
#vtk_auto_conversion(filename, dir)

function y_by_cuts()
   np = pyimport("numpy")
   filename = "y*out"
   data = readdata(filename)

   fig, ax = plt.subplots(1, 1)
   var = "by"; plotinterval = 0.1
   VarIndex_ = findfirst(x->x==lowercase(var), lowercase.(data.head.wnames))
   X = vec(data.x[:,:,1])
   Y = vec(data.x[:,:,2])
   W = vec(data.w[:,:,VarIndex_])

   # Create grid values first.
   xi = range(minimum(X), stop=maximum(X), step=plotinterval)
   yi = range(minimum(Y), stop=maximum(Y), step=plotinterval)
   # Perform linear interpolation of the data (x,y) on grid(xi,yi)
   triang = matplotlib.tri.Triangulation(X,Y)
   interpolator = matplotlib.tri.LinearTriInterpolator(triang, W)
   Xi, Yi = np.meshgrid(xi, yi)
   wi = interpolator(Xi, Yi)

   c = contourf(xi,yi,wi,50)
   fig.set_size_inches(6,6)
   plt.axis("scaled")
   colorbar(c, ax=ax,ticks=np.arange(-100.0, 0.0, 12.5))
   xlabel("x"); ylabel("z")
   title(L"By")
   plt.savefig("test.png")
end

function plot_beta(filename::String)

   data = readdata(filename)

   cutPlaneIndex = 65
   VarIndex_ = 18

   X = @view data.x[:,:,:,1]
   Y = @view data.x[:,:,:,2]
   Z = @view data.x[:,:,:,3]

   fig, ax = plt.subplots(1, 1)
   fig.set_size_inches(3.3,6);

   W = @view data.w[:,:,:,VarIndex_]
   Bx = @view data.w[:,:,:,5]
   By = @view data.w[:,:,:,6]
   Bz = @view data.w[:,:,:,7]

   cut1 = @view X[:,cutPlaneIndex,:]
   cut2 = @view Z[:,cutPlaneIndex,:]
   W    = @view W[:,cutPlaneIndex,:]
   Bx   = @view Bx[:,cutPlaneIndex,:]
   By   = @view By[:,cutPlaneIndex,:]
   Bz   = @view Bz[:,cutPlaneIndex,:]
   PB   = hypot.(Bx, By, Bz)

   c = ax.contourf(cut1, cut2, W./PB)
   fig.colorbar(c,ax=ax)
   ax.axis("scaled")
   title(data.head.wnames[VarIndex_])

   xlabel("x"); ylabel("z")

   dim = [0.125, 0.013, 0.2, 0.045]
   str = @sprintf "it=%d, time=%4.2f" filehead[:it] filehead[:time]
   at = matplotlib.offsetbox.AnchoredText(str,
      loc="lower left", prop=Dict("size"=>8), frameon=true,
      bbox_to_anchor=(0., 1.),
      bbox_transform=ax.transAxes)
      at.patch.set_boxstyle("round,pad=0.,rounding_size=0.2")
      ax.add_artist(at)
   #plt.savefig("beta.png")
end
