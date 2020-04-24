# Cross Polar Cap Potential (CPCP) calculation at z=2.0 cut plane.
#
# Part 2: Doing surface flux integral and line integral to obtain CPCP.
#
# Hongyang Zhou, hyzhou@umich.edu 03/22/2020

using DelimitedFiles, Glob, PyPlot, Statistics, LinearAlgebra
using NumericalIntegration, ConcaveHull, Random

const Rg = 2634000.                  # [m], radius of Ganymede
const e  = 1.60217662e-19            # [C], electron charge
const Δx = 1/32                      # [Rg],grid resolution
const Δ = Δx*Rg                      # [m]
const Δ²= Δ^2                        # [m^2]
const CoefP = 1e3 * 1e-9 * Rg * 1e-3 # [km/s]*[nT]*[Rg] --> [kV]
const CoefJ = 1e-6 /e * 1e-6* 1e-3  # [μA/m^2]/[#/cc] --> [km/s]
const ϵ = Float32(1e-5) # for comparing floating point values

"Integrate total magnetic flux over the closed surface."
function integrate_surface(data; DoPlot=false, verbose=false)

   bz= @view data[:,7] # [nT]
   x = @view data[:,11]
   y = @view data[:,12]

   # I can also use sort instead.
   bzLow = quantile(bz,0.01)
   bzReasonableRange = quantile(bz,0.95) - quantile(bz,0.05)

   ϕ = 0.0f0
   for z in bz
      bzLow - z > bzReasonableRange && continue # skip suspicious values
      ϕ += z
   end

   if verbose
      b̄z = mean(bz)
      Random.seed!(1234)
      xC = [i + ϵ * rand(Float32) for i in x]
      yC = [i + ϵ * rand(Float32) for i in y]
      v = [[xC[i], yC[i]] for i in 1:length(xC)]
      hull = concave_hull(v, 100) # ConcaveHull
      hull_area = area(hull)

      xB = [p[1] for p in hull.vertices]
      yB = [p[2] for p in hull.vertices]

      if DoPlot
         figure(figsize=(4.5,7.5))
         plt.rc("font", family = "serif", size = 14)
         scatter(xC, yC, s=1, c=bz)
         colorbar()
         plot(xB, yB, "r")

         xlabel(L"x\, [R_G]", size = 14)
         ylabel(L"y\, [R_G]", size = 14)
         title(L"B_z\, [nT]", size = 14)
         xlim(-2.,2.)
         #axis("scaled")
         #tight_layout()
      end

      @info "mean bz = $b̄z"
      @info "area = $hull_area"
      @info "estimate: $(b̄z*hull_area)"
      @info "compute: $(ϕ/32/32)"
   end

   ϕ *= Δ²*1e-9 # [T*m^2]
end

"Not to be confused with the built-in extrema function. Faster."
function minmax(x::AbstractArray)
   a = b = first(x)
   @inbounds @simd for i in eachindex(x)
      if x[i] > b
         b = x[i]
      elseif x[i] < a
         a = x[i]
      end
   end
   a, b
end

"Integrate electric field over the x=xMax line and return the difference."
function integrate_line(data; verbose=false)

   xMax = maximum(@view data[:,11])
   index_ = [i for i in 1:size(data)[1] if isapprox(data[i,11], xMax, atol=ϵ)]

   nP = length(index_)

   ρ = data[index_,1]
   ux= data[index_,2]
   uy= data[index_,3]
   uz= data[index_,4]
   bx= data[index_,5]
   by= data[index_,6]
   bz= data[index_,7]
   jx= data[index_,8]
   jy= data[index_,9]
   jz= data[index_,10]
   y = data[index_,12]

   # Get E field: E = -uxB
   Ex = Vector{Float64}(undef,nP)
   Ey = Vector{Float64}(undef,nP)
   Ez = Vector{Float64}(undef,nP)

   # Calculate electron bulk velocity in Hall MHD
   # Transform the hall velocity into [km/s]
   uex = @. ux - jx/ρ * CoefJ
   uey = @. uy - jy/ρ * CoefJ
   uez = @. uz - jz/ρ * CoefJ

   for i = 1:nP
      Ex[i], Ey[i], Ez[i] =
      -cross([uex[i], uey[i], uez[i]], [bx[i], by[i], bz[i]])
   end

   ϕy = cumul_integrate(y, Ey)
   ϕymin, ϕymax = minmax(ϕy)

   potential_line = (ϕymax - ϕymin) * CoefP # [kV]

   verbose && @info potential_line, y[end]-y[1]

   return potential_line, y[end]-y[1]
end

"""
	get_E_integral(filenames; ftype="hall")

Obtain the upstream electric field integral from surface flux integral and
middle line integration.
"""
function get_E_integral(filenames; ftype="hall")

   nfile = length(filenames)
   ϕ = zeros(Float32, nfile)
   potential_M     = Vector{Float32}(undef, nfile)
   potential_Mhalf = Vector{Float32}(undef, nfile-1)
   yLength         = Vector{Float32}(undef, nfile)

   for i = 1:nfile
      @info i
      fieldline = readdlm(filenames[i], ',', Float32, header=true)
      data = fieldline[1]

      ϕ[i] = integrate_surface(data)
      potential_M[i], yLength[i] = integrate_line(data,verbose=true)
   end

   potential_T = diff(ϕ) .* 1e-3 # [kV]

   # Get the potential drop in the middle line at half times
   for i in 1:nfile-1
      potential_Mhalf[i] = 0.5*(potential_M[i] + potential_M[i+1])
   end

   potential_U = potential_T .+ potential_Mhalf # Upstream potential, [kV]

   if ftype == "hall"
      t = deleteat!(collect(0:1200), [355, 457, 505, 568] .+ 1)
   else # PIC
      t = deleteat!(collect(300:1500), [567, 591, 610, 680] .- 300 .+ 1)
   end

   t_half = zeros(Float32, length(t)-1)

   # Get half times
   for i in 1:length(t)-1
      t_half[i] = 0.5*(t[i] + t[i+1])
   end

   outname = joinpath(dirname(filenames[1]),"CPCP_test.txt")
   isfile(outname) && rm(outname) # Clear the previous existing file

   for i = 1:length(t_half)
      open(outname, "a") do io
         writedlm(io, zip(t_half[i], round(potential_U[i],digits=4),
         round(potential_T[i],digits=4), round(potential_Mhalf[i],digits=4)))
      end
   end

   #@info mean(yLength)
end


function plot_surface(filenames)

   for i = 1:length(filenames)
      fig, ax = subplots(1,1)
      fieldline = readdlm(filenames[i], ',', Float32, header=true)
      data = fieldline[1]
      bz= @view data[:,7] # [nT]
      x = @view data[:,11]
      y = @view data[:,12]
      s = ax.scatter(x, y, s=1, c=bz)
      ax.set_xlim(-2.5,2.0)
      ax.set_ylim(-3.1,3.1)
      axis("scaled")
      fig.colorbar(s)
      title("t=$i")
      savefig("figure/ux"*lpad(i,4,"0")*".png")
      close(fig)
   end

end

"Upstream boundary line integral from closed region point data."
function upstream_boundary_integral(data; DoPlot=false)

   x = @view data[:,11]
   y = @view data[:,12]

   xMax = maximum(x)
   yMin, yMax = minmax(y)

   # Get boundary points through concave hull
   Random.seed!(1234)
   xC = [i + ϵ * rand(Float32) for i in x]
   yC = [i + ϵ * rand(Float32) for i in y]
   v = [[xC[i], yC[i]] for i in 1:length(xC)]
   hull = concave_hull(v, 100) # ConcaveHull
   xH = [p[1] for p in hull.vertices]
   yH = [p[2] for p in hull.vertices]

   # Find boundary indexes
   boundary_index_ = Int32[]
   for iH = axes(xH,1)
      for i = axes(x,1)
         if isapprox(xH[iH], x[i], atol=2ϵ) && isapprox(yH[iH], y[i], atol=2ϵ)
            push!(boundary_index_, i)
            break
         end
      end
   end

   xB = data[boundary_index_,11]
   yB = data[boundary_index_,12]

   line_index_ = @. ((xB == xMax) & (yB !== yMin) & (yB !== yMax))

   up_index_ = .!line_index_
   bc_seq_ = sortperm(yB[up_index_])
   ρ  = data[boundary_index_,1][up_index_][bc_seq_]
   ux = data[boundary_index_,2][up_index_][bc_seq_]
   uy = data[boundary_index_,3][up_index_][bc_seq_]
   uz = data[boundary_index_,4][up_index_][bc_seq_]
   bx = data[boundary_index_,5][up_index_][bc_seq_]
   by = data[boundary_index_,6][up_index_][bc_seq_]
   bz = data[boundary_index_,7][up_index_][bc_seq_]
   jx = data[boundary_index_,8][up_index_][bc_seq_]
   jy = data[boundary_index_,9][up_index_][bc_seq_]
   jz = data[boundary_index_,10][up_index_][bc_seq_]
   xUp= xB[up_index_][bc_seq_]
   yUp= yB[up_index_][bc_seq_]

   if DoPlot
      fig, ax = plt.subplots(1,1,figsize=(6.0,7.5))
      plt.rc("font", family = "serif", size = 14)
      scatter(x, y, s=1, c=data[:,7])
      colorbar()
      plot(xB, yB, "C0", linewidth=3)
      plot([xB[end],xB[1]], [yB[end],yB[1]], "C0", linewidth=3)
      plot(xUp, yUp, "r", linewidth=3)
      #scatter(xB, yB, s=5, c="C0")
      #scatter(xUp, yUp, s=5, c="r")
      plot(xUp[1], yUp[1], "k^",markersize=8)
      plot(xUp[end], yUp[end],"kv",markersize=8)
      ax.minorticks_on()
      xlabel(L"x\, [R_G]", size = 14)
      ylabel(L"y\, [R_G]", size = 14)
      title(L"B_z\, [nT]", size = 14)
      xlim(-2.,2.)
      ax.annotate("U", xy=(-1, 2.03), xytext=(-1.5, 2.5),
         arrowprops=Dict("facecolor"=>"black", "shrink"=>0.02))
      ax.annotate("M", xy=(1.3, 0.1), xytext=(1.7, 0.),
         arrowprops=Dict("facecolor"=>"black", "shrink"=>0.02))
      ax.annotate("A", xy=(1.0, 0.0), xytext=(0., 0.), color="w",size=20)
      #=
      figure(figsize=(2.7,4.8))
      scatter(xB, yB, s=5, c="C0")
      scatter(xUp, yUp, s=5, c="r")
      axis("scaled")
      xlabel(L"x\, [R_G]")
      ylabel(L"y\, [R_G]")
      #title("Boundary points")
      tight_layout()
      =#
   end

   nP = length(xUp)

   # Get E field: E = -uxB
   Ex = Vector{Float64}(undef,nP)
   Ey = Vector{Float64}(undef,nP)
   Ez = Vector{Float64}(undef,nP)

   # Calculate electron bulk velocity in Hall MHD
   # Transform the hall velocity into [km/s]
   uex = @. ux - jx/ρ * CoefJ
   uey = @. uy - jy/ρ * CoefJ
   uez = @. uz - jz/ρ * CoefJ

   for i = 1:nP
      Ex[i], Ey[i], Ez[i] =
      -cross([uex[i], uey[i], uez[i]], [bx[i], by[i], bz[i]])
   end

   ϕup = cumul_integrate(xUp, Ex) .+ cumul_integrate(yUp, Ey)

   line_index_ = (xB .== xMax)
   line_seq_ = sortperm(yB[line_index_],rev=true)
   yL = yB[line_index_][line_seq_]
   xL = xB[line_index_][line_seq_]

   ρL  = data[boundary_index_,1][line_index_][line_seq_]
   uxL = data[boundary_index_,2][line_index_][line_seq_]
   uyL = data[boundary_index_,3][line_index_][line_seq_]
   uzL = data[boundary_index_,4][line_index_][line_seq_]
   bxL = data[boundary_index_,5][line_index_][line_seq_]
   byL = data[boundary_index_,6][line_index_][line_seq_]
   bzL = data[boundary_index_,7][line_index_][line_seq_]
   jxL = data[boundary_index_,8][line_index_][line_seq_]
   jyL = data[boundary_index_,9][line_index_][line_seq_]
   jzL = data[boundary_index_,10][line_index_][line_seq_]

   nL = length(xL)

   # Get E field: E = -uxB
   ExL = Vector{Float64}(undef,nL)
   EyL = Vector{Float64}(undef,nL)
   EzL = Vector{Float64}(undef,nL)

   # Calculate electron bulk velocity in Hall MHD
   # Transform the hall velocity into [km/s]
   uexL = @. uxL - jxL/ρL * CoefJ
   ueyL = @. uyL - jyL/ρL * CoefJ
   uezL = @. uzL - jzL/ρL * CoefJ

   for i = 1:nL
      ExL[i], EyL[i], EzL[i] =
      -cross([uexL[i], ueyL[i], uezL[i]], [bxL[i], byL[i], bzL[i]])
   end

   ϕL = cumul_integrate(yL, EyL)

   if DoPlot
      EIntUp = ϕup.*CoefP
      EIntL = ϕL.*CoefP .+ EIntUp[end]
      figure(figsize=(4.5,3.0))
      scatter(1,EIntUp[1],s=30,c="k",marker="^")
      scatter(nP,EIntUp[end],s=30,c="k",marker="v")
      plot(1:nP,EIntUp,"r")
      plot(nP:nP+1,[EIntUp[end],EIntL[1]],"r")
      plot(nP+1:nP+nL,EIntL, "C0")
      xlabel("Points along curve clockwise")
      ylabel("Electric field integral [kV]")
      tight_layout()
   end

   ϕmin, ϕmax = minmax(ϕup)

   potential_line = (ϕmax - ϕmin) * CoefP # [kV]

   if DoPlot
      @info "Boundary curve potential drop = ", potential_line
      @info "Middle line potential drop = ", (maximum(ϕL) - minimum(ϕL))*CoefP
   end
   return potential_line
end

## Potential drop from surface integration and line integration

#filenames = glob("surface_value*txt","dataSurf_Hall")
#get_E_integral(filenames, ftype="hall")

filenames = glob("surface_value*txt","dataSurf_PIC")
#get_E_integral(filenames, ftype="pic")

## Demonstration for one snapshot

i = 281#952
fieldline = readdlm(filenames[i], ',', Float32, header=true)
data = fieldline[1]
#phi1 = integrate_surface(data,DoPlot=true,verbose=true)
#potential_line1, yL1 = integrate_line(data,verbose=true)
upstream_boundary_integral(data,DoPlot=true)


## Upstream boundary curve integration from closed points data
#=
ftype = "pic"
outname = joinpath("dataBoundary_PIC","CPCP_test.txt")
isfile(outname) && rm(outname) # Clear the previous existing file

if ftype == "hall"
   t = deleteat!(collect(0:1200), [355, 457, 505, 568] .+ 1)
else # PIC
   t = deleteat!(collect(300:1500), [567, 591, 610, 680] .- 300 .+ 1)
end

potential_line = zeros(length(filenames))

for i in 1:length(filenames)
   @info i
   filename = filenames[i]
   fieldline = readdlm(filename, ',', Float32, header=true)
   data = fieldline[1]

   potential_line[i] = upstream_boundary_integral(data,DoPlot=false)

   open(outname, "a") do io
      writedlm(io, [t[i] round(potential_line[i],digits=4)])
   end
end
=#

##
# Next step is time-series plotting.
