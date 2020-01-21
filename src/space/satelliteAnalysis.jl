# Static satellite analysis
#
# Hongyang Zhou, hyzhou@umich.edu 12/28/2019

# In : -1.525      -2.269       0.694
# Out: -1.212       2.771       0.785

using DelimitedFiles, Statistics, PyPlot, LinearAlgebra, Distributions
using PyCall
gridspec = pyimport("matplotlib.gridspec")

const μ₀ = 4π*1e-7 # Vacuum permeability, [H/m]
const amu = 1.66054e-27 # [kg]
# Rho Ux Uy Uz Bx By Bz Pe P, NaN for compensating the No index in the 1st col.
const upstream_value = [NaN, 56., 140., 0., 0., -10., -6., -86., 0.2, 3.4]
const γ = 5/3 # adiabatic index

struct Index
   Rho_::Integer
   Ux_::Integer
   Uy_::Integer
   Uz_::Integer
   Bx_::Integer
   By_::Integer
   Bz_::Integer
   Pe_::Integer
   P_::Integer
   B_::UnitRange{Int64}
   U_::UnitRange{Int64}
end

"""Import satellite data. Require one line header presented."""
function read_data(fname)
   f = readdlm(fname, ',', Float32, '\n'; header=true)

   header = f[2][1:end-1]
   data   = f[1]

   satelliteNo = unique(data[:,1]) # number of static satellites

   return header, data, satelliteNo
end

function get_indexes(header::Vector{AbstractString})
   Rho_= findfirst(x->x=="Rho", header) + 1
   P_  = findfirst(x->x=="P", header) + 1
   Pe_ = findfirst(x->x=="Pe", header) + 1
   Bx_ = findfirst(x->x=="Bx", header) + 1
   By_ = Bx_ + 1
   Bz_ = Bx_ + 2
   B_  = Bx_:(Bx_+2)
   Ux_ = findfirst(x->x=="Ux", header) + 1
   Uy_ = Ux_ + 1
   Uz_ = Ux_ + 2
   U_  = Ux_:(Ux_+2)

   id = Index(Rho_, Ux_, Uy_, Uz_, Bx_, By_, Bz_, Pe_, P_, B_, U_)
end

"""
	single_satellite_plot(filename, dir, nShift)

Plot all MHD quantities for a static location during the simulation.
# Arguments
- `filename::String`: Input file generated by ParaView.
- `dir::String`: file directory.
- `nShift::Integer`: shift w.r.t. the satellite starting location.
"""
function single_satellite_plot(filename="satellites_PIC.txt",
   dir="/Users/hyzhou/Documents/Computer/ParaView/data/", nShift=91)

   @assert nShift≥0 "nShift must be non-negative!"

   header, data, satelliteNo = read_data(dir*filename)

   index_ = findall(x->x==satelliteNo[1], data[:,1])

   id = get_indexes(header)
   Rho_, Ux_, Uy_, Uz_, Bx_, By_, Bz_, Pe_, P_, U_, B_ =
   id.Rho_, id.Ux_, id.Uy_, id.Uz_, id.Bx_, id.By_, id.Bz_, id.Pe_, id.P_,
   id.U_, id.B_

   # magnetic pressure, [nPa]
   PB = sum(data[index_ .+ nShift,B_].^2; dims=2) ./ (2*μ₀) ./ (1e9)

   # dynamic pressure, [nPa]
   Pd = sum(data[index_ .+ nShift,Rho_] .* data[index_ .+ nShift,U_].^2; dims=2) * amu * 1e21


   # All quantities at one location
   figure(figsize=(12.5,6))
   plt.rc("font", family="serif", size=14)
   gs = gridspec.GridSpec(4, 1)
   gs.update(wspace=0.0, hspace=0.3) # set the spacing between axes.
   f1 = subplot(get(gs, (0, 0))) # rho

   plot(data[index_ .+ nShift,Rho_], label=L"\rho\ [amu/cc]")
   xlim(0,1200)
   plt.axhline(y=upstream_value[Rho_], color="k", linestyle="--", alpha=0.5)
   grid(true)
   f1.axes.xaxis.set_ticklabels([])
   legend(loc="lower left", bbox_to_anchor=(0.0, 0.8), ncol=4, frameon=false)

   f2 = subplot(get(gs, (1, 0))) # p
   plot(data[index_ .+ nShift,P_], alpha=0.9, label=L"P_i")
   plot(data[index_ .+ nShift,Pe_], alpha=0.9, label=L"P_e")
   plot(Pd, alpha=0.9, label=L"P_d")
   plot(PB, alpha=0.9, label=L"P_B\ [nPa]")
   xlim(0,1200)
   ylim(0,maximum(Pd[60:end])+0.02)
   legend(loc="lower left", bbox_to_anchor=(0.0, 0.8), ncol=4, frameon=false)
   plt.axhline(y=upstream_value[P_], color="k", linestyle="--", alpha=0.5)
   grid(true)
   f2.axes.xaxis.set_ticklabels([])

   f3 = subplot(get(gs, (2, 0))) # U
   plot(data[index_ .+ nShift,Ux_], alpha=0.9, label="Ux")
   plot(data[index_ .+ nShift,Uy_], alpha=0.9, label="Uy")
   plot(data[index_ .+ nShift,Uz_], alpha=0.7, label="Uz  [km/s]")
   xlim(0,1200)
   #ylim(-100,200)#ylim(-70,180)
   legend(loc="lower left", bbox_to_anchor=(0.0, 0.8), ncol=3, frameon=false)
   plt.axhline(y=upstream_value[Ux_], color="k", linestyle="--", alpha=0.5)
   plt.axhline(y=upstream_value[Uy_], color="k", linestyle="--", alpha=0.5)
   plt.axhline(y=upstream_value[Uz_], color="k", linestyle="--", alpha=0.5)
   grid(true)
   f3.axes.xaxis.set_ticklabels([])

   f4 = subplot(get(gs, (3, 0))) # B
   plot(data[index_ .+ nShift,Bx_], alpha=0.9, label="Bx")
   plot(data[index_ .+ nShift,By_], alpha=0.9, label="By")
   plot(data[index_ .+ nShift,Bz_], alpha=0.7, label="Bz  [nT]")
   xlim(0,1200)
   legend(loc="lower left", bbox_to_anchor=(0.0, 0.8), ncol=3, frameon=false)
   plt.axhline(y=upstream_value[Bx_], color="k", linestyle="--", alpha=0.5)
   plt.axhline(y=upstream_value[By_], color="k", linestyle="--", alpha=0.5)
   #plt.axhline(y=upstream_value[Bz_], color="k", linestyle="--", alpha=0.5)
   grid(true)
   xlabel("simulation time [s]")

   suptitle("Near magnetopause, location: "*string(data[nShift+1,11:end])[8:end]*L"R_G", fontsize=14)
   #plt.tight_layout(rect=[0, -0.02, 1, 0.98])
end

"""
	multi_satellite_plot(filename, dir)
# Arguments
- `filename::String`: Input file generated by ParaView.
- `dir::String`: file directory.
- `var::String`: variable for plotting.
"""
function multi_satellite_plot(filename="satellites_PIC.txt",
   dir="/Users/hyzhou/Documents/Computer/ParaView/data/", var="Rho")

   header, data, satelliteNo = read_data(dir*filename)

   var_ = findfirst(x->x==var, header) + 1

   index_ = findall(x->x==satelliteNo[1], data[:,1])

   shift = [i*46 for i in 0:5]
   #shift = [i + 88 for i in 0:5]
   satelliteLoc = data[1 .+ shift,11:13]

   # one variable at 6 locations
   figure(figsize=(10,6))
   for iplot in 1:6
      plt.subplot(3,2,iplot)
      plot(data[index_ .+ shift[iplot],var_])
      plt.axhline(y=upstream_value[var_], color="r", linestyle="--", alpha=0.5)
      plt.title(string(satelliteLoc[iplot,:])[8:end])
      grid(true)
      f1 = plt.gca()
      if iplot ∉ [5,6]
         f1.axes.xaxis.set_ticklabels([])
         #f1.axes.yaxis.set_ticklabels([])
      end
   end

   suptitle("$var at 6 satellite locations", fontsize=14)
   plt.tight_layout()

end

"""
	multi_satellite_contour(filename, dir, DoSave, DoSubtractMean)

Static satellite analysis, contour plots of location and time.
`x` is the satellite position, `y` is the time, and z are values.
# Arguments
- `filename::String`: Input file generated by ParaView.
- `dir::String`: file directory.
- `DoSave::Bool`: Save the plots.
- `DoSubtractMean::Bool`: Subtract the average state for each variable.
"""
function multi_satellite_contour(filename="satellites_PIC.txt",
   dir="/Users/hyzhou/Documents/Computer/ParaView/data/"; DoSave=false,
   DoSubtractMean = true)

   header, data, satelliteNo = read_data(dir*filename)

   # Remove the trailing satellites
   #satelliteNo = satelliteNo[1:end-20]

   index_ = findall(x->x==satelliteNo[1], data[:,1])

   c = Array{Float32, 2}(undef, length(index_), length(satelliteNo))

   # Subtract the average
   for var in header[1:end-3]
      var_ = findfirst(x->x==var, header) + 1
      @show var, var_

      for i in satelliteNo
         c[:,Int(i+1)] = data[index_ .+ Int(i),var_]
      end

      cmean = mean(c, dims=1)

      figure()
      if DoSubtractMean
         contourf(c .- cmean, 50)
      else
         contourf(c, 50)
      end

      plt.set_cmap("plasma")
      colorbar()
      title(var)

      if DoSave
         if occursin("pic",lowercase(filename))
            subname = "pic"
         elseif occursin("hall",lowercase(filename))
            subname = "hall"
         end
         if DoSubtractMean
            savefig(var*"_traj_contour_subtractMean_"*subname*".png")
         else
            savefig(var*"_traj_contour_"*subname*".png")
         end
         close()
      end
   end

end


"""
	smooth(x, n)

Return the moving box average of the vector data `x` with box length 'n'.
One-sided average on the left and right edge.
"""
function smooth(x::Vector, n::Int=100)
   nx = length(x)
   x̄ = zeros(eltype(x),length(x))

   # left points
   for i = 1:n
      x̄[i] = mean(x[1:(i+n)])
   end

   # middle points
   for i = (n+1):(nx-n)
      x̄[i] = mean(x[(i-n):(i+n)])
   end

   # right points
   for i = (nx-n+1):nx
      x̄[i] = mean(x[(i-n):end])
   end

   return x̄
end

"""
	wave_analysis(nShift, DoPlot, filename, dir)

Static satellite wave analysis.
`x` is the satellite position, `y` is the time, and z are values.
# Arguments
- `nShift::Integer`: shift w.r.t. the satellite starting location.
- `DoPlot::Bool`: Plot the perturbations.
- `filename::String`: Input file generated by ParaView.
- `dir::String`: file directory.
"""
function wave_analysis(nShift=115; DoPlot=false, filename="satellites_PIC.txt",
   dir="/Users/hyzhou/Documents/Computer/ParaView/data/", verbose=true)

   header, data, satelliteNo = read_data(dir*filename)

   index_ = findall(x->x==satelliteNo[1], data[:,1])

   verbose && println("satellite location: ",nShift, " ",data[nShift+1,11:end])

   id = get_indexes(header)
   Rho_, Ux_, Uy_, Uz_, Bx_, By_, Bz_, Pe_, P_, U_, B_ =
   id.Rho_, id.Ux_, id.Uy_, id.Uz_, id.Bx_, id.By_, id.Bz_, id.Pe_, id.P_,
   id.U_, id.B_

   ρ = data[index_ .+ nShift,Rho_]
   P = data[index_ .+ nShift,P_]
   Pe= data[index_ .+ nShift,Pe_]
   Ux= data[index_ .+ nShift,Ux_]
   Uy= data[index_ .+ nShift,Uy_]
   Uz= data[index_ .+ nShift,Uz_]
   Bx= data[index_ .+ nShift,Bx_]
   By= data[index_ .+ nShift,By_]
   Bz= data[index_ .+ nShift,Bz_]
   B = @. √(Bx^2 + By^2 + Bz^2)

   # Subtract the background
   ρ̄ = smooth(ρ)
   P̄ = smooth(P)
   P̄e= smooth(Pe)
   Ūx= smooth(Ux)
   Ūy= smooth(Uy)
   Ūz= smooth(Uz)
   B̄x= smooth(Bx)
   B̄y= smooth(By)
   B̄z= smooth(Bz)
   B̄ = smooth(B)

   Δρ = ρ .- ρ̄
   ΔUx= Ux .- Ūx
   ΔUy= Uy .- Ūy
   ΔUz= Uz .- Ūz
   ΔBx= Bx .- B̄x
   ΔBy= By .- B̄y
   ΔBz= Bz .- B̄z
   ΔB = @. sqrt(Bx^2 + By^2 + Bz^2) - sqrt(B̄x^2 + B̄y^2 + B̄z^2)
   ΔP = @. P - P̄ + (Pe - P̄e)

   VA = @. B/√(μ₀*ρ*amu)*1e-15 # [km/s]

   start_ = 100 # remove the starting points in correlation analysis

   # Determine if ΔB and ΔU are correlated --> Alfven wave?
   rUBx = cor((ΔUx./VA)[start_:end],(ΔBx./B)[start_:end])
   rUBy = cor((ΔUy./VA)[start_:end],(ΔBy./B)[start_:end])
   rUBz = cor((ΔUz./VA)[start_:end],(ΔBz./B)[start_:end])

   rAlfven = sort([rUBx, rUBy, rUBz])[2] # Take the intermediate value

   if verbose
      println("correlation between ΔUx and ΔBx = ", round(rUBx, digits=2))
      println("correlation between ΔUy and ΔBy = ", round(rUBy, digits=2))
      println("correlation between ΔUz and ΔBz = ", round(rUBz, digits=2))
   end

   n = length(ρ)
   α = 0.05 # acceptance quantile

   t = abs(rUBx) / √(1-rUBx^2)*√(n-start_-2)
   if t ≤ quantile(TDist(n-start_-2), 1-α)
      verbose &&
      println("ΔUx and ΔBx hypothesis test rejected! No strong evidence of correlation!")
   end

   t = abs(rUBy) / √(1-rUBy^2)*√(n-start_-2)

   if t ≤ quantile(TDist(n-start_-2), 1-α)
      verbose &&
      println("ΔUy and ΔBy hypothesis test rejected! No strong evidence of correlation!")
   end

   t = abs(rUBz) / √(1-rUBz^2)*√(n-start_-2)
   if t ≤ quantile(TDist(n-start_-2), 1-α)
      verbose &&
      println("ΔUz and ΔBz hypothesis test rejected! No strong evidence of correlation!")
   end

   # Determine if ΔB and Δρ are correlated --> fast/slow wave?
   rUρ = cor(ΔB[start_:end], Δρ[start_:end]) # this one is not actually used

   verbose && println("correlation between ΔB and Δρ = ", round(rUρ, digits=2))

   t = abs(rUρ) / √(1-rUρ^2)*√(n-start_-2)

   if t ≤ quantile(TDist(n-start_-2), 1-α)
      verbose &&
      println("ΔB and Δρ hypothesis test rejected! No strong evidence of correlation!")
   end

   # thermal pressure perturbation vs. magnetic pressure perturbation
   ΔPB = @. (Bx^2 + By^2 + Bz^2 - B̄x^2 - B̄y^2 - B̄z^2)/(2*μ₀)*1e-9
   rPBPt = cor(ΔPB[start_:end], ΔP[start_:end])

   verbose && println("correlation between ΔPB and ΔPt = ", round(rPBPt, digits=2))

   if rPBPt > 0.5
      verbose && println("possibly fast magnetosonic wave")
   elseif rPBPt < -0.5
      verbose && println("possibly slow magnetosonic wave")
   end

   if DoPlot
      figure(figsize=(12,3.5))
      plt.rc("font", family="serif", size=14)
      plot(ΔP, label=L"\Delta P_t")
      plot(ΔPB, label=L"\Delta P_B")
      xlim(100,1200)
      ylim(-3.5,3.5)
      legend(loc="best", ncol=2, frameon=false)
      xlabel("simulation time [s]")
      plt.axhline(y=0.0, color="k", linestyle="--", alpha=0.5)
      plt.annotate("cc=$(round(rPBPt,digits=2))", xy=(0.9, 0.05),
         xycoords="axes fraction")
      #plt.title("Hall MHD, Location: "*string(data[nShift+1,11:end])[8:end]*L"R_G")
      plt.title("MHD-EPIC, Location: "*string(data[nShift+1,11:end])[8:end]*L"R_G")
      tight_layout()
   end

   # Energy perturbations
   Ek = @. (((ΔUx*Ux+ΔUy*Uy+ΔUz*Uz)+0.5*(ΔUx^2 + ΔUy^2 + ΔUz^2))*ρ̄ +
      Δρ*(Ūx^2 + Ūy^2 + Ūz^2))*amu*1e12
   Et = @. ΔP*1e-9/(γ-1)
   Eb = @. ((ΔBx*Bx+ΔBy*By+ΔBz*Bz)+0.5*(ΔBx^2 + ΔBy^2 + ΔBz^2))/μ₀*1e-18

   verbose && println("correlation between Ek and Et = ",
      round(cor(Ek[start_:end], Et[start_:end]), digits=2))
   verbose && println("correlation between Ek and Eb = ",
      round(cor(Ek[start_:end], Eb[start_:end]), digits=2))
   verbose && println("correlation between Et and Eb = ",
      round(cor(Et[start_:end], Eb[start_:end]), digits=2))

   #println("Average kinetic energy = ",mean(Ek))
   #println("Average magnetic energy = ",mean(Eb))

   if DoPlot
      figure(figsize=(12,5.5))
      plt.rc("font", family="serif", size=14)
      gs = gridspec.GridSpec(3, 1)
      gs.update(wspace=0.0, hspace=0.0) # set the spacing between axes.
      f1 = subplot(get(gs, (0, 0)))
      plot(ΔUx ./ VA, label=L"\Delta U_x / V_A")
      plot(ΔBx ./ B, label=L"\Delta B_x / B_0")
      plt.axhline(y=0.0, color="k", linestyle="--", alpha=0.5)
      xlim(0,1200)
      legend(loc="best", ncol=2, frameon=false)
      f1.axes.xaxis.set_ticklabels([])
      plt.annotate("cc=$(round(rUBx,digits=2))", xy=(0.9, 0.05),
         xycoords="axes fraction")
      f2 = subplot(get(gs, (1, 0)))
      plot(ΔUy ./ VA, label=L"\Delta U_y / V_A")
      plot(ΔBy ./ B, label=L"\Delta B_y / B_0")
      plt.axhline(y=0.0, color="k", linestyle="--", alpha=0.5)
      xlim(0,1200)
      legend(loc="best", ncol=2, frameon=false)
      f2.axes.xaxis.set_ticklabels([])
      plt.annotate("cc=$(round(rUBy,digits=2))", xy=(0.9, 0.05),
         xycoords="axes fraction")
      subplot(get(gs, (2, 0)))
      plot(ΔUz ./ VA, label=L"\Delta U_z / V_A")
      plot(ΔBz ./ B, label=L"\Delta B_z / B_0")
      plt.axhline(y=0.0, color="k", linestyle="--", alpha=0.5)
      xlim(0,1200)
      legend(loc="best", ncol=2, frameon=false)
      xlabel("simulation time [s]")
      plt.annotate("cc=$(round(rUBz,digits=2))", xy=(0.9, 0.05),
         xycoords="axes fraction")
      #suptitle("Hall MHD, Location: "*string(data[nShift+1,11:end])[8:end]*L"R_G")
      suptitle("MHD-EPIC, Location: "*string(data[nShift+1,11:end])[8:end]*L"R_G")

      figure(figsize=(12,3.5))
      plt.rc("font", family="serif", size=14)
      plot(ΔB, alpha=0.95, label=L"\Delta B")
      plot(Δρ, alpha=0.95, label=L"\Delta \rho")
      xlim(0,1200)
      #ylim(-3,3)
      legend(loc="best", ncol=2, frameon=false)
      xlabel("simulation time [s]")
      plt.title("Location: "*string(data[nShift+1,11:end])[8:end]*L"R_G")
      tight_layout()

      #=
      figure()
      plot(Ek, alpha=0.7, label="Ek")
      plot(Et, alpha=0.7, label="Et")
      plot(Eb, alpha=0.7, label="Eb")
      legend()
      =#
   end

   #=
   figure()
   plot(ΔUx, alpha=0.7)
   plot(ΔP .* -50, alpha=0.7)
   xlabel("delta Ux")
   ylabel("delta P")
   =#

   #=
   # Check the linearity
   figure(figsize=(12,4))
   subplot(131)
   scatter(ΔUx,ΔBx)
   subplot(132)
   scatter(ΔUy,ΔBy)
   subplot(133)
   scatter(ΔUz,ΔBz)

   figure()
   scatter(ΔB,Δρ)
   =#
   return rAlfven, rPBPt
end

function check_wave_type(filename="satellites_PIC.txt",
   dir="/Users/hyzhou/Documents/Computer/ParaView/data/")
   f = readdlm("/Users/hyzhou/Documents/Computer/ParaView/data/G8Traj.csv",
      ',', Float32, '\n'; header=true)
   nS = size(f[1])[1]

   # X, Y, Z, rAlfven, rFastSlow
   fnew = Array{Float32,2}(undef, nS, 5)

   header, data, satelliteNo = read_data(dir*filename)

   if nS != length(satelliteNo)
      error("Number of satellites not right! Check the input files!")
   end

   for iShift = 1:nS
      r1,r2 = wave_analysis(iShift-1; DoPlot=false, filename=filename,
         dir=dir, verbose=false)
      fnew[iShift,1:3] = f[1][iShift,:]
      fnew[iShift,3] = 0.75 # For plotting purpose
      fnew[iShift,4:5] .= r1, r2
      println("iShift = ", iShift)
   end

   # Write to file
   open(dir*"waveAlongTrajG8.csv", "w") do io
      write(io, "\"X\",\"Y\",\"Z\",\"rA\",\"rM\"\n")
      writedlm(io, fnew, ',')
   end

   return fnew
end

#single_satellite_plot("satellites_PIC.txt",
#   "/Users/hyzhou/Documents/Computer/ParaView/data/", 185)
#multi_satellite_plot()
#multi_satellite_contour("satellites_boundary_PIC.txt", DoSubtractMean=false)


nShift = 185
#single_satellite_plot("satellites_Hall.txt",
#   "/Users/hyzhou/Documents/Computer/ParaView/data/", nShift)
wave_analysis(nShift; DoPlot=true, filename="satellites_PIC.txt",
   verbose=true)
#fnew = check_wave_type()
#check_wave_type()

#=
# Cross correlation check
using DSP

figure()
subplot(411)
plot(xcorr(ΔUx, ΔBx))
subplot(412)
plot(xcorr(ΔUy, ΔBy))
subplot(413)
plot(xcorr(ΔUz, ΔBz))
subplot(414)
plot(xcorr(ΔB,Δρ))
=#

#=
using StatsBase
figure()
subplot(411)
plot(crosscor(ΔUx, ΔBx))
subplot(412)
plot(crosscor(ΔUy, ΔBy))
subplot(413)
plot(crosscor(ΔUz, ΔBz))
subplot(414)
plot(crosscor(ΔB,Δρ))
=#
