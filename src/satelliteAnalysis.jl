# Static satellite analysis
##
# 1. static_location_plot
# 2. multi_satellite_plot
# 3. multi_satellite_contour
# 4. wave analysis
# 5. check_wave_type
##
# Galileo G8 In/Out locations in GPhiO cooridinates:
# In : -1.525      -2.269       0.694
# Out: -1.212       2.771       0.785
#
# Eventually this should become a module!
#
# Hongyang Zhou, hyzhou@umich.edu 12/28/2019

using DelimitedFiles, Statistics, PyPlot, LinearAlgebra, Distributions, PyCall
gridspec = pyimport("matplotlib.gridspec")

include("utility.jl")
if !isdefined(Main, :Ganymede)
   include("Ganymede.jl")
   using .Ganymede
end


"""
    static_location_plot(filename, dir, nShift)

Plot all MHD quantities for a static location during the simulation.
# Arguments
- `filename::String`: Input file generated by ParaView.
- `dir::String`: file directory.
- `nShift::Integer`: shift w.r.t. the satellite starting location.
"""
function static_location_plot(filename="satellites_PIC.txt",
   dir="/Users/hyzhou/Documents/Computer/ParaView/data/", nShift=91)

   @assert nShift≥0 "nShift must be non-negative!"

   header, data, satelliteNo = read_satellite_data(dir*filename)

   index_ = findall(x->x==satelliteNo[1], data[:,1])

   id = getIndex(header)
   Rho_, Ux_, Uy_, Uz_, Bx_, By_, Bz_, Pe_, P_, U_, B_ =
   id.Rho_+1, id.Ux_+1, id.Uy_+1, id.Uz_+1, id.Bx_+1, id.By_+1, id.Bz_+1,
   id.Pe_+1, id.P_+1, id.U_.+1, id.B_.+1

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
   plt.axhline(y=upstream_value[Rho_-1], color="k", linestyle="--", alpha=0.5)
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
   plt.axhline(y=upstream_value[P_-1], color="k", linestyle="--", alpha=0.5)
   grid(true)
   f2.axes.xaxis.set_ticklabels([])

   f3 = subplot(get(gs, (2, 0))) # U
   plot(data[index_ .+ nShift,Ux_], alpha=0.9, label="Ux")
   plot(data[index_ .+ nShift,Uy_], alpha=0.9, label="Uy")
   plot(data[index_ .+ nShift,Uz_], alpha=0.7, label="Uz  [km/s]")
   xlim(0,1200)
   #ylim(-100,200)#ylim(-70,180)
   legend(loc="lower left", bbox_to_anchor=(0.0, 0.8), ncol=3, frameon=false)
   plt.axhline(y=upstream_value[Ux_-1], color="k", linestyle="--", alpha=0.5)
   plt.axhline(y=upstream_value[Uy_-1], color="k", linestyle="--", alpha=0.5)
   plt.axhline(y=upstream_value[Uz_-1], color="k", linestyle="--", alpha=0.5)
   grid(true)
   f3.axes.xaxis.set_ticklabels([])

   f4 = subplot(get(gs, (3, 0))) # B
   plot(data[index_ .+ nShift,Bx_], alpha=0.9, label="Bx")
   plot(data[index_ .+ nShift,By_], alpha=0.9, label="By")
   plot(data[index_ .+ nShift,Bz_], alpha=0.7, label="Bz  [nT]")
   xlim(0,1200)
   legend(loc="lower left", bbox_to_anchor=(0.0, 0.8), ncol=3, frameon=false)
   plt.axhline(y=upstream_value[Bx_-1], color="k", linestyle="--", alpha=0.5)
   plt.axhline(y=upstream_value[By_-1], color="k", linestyle="--", alpha=0.5)
   #plt.axhline(y=upstream_value[Bz_], color="k", linestyle="--", alpha=0.5)
   grid(true)
   xlabel("simulation time [s]")

   suptitle("Near magnetopause, location: "*
      string(data[nShift+1,11:end])[8:end]*L"R_G", fontsize=14)
   #plt.tight_layout(rect=[0, -0.02, 1, 0.98])
end

"""
    multi_satellite_plot(filename, dir)

One variable `var` plotted at 6 locations.
"""
function multi_satellite_plot(filename="satellites_PIC.txt",
   dir="/Users/hyzhou/Documents/Computer/ParaView/data/", var="Rho")

   header, data, satelliteNo = read_satellite_data(dir*filename)

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
      plt.axhline(y=upstream_value[var_-1], color="r",linestyle="--", alpha=0.5)
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
    multi_satellite_contour(filename, dir; DoSave=false,
       DoSubtractMean = true, nLead=10, nTrail=10)

Static satellite analysis, contour plots of all the variables as a function of
location and time.
`x` is the satellite position, `y` is the time, and `z` are values.
# Arguments
- `filename::String`: Input file generated by ParaView.
- `dir::String`: file directory.
- `plane::Char`: {'y','z'} cut plane.
- `DoSave::Bool`: Save the plots.
- `DoSubtractMean::Bool`: Subtract the average state for each variable.
- `nLead::Int`: the number of points removing from the start.
- `nTrail::Int`: the number of points removing from the end.

Issue: nLead must be ≥ 1, which should not be the case!
"""
function multi_satellite_contour(filename="satellites_PIC.txt",
   dir="/Users/hyzhou/Documents/Computer/ParaView/data/"; plane='y',
   DoSave=false, DoSubtractMean = true, nLead=10, nTrail=10)

   header, data, satelliteNo = read_satellite_data(dir*filename)

   # Remove the leading and trailing satellites
   satelliteNo = satelliteNo[1+nLead:end-nTrail]

   index_ = findall(x->x==0.0f0, data[:,1])

   if length(index_) ≤ 2
      @error "number of snapshots must be equal or larger than 2!"
   end

   c = Array{Float32, 2}(undef, length(index_), length(satelliteNo))

   crange = (-9.5,9.5)
   clength= 40
   vdisp = range(crange..., length=11)
   vplot = range(crange..., length=clength)

   # Subtract the average
   for var in header[1:end-3]
      var_ = findfirst(x->x==var, header) + 1
      @show var, var_

      for i in eachindex(satelliteNo)
         c[:,i] = data[index_ .+ Int(satelliteNo[i]),var_]
      end

      figure(figsize=(4,8))
      if DoSubtractMean
         cmean = mean(c, dims=1)
         if plane == 'y'
            contourf(data[Int.(satelliteNo),end], 1:length(index_), c .- cmean, vplot)
         elseif plane == 'z'
            contourf(data[Int.(satelliteNo),end-1], 1:length(index_), c .- cmean, vplot, extend="both")
         end
         #contourf(data[Int.(satelliteNo),end], 1:length(index_), c .- cmean, 50)
         #contourf(c .- cmean)
      else
         contourf(data[Int.(satelliteNo),end], 1:length(index_),c,50)
      end

      if plane == 'y'
         xlabel(L"z\ [R_G]")
      elseif plane == 'z'
         xlabel(L"y\ [R_G]")
      end
      ylabel("simulation time [s]")
      plt.set_cmap("seismic")
      #colorbar()
      colorbar(boundaries=vplot, ticks=vdisp)
      title(var)
      tight_layout()

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


function satellite_p_contour(filename="satellites_y0_PIC.txt",
   dir="/Users/hyzhou/Documents/Computer/ParaView/data/"; plane='y',
   DoSubtractMean = true, nLead=10, nTrail=10, No=1)

   header, data, satelliteNo = read_satellite_data(dir*filename)

   # Remove the leading and trailing satellites
   satelliteNo = satelliteNo[1+nLead:end-nTrail]

   index_ = findall(x->x==0.0f0, data[:,1])

   p = Array{Float32, 2}(undef, length(index_), length(satelliteNo))
   pmean = similar(p)

   crange = (-9.5,9.5)
   clength= 40
   vdisp = range(crange..., length=11)
   vplot = range(crange..., length=clength)

   pe_ = 9  # index for electron pressure
   pi_ = 10 # index for ion pressure

   for i in eachindex(satelliteNo)
      p[:,i] = data[index_ .+ Int(satelliteNo[i]),pe_] .+
               data[index_ .+ Int(satelliteNo[i]),pi_]
      pmean[:,i] = sma(p[:,i],100) # box average over +-100s for each location
   end

   t = 1:length(index_)

   if plane == 'y'
      z = data[Int.(satelliteNo),end]
   elseif plane == 'z'
      y = data[Int.(satelliteNo),end-1]
   end

   DN = matplotlib.colors.DivergingNorm

   fig, ax = subplots(figsize=(3.5,7))
   plt.rc("font", family="serif", size=14)

   if DoSubtractMean
      if plane == 'y'
         z = data[Int.(satelliteNo),end]
         contourf(z, t, p .- pmean, norm=DN(0), vplot)
      elseif plane == 'z'
         y = data[Int.(satelliteNo),end-1]
         contourf(y, 1:length(index_), p .- pmean, norm=DN(0), vplot,
         extend="both")
      end

   else
      contourf(data[Int.(satelliteNo),end], t,c,50)
   end

   if plane == 'y'
      #zpickUp_ = [index for (index, value) in enumerate(z) if value > 0.2]
      #zUpMean = mean((p.-pmean)[:,zpickUp_], dims=2)

      #zpickDn_ = [index for (index, value) in enumerate(z) if value < -0.2]
      #zDnMean = mean((p.-pmean)[:,zpickDn_], dims=2)

      peakUp_index = Int[]
      peakDn_index = Int[]
      σUp = std(p[:,59]) # at z=0.5
      σDn = std(p[:,8]) # at z=-0.5
      @info σUp, σDn
      tGap = 9 # peaks must be differed by a time range to be picked
      pmean = mean(p, dims=1) # Averaged pressure at each location over time

      for i = 1:size(p,1)
         if p[i,59] - pmean[59] > 1.0σUp
         #if zUpMean[i] > 1.5σUp
            if isempty(peakUp_index)
               append!(peakUp_index, i)
            elseif i - peakUp_index[end] > tGap
               append!(peakUp_index, i)
            end
         end

         if p[i,8] - pmean[8] > 1.0σDn
         #if zDnMean[i] > 1.5σDn
            if isempty(peakDn_index)
               append!(peakDn_index, i)
            elseif i - peakDn_index[end] > tGap
               append!(peakDn_index, i)
            end
         end
      end

      @info "number of FTE = $(length(peakUp_index)+length(peakDn_index))"

      zUp = fill(data[Int.(satelliteNo),end][59], size(peakUp_index))
      zDn = fill(data[Int.(satelliteNo),end][8], size(peakDn_index))
      plot(zUp, peakUp_index, linestyle="", marker="P", markersize=6,
      markerfacecolor="gray", markeredgecolor="None")

      plot(zDn, peakDn_index, linestyle="", marker="P", markersize=6,
      markerfacecolor="k", markeredgecolor="None")

      xlabel(L"z\ [R_G]")
   elseif plane == 'z'
      xlabel(L"y\ [R_G]")
   end

   plt.set_cmap("seismic")
   # create an axes on the right side of ax. The width of cax will be 5%
   # of ax and the padding between cax and ax will be fixed at 0.05 inch.
   axes_grid1 = pyimport("mpl_toolkits.axes_grid1")
   divider = axes_grid1.make_axes_locatable(ax)
   cax = divider.append_axes("right", size="5%", pad=0.05)

   plt.colorbar(cax=cax, boundaries=vplot, ticks=vdisp)
   ax.set_title(L"\Delta P_t\ [nPa]")

   ax.annotate("($('a'+No-1))", xy=(-0.17, 1.0), xycoords="axes fraction")

   if occursin("pic",lowercase(filename))
      ax.set_ylabel("MHD-EPIC, simulation time [s]")
   else
      ax.set_ylabel("Hall MHD, simulation time [s]")
   end

   tight_layout()

   if plane == 'y'
      return peakUp_index, peakDn_index
   else
      return nothing
   end
end

"""
    satellite_p_contour_test(filename, dir; plane='y', nLead=10, nTrail=10,
       No=1, DoSubtractMean=true)

Test on different methods of obtaining the FTEs.
"""
function satellite_p_contour_test(filename="satellites_y0_PIC.txt",
   dir="/Users/hyzhou/Documents/Computer/ParaView/data/"; plane='y',
   DoSubtractMean=true, nLead=10, nTrail=10, No=1)

   header, data, satelliteNo = read_satellite_data(dir*filename)

   # Remove the leading and trailing satellites
   satelliteNo = satelliteNo[1+nLead:end-nTrail]
   # Starting index for the first satellite
   index_ = findall(x->x==0.0f0, data[:,1])
   # Pressure at each time for each satellite location
   p = Array{Float32, 2}(undef, length(index_), length(satelliteNo))
   #pmean = similar(p)

   crange = (-9.5,9.5)
   clength= 40
   vdisp = range(crange..., length=11)
   vplot = range(crange..., length=clength)

   # Subtract the average
   var_ = 10

   for i in eachindex(satelliteNo)
      p[:,i] = data[index_ .+ Int(satelliteNo[i]),var_]
      #pmean[:,i] = sma(p[:,i],100) # box average over +-100s for each location
   end
   pmean = mean(p, dims=1) # Averaged pressure at each location over time

   z = data[Int.(satelliteNo),end] # z coordinates of satellites

   zpickUp_ = [index for (index, value) in enumerate(z) if value > 0.2]
   zUpMean = mean((p.-pmean)[:,zpickUp_], dims=2)

   zpickDn_ = [index for (index, value) in enumerate(z) if value < -0.2]
   zDnMean = mean((p.-pmean)[:,zpickDn_], dims=2)

   peakUp_index = Int[]
   peakDn_index = Int[]
   σUp = std(p[:,59]) # may be changed
   σDn = std(p[:,8])  # may be changed
   tGap = 10 # peaks must be differed by a time range to be picked
   for i in eachindex(index_)
      #if c[i,59] - cmean[59] > 1.5σUp
      if zUpMean[i] > 1.5σUp
         if isempty(peakUp_index)
            append!(peakUp_index, i)
         elseif i - peakUp_index[end] > tGap
            append!(peakUp_index, i)
         end
      end

      #if c[i,8] - cmean[8] > 1.5σDn
      if zDnMean[i] > 1.5σDn
         if isempty(peakDn_index)
            append!(peakDn_index, i)
         elseif i - peakDn_index[end] > tGap
            append!(peakDn_index, i)
         end
      end
   end

   zMean = 0.5.*(zUpMean.+zDnMean)
   figure()
   plot(zMean)
   tight_layout()
   #[plt.axvline(x=i, color="black", linestyle="--") for i in peakUp_index]

   #figure()
   #plot(zDnMean)
   #[plt.axvline(x=i, color="black", linestyle="--") for i in peakDn_index]

   t = 1:length(index_)
   Δp = p[:,59].-pmean[59]

   fig, ax = subplots(figsize=(8,3))
   plot(t,Δp)
   fill_between(t, 0, 1, where=Δp .> 1.2σUp,
      color="green", alpha=0.5, transform=ax.get_xaxis_transform())
   title("ΔP at z=0.5")
   tight_layout()

   #return zMean
   return Δp
end


"""
    wave_plot(nShift; DoPlot=false, filename, dir, iPlot=1, verbose)

Static satellite wave analysis.
`x` is the satellite position, `y` is the time, and `z` are values.
# Arguments
- `nShift::Integer`: shift w.r.t. the satellite starting location.
- `DoPlot::Bool`: Plot the perturbations.
- `filename::String`: Input file generated by ParaView.
- `dir::String`: file directory.
- `iPlot::Integer=1`: label tag index in alphabetical order starting from (a).
- `verbose::Bool=false`: Display information.
"""
function wave_plot(nShift=115; DoPlot=false, filename="satellites_PIC.txt",
   dir="/Users/hyzhou/Documents/Computer/ParaView/data/", iPlot=1, verbose=true)

   header, data, satelliteNo = read_satellite_data(dir*filename)

   index_ = findall(x->x==satelliteNo[1], data[:,1])

   verbose && @info "satellite location: ",nShift, " ",data[nShift+1,11:end]

   id = getIndex(header)
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
   ρ̄ = sma(ρ)
   P̄ = sma(P)
   P̄e= sma(Pe)
   Ūx= sma(Ux)
   Ūy= sma(Uy)
   Ūz= sma(Uz)
   B̄x= sma(Bx)
   B̄y= sma(By)
   B̄z= sma(Bz)
   B̄ = sma(B)

   Δρ = ρ .- ρ̄
   ΔUx= Ux .- Ūx
   ΔUy= Uy .- Ūy
   ΔUz= Uz .- Ūz
   ΔBx= Bx .- B̄x
   ΔBy= By .- B̄y
   ΔBz= Bz .- B̄z
   ΔB = @. hypot(Bx, By, Bz) - hypot(B̄x, B̄y, B̄z)
   ΔP = @. P - P̄ + (Pe - P̄e)

   VA = @. B/√(μ₀*ρ*amu)*1e-15 # [km/s]

   start_ = 100 # remove the starting points in correlation analysis

   # Determine if ΔB and ΔU are correlated --> Alfven wave?
   rUBx = cor((ΔUx./VA)[start_:end],(ΔBx./B)[start_:end])
   rUBy = cor((ΔUy./VA)[start_:end],(ΔBy./B)[start_:end])
   rUBz = cor((ΔUz./VA)[start_:end],(ΔBz./B)[start_:end])

   rAlfven = sort([rUBx, rUBy, rUBz])[2] # Take the intermediate value

   if verbose
      @info "correlation between ΔUx and ΔBx = ", round(rUBx, digits=2)
      @info "correlation between ΔUy and ΔBy = ", round(rUBy, digits=2)
      @info "correlation between ΔUz and ΔBz = ", round(rUBz, digits=2)
   end

   n = length(ρ)
   α = 0.05 # acceptance quantile

   t = abs(rUBx) / √(1-rUBx^2)*√(n-start_-2)
   if t ≤ quantile(TDist(n-start_-2), 1-α)
      verbose && @info "ΔUx and ΔBx hypothesis test rejected!"*
         " No strong evidence of correlation!"
   end

   t = abs(rUBy) / √(1-rUBy^2)*√(n-start_-2)

   if t ≤ quantile(TDist(n-start_-2), 1-α)
      verbose && @info "ΔUy and ΔBy hypothesis test rejected!"*
         " No strong evidence of correlation!"
   end

   t = abs(rUBz) / √(1-rUBz^2)*√(n-start_-2)
   if t ≤ quantile(TDist(n-start_-2), 1-α)
      verbose && @info "ΔUz and ΔBz hypothesis test rejected!"*
         " No strong evidence of correlation!"
   end

   # Determine if ΔB and Δρ are correlated --> fast/slow wave?
   rUρ = cor(ΔB[start_:end], Δρ[start_:end]) # this one is not actually used

   verbose && @info "correlation between ΔB and Δρ = ", round(rUρ, digits=2)

   t = abs(rUρ) / √(1-rUρ^2)*√(n-start_-2)

   if t ≤ quantile(TDist(n-start_-2), 1-α)
      verbose && @info "ΔB and Δρ hypothesis test rejected!"*
         " No strong evidence of correlation!"
   end

   # thermal pressure perturbation vs. magnetic pressure perturbation
   ΔPB = @. (Bx^2 + By^2 + Bz^2 - B̄x^2 - B̄y^2 - B̄z^2)/(2*μ₀)*1e-9
   rPBPt = cor(ΔPB[start_:end], ΔP[start_:end])

   verbose && @info "correlation between ΔPB and ΔPt = ", round(rPBPt, digits=2)

   if rPBPt > 0.5
      verbose && @info "possibly fast magnetosonic wave"
   elseif rPBPt < -0.5
      verbose && @info "possibly slow magnetosonic wave"
   end

   if DoPlot
      figure(figsize=(10.0,2.0))
      ax = subplot(111)
      plt.rc("font", family="serif", size=12)
      plot(ΔP, label=L"\Delta P_t")
      plot(ΔPB, label=L"\Delta P_B")
      xlim(100,1200)
      #ylim(-0.16,0.16)
      ylim(-3.5,3.5)
      legend(loc=(0.02,-0.05), ncol=2, frameon=false)
      xlabel("simulation time [s]")
      plt.axhline(y=0.0, color="k", linestyle="--", alpha=0.5)
      plt.annotate("cc=$(round(rPBPt,digits=2))", xy=(0.9, 0.05),
         xycoords="axes fraction")
      if occursin("PIC", filename)
         titletext = "MHD-EPIC, Location: "*
            string(data[nShift+1,11:end])[8:end]*L"R_G"
      else
         titletext = "Hall MHD, Location: "*
            string(data[nShift+1,11:end])[8:end]*L"R_G"
      end
      plt.title(titletext)
      ax.annotate("($('a'+iPlot-1))", xy=(0.0, 1.1), xycoords="axes fraction")
      ax.minorticks_on()
      ax.tick_params(which="both",top=true, right=true)
      ax.tick_params(which="both", direction="in")
      tight_layout()
   end

   # Energy perturbations
   Ek = @. (((ΔUx*Ux+ΔUy*Uy+ΔUz*Uz)+0.5*(ΔUx^2 + ΔUy^2 + ΔUz^2))*ρ̄ +
      Δρ*(Ūx^2 + Ūy^2 + Ūz^2))*amu*1e12
   Et = @. ΔP*1e-9/(γ-1)
   Eb = @. ((ΔBx*Bx+ΔBy*By+ΔBz*Bz)+0.5*(ΔBx^2 + ΔBy^2 + ΔBz^2))/μ₀*1e-18

   verbose && @info "correlation between Ek and Et = ",
      round(cor(Ek[start_:end], Et[start_:end]), digits=2)
   verbose && @info "correlation between Ek and Eb = ",
      round(cor(Ek[start_:end], Eb[start_:end]), digits=2)
   verbose && @info "correlation between Et and Eb = ",
      round(cor(Et[start_:end], Eb[start_:end]), digits=2)

   #println("Average kinetic energy = ",mean(Ek))
   #println("Average magnetic energy = ",mean(Eb))

   if DoPlot
      figure(figsize=(12.0,3.5))
      plt.rc("font", family="serif", size=12)
      gs = gridspec.GridSpec(3, 1)
      gs.update(wspace=0.0, hspace=0.0) # set the spacing between axes.
      f1 = subplot(get(gs, (0, 0)))
      plot(ΔUx ./ VA, label=L"\Delta U_x / V_A")
      plot(ΔBx ./ B, label=L"\Delta B_x / B_0")
      plt.axhline(y=0.0, color="k", linestyle="--", alpha=0.5)
      xlim(100,1200)
      #legend(loc=(0.36,-0.06), ncol=2, frameon=false)
      legend(loc=(0.00,0.6), ncol=2, frameon=false)
      f1.axes.xaxis.set_ticklabels([])
      f1.minorticks_on()
      f1.tick_params(which="both",top=true, right=true)
      f1.tick_params(which="both", direction="in")
      plt.annotate("cc=$(round(rUBx,digits=2))", xy=(0.885, 0.05),
         xycoords="axes fraction")

      f2 = subplot(get(gs, (1, 0)))
      plot(ΔUy ./ VA, label=L"\Delta U_y / V_A")
      plot(ΔBy ./ B, label=L"\Delta B_y / B_0")
      plt.axhline(y=0.0, color="k", linestyle="--", alpha=0.5)
      xlim(100,1200)
      #legend(loc=(0.36,-0.06), ncol=2, frameon=false)
      legend(loc=(0.00,0.6), ncol=2, frameon=false)
      f2.axes.xaxis.set_ticklabels([])
      f2.minorticks_on()
      f2.tick_params(which="both",top=true, right=true)
      f2.tick_params(which="both", direction="in")
      plt.annotate("cc=$(round(rUBy,digits=2))", xy=(0.885, 0.05),
         xycoords="axes fraction")

      f3 = subplot(get(gs, (2, 0)))
      plot(ΔUz ./ VA, label=L"\Delta U_z / V_A")
      plot(ΔBz ./ B, label=L"\Delta B_z / B_0")
      plt.axhline(y=0.0, color="k", linestyle="--", alpha=0.5)
      xlim(100,1200)
      #legend(loc=(0.36,-0.06), ncol=2, frameon=false)
      legend(loc=(0.00,0.06), ncol=2, frameon=false)

      #xlabel("simulation time [s]", x=0.5, y=0.5)
      plt.annotate("simulation time [s]", xy=(0.43, -0.36),
         xycoords="axes fraction")
      f3.minorticks_on()
      f3.tick_params(which="both",top=true, right=true)
      f3.tick_params(which="both", direction="in")
      plt.annotate("cc=$(round(rUBz,digits=2))", xy=(0.885, 0.05),
         xycoords="axes fraction")

      if occursin("PIC", filename)
         titletext = "MHD-EPIC, Location: "*
            string(data[nShift+1,11:end])[8:end]*L"R_G"
      else
         titletext = "Hall MHD, Location: "*
            string(data[nShift+1,11:end])[8:end]*L"R_G"
      end
      suptitle(titletext, x=0.5, y=0.94)
      f1.annotate("($('a'+iPlot-1))", xy=(0.0, 1.1), xycoords="axes fraction")
   end

   if DoPlot
      figure(figsize=(12.0,2.5))
      plt.rc("font", family="serif", size=12)
      plot(ΔB, alpha=0.95, label=L"\Delta B")
      plot(Δρ, alpha=0.95, label=L"\Delta \rho")
      xlim(0,1200)
      #ylim(-3,3)
      legend(loc="upper left", ncol=2, frameon=false)
      xlabel("simulation time [s]")
      plt.title("Location: "*string(data[nShift+1,11:end])[8:end]*L"R_G")
      tight_layout()
   end

   return rAlfven, rPBPt
end

function check_wave_type(filename="satellites_PIC.txt",
   dir="/Users/hyzhou/Documents/Computer/ParaView/data/")

   f = readdlm("/Users/hyzhou/Documents/Computer/ParaView/data/G8Traj.csv",
      ',', Float32, '\n'; header=true)
   nS = size(f[1])[1]

   # X, Y, Z, rAlfven, rFastSlow
   fnew = Array{Float32,2}(undef, nS, 5)

   header, data, satelliteNo = read_satellite_data(dir*filename)

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

## Execution
#static_location_plot("satellites_PIC.txt",
#   "/Users/hyzhou/Documents/Computer/ParaView/data/", 185)
#multi_satellite_plot()
#multi_satellite_contour("satellites_y0_PIC.txt", DoSubtractMean=true)
#multi_satellite_contour("satellites_boundary_PIC.txt", plane='z', DoSubtractMean=true)

peakUp, peakDn = satellite_p_contour("satellites_y0_PIC.txt"; No=2, plane='y')
#satellite_p_contour("satellites_boundary_Hall.txt"; No=3, plane='z')

#pMeanHall = satellite_p_contour_test("satellites_y0_Hall.txt"; No=1)
#pMeanPIC = satellite_p_contour_test("satellites_y0_PIC.txt"; No=1)

#nShift = 130
#static_location_plot("satellites_Hall.txt",
#   "/Users/hyzhou/Documents/Computer/ParaView/data/", nShift)
#wave_plot(nShift; DoPlot=true, filename="satellites_Hall.txt", iPlot=1,
#   verbose=true)
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
