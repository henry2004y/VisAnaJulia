# Final part for the global upstream reconnection efficiency estimation.
# Procedures:
# 1. Read in the processed reconnection rate data from Hall MHD and MHD-EPIC
# 2. Visualize the time series
# 3. Periodogram analysis with FFT, Welch method and Savitzky-Golay filter.
#
# Note that in the common log scale, the FFT doesn't show any significant peaks!
#
# Hongyang Zhou, hyzhou@umich.edu 04/14/2020

using FFTW, PyPlot, DelimitedFiles
using DSP, Statistics

if !isdefined(Main, :Ganymede)
   include("Ganymede.jl")
   using .Ganymede: Rg, upstream_value
end

include("Utility.jl")

# Background
# G8 B= -10 -6 -86 [nT], U=140 0 0 [km/s]
Uxbk = upstream_value[2]
Bybk = upstream_value[6]
Bzbk = upstream_value[7]
TiltedAngle = atan(Bybk/Bzbk)

#dir = "/Users/hyzhou/Documents/Computer/ParaView/data/"
dirHallSurface = "/Users/hyzhou/Documents/Computer/ParaView/scripts/ReconnectionRate/dataSurf_Hall/"
dirPICSurface = "/Users/hyzhou/Documents/Computer/ParaView/scripts/ReconnectionRate/dataSurf_PIC/"

##
# Get upstream theoretical maximum potential drop
Δy = 4.1 # or 5.64 [Rg], from z=2 cut estimation

# Take clock angle in yz plane into consideration, 179 [kV]
Potential_bk = Δy * cos(TiltedAngle) * Rg *
   Uxbk*1e3 * abs(Bzbk*cos(TiltedAngle)+Bybk*sin(TiltedAngle))*1e-9 *
   1e-3 # [kV]

#data_hall = readdlm(dir*"CPCP_hall.txt", comments=true)
#data_pic  = readdlm(dir*"CPCP_pic.txt", comments=true)
data_hall = readdlm(dirHallSurface*"CPCP_test.txt", comments=true)
data_pic  = readdlm(dirPICSurface*"CPCP_test.txt", comments=true)
data_pic[:,1] .-= 300.


if size(data_hall) != size(data_pic)
   @warn "Hall ($(size(data_hall)[1])) and MHD-EPIC ($(size(data_pic)[1])) " *
   "have different snapshots!"
end

Fs = 1                 # Sampling frequency
T = 1/Fs               # Sampling period
L = size(data_hall)[1] # Length of signal
t = (0:L-1)*T          # Time vector

#missingtime = setdiff(0:1200,data_hall[:,1]) # Check missing frames

e_hall = data_hall[:,2] ./ Potential_bk
e_pic  = data_pic[:,2] ./ Potential_bk

r_hall_mean = mean(data_hall[101:end,2])
r_hall_std = std(data_hall[101:end,2])
r_pic_mean = mean(data_pic[:,2])
r_pic_std = std(data_pic[:,2])

hall_mean = fill(r_hall_mean, size(e_hall))
pic_mean = fill(r_pic_mean, size(e_pic))

##
fig, ax = plt.subplots(3,1,figsize=(10.0,8.0))
plt.rc("font", family = "serif", size = 13)

ax[1].plot(data_hall[:,1], data_hall[:,2], label="Hall MHD")
ax[1].plot(data_hall[:,1], hall_mean, "k--", linewidth=0.8, label="Hall mean")
ax[1].plot(data_hall[end,1]:data_hall[end,1]+19, fill(r_hall_mean, 20), "k--",
   linewidth=0.8)
ax[1].errorbar(data_hall[end,1]+20., r_hall_mean, yerr=r_hall_std, ecolor="C0",
   capsize=3)

ax[1].set_xlim(-10.,1230.)
ax[1].set_ylim(30,135)
ax[1].set_xlabel("Simulation time [s]")
ax[1].set_ylabel("Reconnection rate [kV]")
ax[1].legend(loc="lower left", bbox_to_anchor=(0.07, 0.9), ncol=2,
   frameon=false)

ax[1].minorticks_on()
ax[1].tick_params(which="both", direction="in")
ax[1].annotate("(a)", xy=(-0.03, 1.01), xycoords="axes fraction", fontsize=14)

ax[2].plot(data_pic[:,1], data_pic[:,2], "C1", alpha=1.0, label="MHD-EPIC")
ax[2].plot(data_pic[:,1], pic_mean, "k--", linewidth=0.8, label="MHD-EPIC mean")
ax[2].plot(data_pic[end,1]:data_pic[end,1]+19, fill(r_pic_mean, 20), "k--",
   linewidth=0.8)
ax[2].errorbar(data_pic[end,1]+20, r_pic_mean, yerr=r_pic_std, ecolor="C1",
   capsize=3)

ax[2].set_xlim(-10.,1230.)
ax[2].set_ylim(35,130)
ax[2].set_xlabel("Simulation time [s]")
ax[2].set_ylabel("Reconnection rate [kV]")
ax[2].legend(loc="lower left", bbox_to_anchor=(0.07, 0.9), ncol=2,
      frameon=false)

ax[2].minorticks_on()
ax[2].tick_params(which="both", direction="in")
ax[2].annotate("(b)", xy=(-0.03, 1.01), xycoords="axes fraction", fontsize=14)

# FTE
#=
peak_pic = [34, 128, 230, 278, 347, 385, 434, 473, 656, 716, 727,
1016, 1083, 1116, 1151, 1181] .+ 10

peak_hall = [112, 153, 210, 352, 385, 412, 450, 462, 488, 515, 564, 636,
686, 755, 780, 797, 886, 930, 998, 1066, 1087, 1160] .+ 15
=#

# Check downwards FTEs!!!

peak_hall = [115, 153, 213, 355, 384, 413, 450, 488, 514, 559, 580, 639,
683, 712, 756, 774, 784, 797, 886, 927, 998, 1067, 1090, 1163, 1190]

peak_hall_Dn = [13, 41, 59, 112, 135, 193, 300, 580, 696,
758, 810, 845, 929, 981, 1001, 1033, 1048, 1059, 1078, 1170]

#append!(peak_hall,[13, 41, 59, 112, 135, 193, 300, 580, 696,
#758, 810, 845, 929, 981, 1001, 1033, 1048, 1059, 1078, 1170])

peak_pic = [35, 126, 230, 278, 345, 385, 436, 473, 653, 706, 715,
1010, 1085, 1118, 1150, 1183]
peak_pic_Dn = [7, 247, 262, 313, 336, 410, 477, 484, 499, 552, 617, 990, 1001,
1029, 1039, 1053, 1067] # downwards peaks

#append!(peak_pic,[7, 247, 262, 313, 336, 410, 477, 484, 499, 552, 617, 990, 1001,
#1029, 1039, 1053, 1067])

for t in peak_hall
   #=
   if t < 600
      ax[1].plot(t, e_hall[t], linestyle="", marker="o",
         markerfacecolor="k", markeredgecolor="None")
   else # there are 3 missing frames in the Hall outputs around 900s
      ax[1].plot(t, e_hall[t-3], linestyle="", marker="o",
         markerfacecolor="k", markeredgecolor="None")
   end
   =#
   #ax[1].plot(data_hall[t], e_hall[t], linestyle="", marker="o",
   #   markerfacecolor="k", markeredgecolor="None")
   ax[1].axvline(x=t, ls="-.", linewidth=0.8, color="r")
end

for t in peak_hall_Dn
   ax[1].axvline(x=t, ls=":", linewidth=0.8, color="g")
end

for t in peak_pic
   #ax[1].plot(t+300, e_pic[t], linestyle="", marker="o",
   #   markerfacecolor="r", markeredgecolor="None")
   #ax[1].plot(data_pic[t], e_pic[t], linestyle="", marker="o",
   #   markerfacecolor="r", markeredgecolor="None")
   ax[2].axvline(x=t, ls="-.", linewidth=0.8, color="r")
end

for t in peak_pic_Dn
   ax[2].axvline(x=t, ls=":", linewidth=0.8, color="g")
end

#include("satelliteAnalysis.jl")
#peak_hall = satellite_p_contour("satellites_y0_Hall.txt"; No=1, plane='y')
#peak_pic  = satellite_p_contour("satellites_y0_PIC.txt"; No=2, plane='y')
#=
# Check missing frames
MissHall = [355, 457, 505, 568]
MissPIC = [567, 591, 610, 680]
tHall = Set(0:1200)
tPIC = Set(0:1200)
for value in MissHall
   delete!(tHall,value)
end
for value in MissPIC
   delete!(tPIC,value)
end
tHall = sort(collect(tHall))
tPIC = sort(collect(tPIC))

ax1 = ax[1].twinx()  # instantiate a second axes that shares the same x-axis

ax1.set_ylabel("p perturbation", color="C2")  # we already handled the x-label with ax1
ax1.plot(tHall, pMeanHall, color="C2")
ax1.tick_params(axis="y", labelcolor="C2")

ax2 = ax[2].twinx()  # instantiate a second axes that shares the same x-axis

ax2.set_ylabel("p perturbation", color="C2")  # we already handled the x-label with ax1
ax2.plot(tPIC, pMeanPIC, color="C2")
ax2.tick_params(axis="y", labelcolor="C2")
=#
##

#Yhall = welch_pgram(data_hall[121:end,2]; onesided=true, fs=Fs)
#Ypic = welch_pgram(data_pic[:,2]; onesided=true, fs=Fs)
Yhall = periodogram(data_hall[121:end,2], onesided=true, nfft=L-120, fs=Fs)
Ypic = periodogram(data_pic[:,2], onesided=true, nfft=L, fs=Fs)

#include("/Users/hyzhou/Documents/Computer/Julia/FFT/savitzkyGolay.jl")
#Yhall_fft = savitzky_golay_filter(data_hall[61:end,2],201,1)
#Ypic_fft  = savitzky_golay_filter(data_pic[:,2],201,1)
#=
Yhall_fft = fft(data_hall[121:end,2])
Ypic_fft = fft(data_pic[:,2])

P2_hall = @. abs(Yhall_fft/(L-120))
P1_hall = P2_hall[1:floor(Int,(L-120)/2+1)]
P1_hall[2:end-1] .= 2.0 .* P1_hall[2:end-1]

f_hall = Fs*(0:((L-120)/2))/(L-120)

P2_pic = @. abs(Ypic_fft/L)
P1_pic = P2_pic[1:floor(Int,L/2+1)]
P1_pic[2:end-1] .= 2.0 .* P1_pic[2:end-1]

f_pic = Fs*(0:(L/2))/L
=#
#ax[3].plot(1 ./ f_hall[2:end], P1_hall[2:end], label="Hall MHD")
#ax[3].plot(1 ./ f_pic[2:end], P1_pic[2:end], label="MHD-EPIC")
# Note: this is PSD in linear scale!
ax[3].plot(1 ./ Yhall.freq[1:end], Yhall.power[1:end], label="Hall MHD")
ax[3].plot(1 ./ Ypic.freq[1:end],  Ypic.power[1:end],  label="MHD-EPIC")
ax[3].set_xlabel("Period [s]")
ax[3].set_ylabel("PSD")
ax[3].grid("on")
ax[3].set_xlim(0,150)
#plt.yscale("log")
#plt.xscale("log")
ax[3].legend()
ax[3].minorticks_on()
ax[3].tick_params(which="both", direction="in")
ax[3].annotate("(c)", xy=(-0.03, 1.0), xycoords="axes fraction", fontsize=14)
tight_layout()

fig.savefig(
   "/Users/hyzhou/Documents/research/paper/Ganymede2/test figure/ReconnectionEfficiency_FFT_surfaceIntegral.png",
   dpi=150.0)

@info "Hall CPCP"
println("mean = ", mean(data_hall[:,2]), " ", r_hall_mean/Potential_bk)
println("std  = ", std(data_hall[:,2]), " ", r_hall_std/Potential_bk)
@info "PIC CPCP"
println("mean = ", mean(data_pic[:,2]), " ", r_pic_mean/Potential_bk)
println("std  = ", std(data_pic[:,2]), " ", r_pic_std/Potential_bk)
