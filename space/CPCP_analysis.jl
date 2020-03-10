using FFTW, PyPlot, DelimitedFiles
using DSP, Statistics

#include("/Users/hyzhou/Documents/Computer/Julia/FFT/savitskyGolay.jl")

if !isdefined(Main, :Ganymede)
   include("Ganymede.jl")
   using .Ganymede: Rg, upstream_value
end

# Background
# G8 B= -10 -6 -86 [nT], U=140 0 0 [km/s]
Uxbk = upstream_value[2]
Bybk = upstream_value[6]
Bzbk = upstream_value[7]
TiltedAngle = atan(Bybk/Bzbk)

dir = "/Users/hyzhou/Documents/Computer/ParaView/data/"

##
# Get upstream theoretical maximum potential drop
filename = dir*"boundary.txt"
bc = readdlm(filename, ',', header=true)
data = bc[1]

#Δy =  maximum(data[:,2]) - minimum(data[:,2])
Δy = 5.7 # [Rg]

# Take clock angle in yz plane into consideration
Potential_bk = Δy * cos(TiltedAngle) * Rg *
   Uxbk*1e3 * abs(Bzbk*cos(TiltedAngle)+Bybk*sin(TiltedAngle))*1e-9 *
   1e-3 # [kV]

Fs = 1           # Sampling frequency
T = 1/Fs         # Sampling period
L = 1197         # Length of signal
t = (0:L-1)*T    # Time vector

data_hall = readdlm(dir*"CPCP_hall.txt")
data_pic  = readdlm(dir*"CPCP_pic.txt")

# find the missing data and pad to them!!!

Yhall = welch_pgram(data_hall[:,2], 180, 100; onesided=true, fs=Fs)
Ypic = welch_pgram(data_pic[:,2], 200, 100; onesided=true, fs=Fs)

e_hall = data_hall[:,2] ./ Potential_bk
e_pic  = data_pic[:,2] ./ Potential_bk
e_hall_mean = mean(e_hall)
e_hall_std = std(e_hall)
e_pic_mean = mean(e_pic)
e_pic_std = std(e_pic)

hall_mean = fill(e_hall_mean, size(e_hall))
pic_mean = fill(e_pic_mean, size(e_pic))

##
fig, ax = plt.subplots(2,1,figsize=(10.0,5.5))
plt.rc("font", family = "serif", size = 13)

ax[1].plot(data_hall[:,1], e_hall, label="Hall MHD")
ax[1].plot(data_pic[:,1].+300, e_pic, alpha=0.8, label="MHD-EPIC")
ax[1].plot(data_hall[:,1], hall_mean, "C0--", linewidth=0.8, label="Hall mean")
ax[1].plot(-20.0:-1.0, fill(e_hall_mean, 20), "C0--", linewidth=0.8)
ax[1].errorbar(data_hall[1,1]-20., e_hall_mean, yerr=e_hall_std, ecolor="C0",
   capsize=3)
ax[1].plot(data_pic[:,1].+300, pic_mean, "C1--", linewidth=0.8, label="PIC mean")
ax[1].plot(data_pic[end,1]+300:data_pic[end,1]+319, fill(e_pic_mean, 20), "C1--",
   linewidth=0.8)
ax[1].errorbar(data_pic[end,1]+320., e_pic_mean, yerr=e_pic_std, ecolor="C1",
   capsize=3)

ax[1].set_xlim(-50.,1550.)
ax[1].set_xlabel("Simulation time [s]")
ax[1].set_ylabel("Reconnection efficiency")
ax[1].legend(loc="lower left", bbox_to_anchor=(0.07, 0.80), ncol=4,
   frameon=false)
#=
ax[1].annotate("Hall mean = $(round(e_hall_mean, digits=2))\n"*
   "Hall std = $(round(e_hall_std, digits=2))\n"*
   "PIC mean = $(round(e_pic_mean, digits=2))\n"*
   "PIC std = $(round(e_pic_std, digits=2))",
   xy=(1270, 0.535), xycoords="data")
=#
ax[1].minorticks_on()
ax[1].tick_params(which="both", direction="in")
ax[1].annotate("(a)", xy=(-0.03, 1.0), xycoords="axes fraction", fontsize=14)
tight_layout()

# FTE
peak_pic = [34, 128, 230, 278, 347, 385, 434, 473, 656, 716, 727,
1016, 1083, 1116, 1151, 1182] .+ 10

peak_hall = [112, 153, 210, 352, 385, 412, 450, 462, 488, 515, 564, 636,
686, 755, 780, 797, 886, 930, 998, 1066, 1087, 1160] .+ 15

for t in peak_pic
   ax[1].plot(t+300, e_pic[t], linestyle="", marker="o",
      markerfacecolor="r", markeredgecolor="None")
end

for t in peak_hall
   if t < 600
      ax[1].plot(t, e_hall[t], linestyle="", marker="o",
         markerfacecolor="k", markeredgecolor="None")
   else # there are 3 missing frames in the Hall outputs around 900s
      ax[1].plot(t, e_hall[t-3], linestyle="", marker="o",
         markerfacecolor="k", markeredgecolor="None")
   end
end

#include("satelliteAnalysis.jl")
#peak_hall = satellite_p_contour("satellites_y0_Hall.txt"; No=1, plane='y')
#peak_pic  = satellite_p_contour("satellites_y0_PIC.txt"; No=2, plane='y')

##
#fig = figure(figsize=(10.0,3.7))
#plt.rc("font", family="serif", size=12)
ax[2].plot(1 ./ Yhall.freq[2:end], Yhall.power[2:end], label="Hall MHD")
ax[2].plot(1 ./ Ypic.freq[2:end], Ypic.power[2:end], label="MHD-EPIC")
ax[2].set_xlabel("Periods [s]")
ax[2].set_ylabel("Power")
ax[2].grid("on")
#plt.yscale("log")
#plt.xscale("log")
ax[2].legend()
ax[2].tick_params(which="both", direction="in")
ax[2].annotate("(b)", xy=(-0.03, 1.0), xycoords="axes fraction", fontsize=14)
tight_layout()

fig.savefig(
   "/Users/hyzhou/Documents/research/paper/Ganymede2/test figure/CPCP_FFT.png",
   dpi=150.0)

@info "Hall CPCP"
println("mean = ", mean(data_hall[:,2]))
println("std  = ", std(data_hall[:,2]))
@info "PIC CPCP"
println("mean = ", mean(data_pic[:,2]))
println("std  = ", std(data_pic[:,2]))
