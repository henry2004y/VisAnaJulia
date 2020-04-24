#
using PyPlot, DelimitedFiles, DSP, Statistics, FFTW

if !isdefined(Main, :Ganymede)
   include("/Users/hyzhou/Documents/Computer/Julia/BATSRUS/VisAnaJulia/space/Ganymede.jl")
   using .Ganymede: Rg, upstream_value
end

function smooth(x::Vector, n=100)
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

# Background
# G8 B= -10 -6 -86 [nT], U=140 0 0 [km/s]
Uxbk = upstream_value[2]
Bybk = upstream_value[6]
Bzbk = upstream_value[7]
TiltedAngle = atan(Bybk/Bzbk)

Δy = 5.64 # [Rg]

# Take clock angle in yz plane into consideration
Potential_bk = Δy * cos(TiltedAngle) * Rg *
   Uxbk*1e3 * abs(Bzbk*cos(TiltedAngle)+Bybk*sin(TiltedAngle))*1e-9 *
   1e-3 # [kV]

dir = "/Users/hyzhou/Documents/Computer/ParaView/data/"
dirHallBoundary = "dataBoundary_Hall/"
dirPICBoundary = "dataBoundary_PIC/"
dirHallSurface = "dataSurf_Hall/"
dirPICSurface = "dataSurf_PIC/"

data_hall = readdlm(dir*"CPCP_hall.txt", comments=true)
data_pic  = readdlm(dir*"CPCP_pic.txt", comments=true)

data_hall_b = readdlm(dirHallBoundary*"CPCP_test.txt", comments=true)
data_pic_b  = readdlm(dirPICBoundary*"CPCP_test.txt", comments=true)

data_hall_s = readdlm(dirHallSurface*"CPCP_test.txt", comments=true)
data_pic_s  = readdlm(dirPICSurface*"CPCP_test.txt", comments=true)

# Check missing frames
#missingtime = setdiff(0:1200,data_hall[:,1])
# hall: [355, 457, 505, 568]

#missingtime = setdiff(300:1500,data_pic[:,1])
# pic: [567, 591, 610, 680]

Fs = 1                 # Sampling frequency
T = 1/Fs               # Sampling period
L = size(data_hall_b)[1] # Length of signal
t = (0:L-1)*T          # Time vector

#missingtime = setdiff(0:1200,data_hall[:,1]) # Check missing frames

e_hall = data_hall[:,2] ./ Potential_bk
e_pic  = data_pic[:,2] ./ Potential_bk
e_hall_mean = mean(e_hall)
e_hall_std = std(e_hall)
e_pic_mean = mean(e_pic)
e_pic_std = std(e_pic)

e_hall_b = data_hall_b[:,2] ./ Potential_bk
e_pic_b  = data_pic_b[:,2] ./ Potential_bk
e_hall_mean_b = mean(e_hall_b)
e_hall_std_b = std(e_hall_b)
e_pic_mean_b = mean(e_pic_b)
e_pic_std_b = std(e_pic_b)


e_hall_s = data_hall_s[:,2] ./ Potential_bk
e_pic_s  = data_pic_s[:,2] ./ Potential_bk

#=
# Strange points: hall -> 13s, 14s, 415s, PIC -> 1255s, 1256s
# Remove strange points
for i = 2:length(e_hall_s)
   if e_hall_s[i] > 1
      e_hall_s[i] = 1.
   elseif e_hall_s[i] < 0
      e_hall_s[i] = 0.
   end
   if e_hall_s[i] - e_hall_s[i-1] > 0.2
      e_hall_s[i] = 0.5 * (e_hall_s[i] + e_hall_s[i-1])
   end
end


for i = 2:length(e_pic_s)
   if e_pic_s[i] > 1 || e_pic_s[i] < 0
      e_pic_s[i] = 0.0
   end
   if e_pic_s[i] - e_pic_s[i-1] > 0.2
      e_pic_s[i] = 0.5 * (e_pic_s[i] + e_pic_s[i-1])
   end
end
=#

e_hall_s = smooth(e_hall_s, 2)
e_pic_s = smooth(e_pic_s, 2)

e_hall_mean_s = mean(e_hall_s)
e_hall_std_s = std(e_hall_s)
e_pic_mean_s = mean(e_pic_s)
e_pic_std_s = std(e_pic_s)


fig, ax = plt.subplots(2,1,figsize=(10.0,5.5))
plt.rc("font", family = "serif", size = 13)

#ax[1].plot(data_hall[:,1], e_hall, label="Hall MHD Old")
ax[1].plot(data_hall_b[:,1], e_hall_b, label="Hall MHD")
ax[1].plot(data_hall_s[:,1], e_hall_s, alpha=0.8, label="Hall MHD S")

ax[1].legend()

#ax[2].plot(data_pic[:,1].+300, e_pic, label="MHD-EPIC Old")
ax[2].plot(data_pic_b[:,1], e_pic_b, label="MHD-EPIC")
ax[2].plot(data_pic_s[:,1], e_pic_s, alpha=0.8, label="MHD-EPIC S")

ax[2].legend()


Yhall = welch_pgram(data_hall_s[:,2], 180, 100; onesided=true, fs=Fs)
Ypic = welch_pgram(data_pic_s[:,2], 180, 100; onesided=true, fs=Fs)

Yhall_fft = fft(data_hall_s[:,2])
Ypic_fft = fft(data_pic_s[:,2])

P2_hall = @. abs(Yhall_fft/L)
P1_hall = P2_hall[1:floor(Int,L/2+1)]
P1_hall[2:end-1] .= 2 .* P1_hall[2:end-1]

f_hall = Fs*(0:(L/2))/L

P2_pic = @. abs(Ypic_fft/L)
P1_pic = P2_pic[1:floor(Int,L/2+1)]
P1_pic[2:end-1] .= 2 .* P1_pic[2:end-1]

f_pic = Fs*(0:(L/2))/L

fig, ax = plt.subplots(2,1,figsize=(10.0,5.5))
#plt.rc("font", family="serif", size=12)
ax[1].plot(1 ./ f_hall[2:end], P1_hall[2:end], label="Hall MHD")
ax[2].plot(1 ./ f_pic[2:end], P1_pic[2:end], label="MHD-EPIC")
ax[2].set_xlabel("Periods [s]")
#ax[2].set_ylabel("Power")
ax[1].grid("on")
ax[2].grid("on")
#plt.yscale("log")
#plt.xscale("log")
ax[2].legend()
ax[2].tick_params(which="both", direction="in")
ax[2].annotate("(b)", xy=(-0.03, 1.0), xycoords="axes fraction", fontsize=14)
tight_layout()


fig, ax = plt.subplots(2,1,figsize=(10.0,5.5))
#plt.rc("font", family="serif", size=12)
ax[1].plot(1 ./ Yhall.freq[2:end], Yhall.power[2:end], label="Hall MHD")
ax[2].plot(1 ./ Ypic.freq[2:end], Ypic.power[2:end], label="MHD-EPIC")
ax[2].set_xlabel("Periods [s]")
ax[2].set_ylabel("Power")
ax[1].grid("on")
ax[2].grid("on")
#plt.yscale("log")
#plt.xscale("log")
ax[2].legend()
ax[2].tick_params(which="both", direction="in")
ax[2].annotate("(b)", xy=(-0.03, 1.0), xycoords="axes fraction", fontsize=14)
tight_layout()



dirHallBoundary = "dataBoundary_Hall/"
dirPICBoundary = "dataBoundary_PIC/"

data_hall_b = readdlm(dirHallBoundary*"CPCP_test.txt", comments=true)
data_pic_b  = readdlm(dirPICBoundary*"CPCP_test.txt", comments=true)

data_test_hall = readdlm("dataSurf_Hall/CPCP_test.txt", comments=true)
data_test_pic = readdlm("dataSurf_PIC/CPCP_test.txt", comments=true)

figure()
plot(data_test_hall[:,1],data_test_hall[:,2], label="1:total of 2 and 3 ")
plot(data_test_hall[:,1],data_test_hall[:,3], label="2:time derivative of surface integral")
plot(data_test_hall[:,1],data_test_hall[:,4], label="3:middle line")
plot(data_hall_b[:,1], data_hall_b[:,2], label="4:upstream curve integral")
xlabel("time [s]")
ylabel("potential [kV]")
legend()
title("Hall MHD")
tight_layout()


figure()
plot(data_test_pic[:,1],data_test_pic[:,2], label="1:total of 2 and 3 ")
plot(data_test_pic[:,1],data_test_pic[:,3], label="2:time derivative of surface integral")
plot(data_test_pic[:,1],data_test_pic[:,4], label="3:middle line")
plot(data_pic_b[:,1], data_pic_b[:,2], label="4:upstream curve integral")
xlabel("time [s]")
ylabel("potential [kV]")
legend()
title("MHD-EPIC")
tight_layout()
