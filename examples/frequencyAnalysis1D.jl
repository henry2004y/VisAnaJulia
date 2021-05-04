# Frequency analysis of 1D data.
#
#
# Hongyang Zhou, hyzhou@umich.edu 12/10/2019

using MAT, Statistics, PyPlot
using SpaceAnalysis

"""
    B_spectrum(filename, varname, vartype, Fs, Brange)

Return single-side spectrum amplitude for magnetic field read from MAT file
`filename`. `varname` ∈ ["Bobs", "Bsim"]. `vartype` ∈ ["x", "y", "z", "mag"].

Note: MAT package seems not able to read timeseries data. A workaround is by
identifying the starting and ending array index in magnetopause crossings.
"""
function B_spectrum(filename, varname="Bobs", vartype="mag", Fs=1/0.3,
   Brange=1500:6000)

   data = matread(filename)
   B = data[varname]
   if vartype == "mag"
      Btmp = hypot.(B[Brange,1], B[Brange,2], B[Brange,3])
   elseif vartype == "x"
      Btmp = B[Brange,1]
   elseif vartype == "y"
      Btmp = B[Brange,2]
   elseif vartype == "z"
      Btmp = B[Brange,3]
   end

   spectrum(Btmp, Fs)
end


filenames = ["B_G8.mat", "B_G8_PIC.mat", "B_G8_PIC_AMR2.mat", "B_G8_Hall_AMR2.mat"]

fObs, pObs   = B_spectrum(filenames[1], "Bobs", "mag")
fHall, pHall = B_spectrum(filenames[1], "Bsim", "mag", 1.0, 1573:2959)
fHall2, pHall2 = B_spectrum(filenames[4], "Bsim", "mag", 1.0, 1573:2959)
fPIC, pPIC   = B_spectrum(filenames[2], "Bsim", "mag", 1.0, 1573:2959)
fPIC2, pPIC2 = B_spectrum(filenames[3], "Bsim", "mag", 1.0, 1573:2959)

nSmooth = 19
nSHalf = (nSmooth-1) ÷ 2

figure()
plot(fObs[2:end], sma(pObs[2:end].^2, nSmooth),"k")
plot(fHall[2:end], sma(pHall[2:end].^2, nSmooth))
plot(fHall2[2:end], sma(pHall2[2:end].^2, nSmooth))
plot(fPIC[2:end], sma(pPIC[2:end].^2, nSmooth),"#d62728")
plot(fPIC2[2:end], sma(pPIC2[2:end].^2, nSmooth))

plt.yscale("log")
xlim(0.0, 0.3)
grid(true)
legend(["Observation", "Hall MHD fine grid", "Hall MHD coarse grid",
   "PIC fine grid", "PIC coarse grid"])
title("Single-Sided Amplitude Spectrum of B(t)")
xlabel("f (Hz)")
ylabel("|P1(f)|")