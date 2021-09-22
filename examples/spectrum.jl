# Spectra and time series plots from virtual satellite data.
#
# t rho vx vy p bz ex ey
#
# Hongyang Zhou, hyzhou@umich.edu

using SignalAnalysis, DelimitedFiles, Glob
using SignalAnalysis.Units
using Plots, Measures

files = glob("satellite*.csv")

samplerate = 2Hz
color = :dense
left_margin = 5mm

const ax = fill(plot(),10,1)

for file in files
   println("file: ", file)
   data = readdlm(file, ','; header=true)
   t0 = data[1][1]
   p = signal(data[1][:,5] .* 1e9, samplerate) # [nPa]
   bz = signal(data[1][:,6] .* 1e9, samplerate)# [nT]
   ex = signal(data[1][:,7] .* 1e3, samplerate)# [mV/m]
   ey = signal(data[1][:,8] .* 1e3, samplerate)# [mV/m]

   l = @layout[[a1;a2;a3;a4;a5] grid(2,2){0.7w}]
   #ax[1] = plot(title="Power Spectra at [$(file[end-5:end-4]), 0, 0]Re",
   #   framestyle=nothing,showaxis=false,xticks=false,yticks=false,margin=0mm)

   x = parse(Int, file[end-5:end-4])
   ax[2] = plot(data[1][:,1], data[1][:,2] ./ 1e6; label="density", title="[$x, 0, 0]Re",
      ylabel="[amu/cc]", left_margin)
   ax[3] = plot(data[1][:,1], data[1][:,3:4] ./ 1e3; label=["Vx" "Vy"],
      ylabel="[km/s]", left_margin)
   ax[4] = plot(data[1][:,1], samples(p); label="p", ylabel="[nPa]", left_margin)
   ax[5] = plot(data[1][:,1], samples(bz); label="bz", ylabel="[nT]", left_margin)
   ax[6] = plot(data[1][:,1], data[1][:,7:8] .* 1e3; label=["Ex" "Ey"],
      ylabel="[mV/m]", xlabel="Time [s]", left_margin)

   # magnitude units in [db], or 10log10(power)
   ax[7] = specgram(p; nfft=128, t0, crange=100, color ,title="P")
   ax[8] = specgram(bz; nfft=128, t0, crange=100, color, title="Bz")
   ax[9] = specgram(ex; nfft=128, t0, crange=100, color, title="Ex")
   ax[10] = specgram(ey; nfft=128, t0, crange=100, color, title="Ey")
   plot(ax[2:10]..., layout=l, size=(1500,1000), legend=true)
   savefig("spectrum_"*file[end-5:end-4]*"Re.png")
end