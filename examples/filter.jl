# Filters provided by DSP.jl

using PyPlot, DSP

# generate data
sin_noise(x) = sin.(x) .+ rand(length(x))
t = range(0, 30, length=6001) # time, [s]

fs = 1 / (t[2] - t[1])

y = sin_noise(t)

# create filter
responsetype = Lowpass(2; fs)
#responsetype = Bandpass(1, 5; fs)
designmethod = FIRWindow(hanning(64))
y_filtered = filt(digitalfilter(responsetype, designmethod), y)

plot(t, y)
plot(t, y_filtered)