# Script for generating upstream MHD fluctuations.
#
# Hongyang Zhou, hyzhou@umich.edu

using PyPlot, SpaceAnalysis


## Background values

n0 = 2e6 # []
T0 = 5e5
V0 = [-600., 0.0, 0.0]*1e3 # [m/s] 
B0 = [0., 5., 0.]*1e-9 # [T]

varBG = BackgroundVariable(n0, T0, V0..., B0...)

fsignal = 1/150    # perturbation frequency, [Hz]
signal  = :density # perturbation type
dv      = 5e4      # perturbation magnitude in velocity, [m/s]
dn      = 1        # perturbation magnitude in density, [amu/cc]
dir     = "y"      # perturbed direction(s)
fsample = 2        # sampling frequency, [Hz]
tstart  = 400      # [s]
tend    = 1000     # [s]


s = generate_signal(varBG; fsignal, signal, dv, fsample, tstart, tend, dir)

save_signal(s)

## Visualization

fig, ax1 = plt.subplots(nrows=1, ncols=1)

if signal == :alfven
   color = "tab:red"
   ax1.set_xlabel("time [s]")
   ax1.set_ylabel("By [nT]"; color)
   ax1.plot(s.t, s.By*1e9; color)
   ax1.tick_params(axis="y", labelcolor=color)

   ax2 = ax1.twinx()
   color = "tab:blue"
   ax2.set_ylabel("Uy [km/s]"; color)
   ax2.plot(s.t, s.Vy/1e3, "--"; color)
   ax2.tick_params(axis="y", labelcolor=color)
elseif signal == :density
   color = "tab:red"
   ax1.set_xlabel("time [s]")
   ax1.set_ylabel("density [amu/cc]"; color)
   ax1.plot(s.t, s.n*1e-6; color)
   ax1.tick_params(axis="y", labelcolor=color)
elseif signal in (:slow, :fast)
   color = "tab:red"
   ax1.set_xlabel("time [s]")
   ax1.set_ylabel("By [nT]"; color)
   ax1.plot(s.t, s.By*1e9; color)
   ax1.tick_params(axis="y", labelcolor=color)

   ax2 = ax1.twinx()
   color = "tab:blue"
   ax2.set_ylabel("density [amu/cc]"; color)
   ax2.plot(s.t, s.n*1e-6; color)
   ax2.tick_params(axis="y", labelcolor=color)
end