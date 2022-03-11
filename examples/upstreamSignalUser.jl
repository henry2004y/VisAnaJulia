# User-defined upstream conditions.
#
# Example:
# t0 = 150 s
# period = 20 s
# n varies from 0.8 to 1.8 amu/cc
# Bz varies from 4.33 to 8.33 nT

using SpaceAnalysis

n0 = 1e6
T0 = 5e5
Vx0 = -750e3
Vy0 = 0.0
Vz0 = -1299e3
Bx0 = 2.5e-9
By0 = 0.0
Bz0 = 4.33e-9

tstart = 150.0
tend = 500.0
fsample = 10

len = floor(Int, (tend - tstart) * fsample) + 1
t = range(tstart, tend, length=len)

period = 20 # [s]
ϕ1 = asin(0.0) - 2π/period*tstart
ϕ2 = asin((4.33 - 5)/3) - 2π/period*tstart

n = @. (0.5*sin(2π/period*t+ϕ1) + 1.0)*1e6
T = fill(T0, len)
Vx, Vy, Vz = fill(Vx0, len), fill(Vy0, len), fill(Vz0, len)
Bx, By= fill(Bx0, len), fill(By0, len)
Bz = @. (3*sin(2π/period*t+ϕ2) + 5)*1e-9

s = BoundaryVariable(t, n, T, Vx, Vy, Vz, Bx, By, Bz)

save_signal(s; file="1dshockupstream.dat")