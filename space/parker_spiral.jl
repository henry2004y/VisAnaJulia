# Parker spiral, heliospheric current sheet
# Reference:
# https://farside.ph.utexas.edu/teaching/plasma/lectures1/node68.html
# https://demonstrations.wolfram.com/TheInterplanetaryMagneticFieldParkerSpiral/
#
# Note: in the second link with Mathematica, many equations assume that θ=π/2
# so that sin(θ) term is missing!
#
# Hongyang Zhou, hyzhou@umich.edu

using PyPlot

include("constants.jl")

const Vr = 300e3 # [m/s]
const rmin = 1.0*Rsun # inner radius, [m]
const rmax = 5.3*AU # outer radius, [m]

# 2D ecliptic plane
φ₀ = 0
r = range(rmin, rmax, length=300)
φ = @. 2π *( φ₀ + rmin*ωsun/Vr*(r/rmin - log(r) - 1 + log(rmin)) )

fig = figure("Ecliptic plane",figsize=(10,10))
ax = PyPlot.axes(polar="true")
ax.plot(φ, r./AU, label="Parker spiral")

rCircle = ones(100)
angle = range(0, 2π, length=100)

ax.plot(angle, rCircle, "k--", label="Earth")
ax.plot(angle, rCircle.*dMercury, "--", label="Mercury")
ax.plot(angle, rCircle.*dVenus, "--", label="Venus")
ax.plot(angle, rCircle.*dMars, "--", label="Mars")
ax.plot(angle, rCircle.*dJupiter, "--", label="Jupiter")

ax.set_rmax(rmax/AU)
ax.set_rlabel_position(-22.5)  # Move radial labels away from plotted line
ax.grid(true)
legend()

ax.set_title("Interplanetary Magnetic Field In the Ecliptic Plane", va="bottom")

# 3D
nr, nθ = 10, 50
r_range = range(rmin, rmax, length=nr)
θ_range = range(π/16, π*15/16, length=nθ)

# meshgrid like operation
r = r_range' .* ones(nθ)
θ = ones(nr)' .* θ_range

fig = plt.figure("3D", figsize=(10,6))
using3D()
ax = fig.gca(projection="3d")

for i = 1:1
   φ₀ = i*π/4
   φ = @. 2π *( φ₀ + rmin*ωsun/Vr*(r/rmin - log(r) - 1 + log(rmin)) )

   X = @. r * sin(θ) * cos(φ) / AU
   Y = @. r * sin(θ) * sin(φ) / AU
   Z = @. r * cos(θ) / AU

   ax.scatter(X, Y, Z)
   #ax.plot_surface(X, Y, Z, alpha=0.3)
end

xlabel("x [AU]")
ylabel("y [AU]")
zlabel("z [AU]")