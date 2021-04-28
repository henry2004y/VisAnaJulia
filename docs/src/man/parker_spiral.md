# Interplanetary Magnetic Field

The [interplanetary magnetic field](https://en.wikipedia.org/wiki/Interplanetary_magnetic_field) (IMF), also called Parker spiral, is the component of the solar magnetic field that is dragged out from the solar corona by the solar wind flow to fill the Solar System. Depending on the polarity of the photospheric footpoint, the heliospheric magnetic field spirals inward or outward; the magnetic field follows the same shape of spiral in the northern and southern parts of the heliosphere, but with opposite field direction. These two magnetic domains are separated by a current sheet (an electric current that is confined to a curved plane). This heliospheric current sheet has a shape similar to a twirled ballerina skirt, and changes in shape through the solar cycle as the Sun's magnetic field reverses about every 11 years. 

References:
[Richard Fitzpatrick](https://farside.ph.utexas.edu/teaching/plasma/lectures1/node68.html)
[Mathematica](https://demonstrations.wolfram.com/TheInterplanetaryMagneticFieldParkerSpiral/)[^1]

[^1]: many equations in the Mathematic page assume that θ=π/2, so that sin(θ) term is missing!

```@example 1
using PyPlot

# Solar
const Tsolarrotate = 25.05*24*3600 # [s]
const ωsun = 2π / Tsolarrotate # [rad/s]
const Rsun = 696340e3 # solar radius, [m]
const AU = 1.496e11 # sun-Earth distance, [m]

# Distance from sun
# https://www.enchantedlearning.com/subjects/astronomy/planets/
const dMercury = 0.39  # [AU], 57.91e9 [m]
const dVenus   = 0.723 # [AU], 108.2e9 [m]
const dMars    = 1.524 # [AU], 227.9e9 [m]
const dJupiter = 5.203 # [AU], 778.5e9 [m]
const dSaturn  = 9.539 # [AU], 1434e9  [m]
const dUranus  = 19.18 # [AU], 2871e9  [m]
const dNeptune = 30.06 # [AU], 4495e9  [m]

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
```

```@example 1
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
```