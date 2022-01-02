# Commonly used constants

const q  = 1.6021765e-19    # electron charge, [C]
const me = 9.10938356e-31   # electron mass, [kg]
const mp = 1.6726219e-27    # proton mass, [kg]
const qₑ = -q
const mₑ = me
const qᵢ = q
const mᵢ = mp
const amu= 1.66054e-27      # Atomic Mass Unit, 1/12 of C12, [kg]
const ϵ0 = 8.8542e-12       # vacuum permittivity, [F/m]
const ϵ₀ = ϵ0
const c  = 3e8              # speed of light, [m/s]
const μ₀ = 4π*1e-7          # Vacuum permeability, [H/m]
const kB = 1.38064852e-23   # Boltzmann constant, [m²kg/(s²K)]

const Rₑ = 6.38e6           # Earth radius, [m]
const Rm = 2444000.         # Mercury's radius, [m]

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
