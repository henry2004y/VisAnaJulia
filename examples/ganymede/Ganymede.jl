module Ganymede

using SpaceAnalysis: me, mp, ϵ0, q 

export Rg, μ₀, amu, upstream_value, γ
export Index, getIndex

const ionElectronMassRatio = 100 # ion-electron mass ratio in PIC simulation
const ionMass = 14 # single fluid ion mass in Ganymede's simulation
const Rg = 2634000. # Ganymede's radius, [m]
const dx = 1/100 # grid resolution for stream tracing, [Rg]
# Rho Ux Uy Uz Bx By Bz Pe P, NaN for compensating the No index in the 1st col.
const upstream_value = [56., 140., 0., 0., -10., -6., -86., 0.2, 3.4]
const γ = 5/3 # adiabatic index

const Z = 1 # Charge state
const ni = 4e6 # upstream ion number density for G8 flyby, [amu/m^3]
const ω_pi = √(ni*q^2/(ϵ0*ionMass*mp))*Z # ion plasma frequency, [rad/s]

# ion inertial length at Ganymede's upstream, [Rg]
const dᵢ = c / (√(ni*q^2/(ϵ0*ionMass*mp))*Z) / Rg

struct Index
   Rho_::Integer
   Ux_::Integer
   Uy_::Integer
   Uz_::Integer
   Bx_::Integer
   By_::Integer
   Bz_::Integer
   Pe_::Integer
   P_::Integer
   B_::UnitRange{Integer}
   U_::UnitRange{Integer}
end

function getIndex(header::Vector{AbstractString})
   Rho_= findfirst(x->x=="Rho", header)
   P_  = findfirst(x->x=="P", header)
   Pe_ = findfirst(x->x=="Pe", header)
   Bx_ = findfirst(x->x=="Bx", header)
   By_ = Bx_ + 1
   Bz_ = Bx_ + 2
   B_  = Bx_:(Bx_+2)
   Ux_ = findfirst(x->x=="Ux", header)
   Uy_ = Ux_ + 1
   Uz_ = Ux_ + 2
   U_  = Ux_:(Ux_+2)

   id = Index(Rho_, Ux_, Uy_, Uz_, Bx_, By_, Bz_, Pe_, P_, B_, U_)
end


end
