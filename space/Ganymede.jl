module Ganymede

export Rg, μ₀, amu, upstream_value, γ
export Index, getIndex

const Rg = 2634000.0 # radius, [m]
const μ₀ = 4π*1e-7 # Vacuum permeability, [H/m]
const amu = 1.66054e-27 # [kg]
# Rho Ux Uy Uz Bx By Bz Pe P, NaN for compensating the No index in the 1st col.
const upstream_value = [56., 140., 0., 0., -10., -6., -86., 0.2, 3.4]
const γ = 5/3 # adiabatic index

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
