module Ganymede

export μ₀, amu, upstream_value, γ
export Index, getIndex

const μ₀ = 4π*1e-7 # Vacuum permeability, [H/m]
const amu = 1.66054e-27 # [kg]
# Rho Ux Uy Uz Bx By Bz Pe P, NaN for compensating the No index in the 1st col.
const upstream_value = [NaN, 56., 140., 0., 0., -10., -6., -86., 0.2, 3.4]
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
   Rho_= findfirst(x->x=="Rho", header) + 1
   P_  = findfirst(x->x=="P", header) + 1
   Pe_ = findfirst(x->x=="Pe", header) + 1
   Bx_ = findfirst(x->x=="Bx", header) + 1
   By_ = Bx_ + 1
   Bz_ = Bx_ + 2
   B_  = Bx_:(Bx_+2)
   Ux_ = findfirst(x->x=="Ux", header) + 1
   Uy_ = Ux_ + 1
   Uz_ = Ux_ + 2
   U_  = Ux_:(Ux_+2)

   id = Index(Rho_, Ux_, Uy_, Uz_, Bx_, By_, Bz_, Pe_, P_, B_, U_)
end


end
