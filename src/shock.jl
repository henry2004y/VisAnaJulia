# Shock calculations based on MHD conservations.

using LinearAlgebra
using NonlinearSolve
using Vlasiator: μ₀, kB, mᵢ

struct State{T<:Real}
   n::T
   V::Vector{T}
   B::Vector{T}
   T::T
   p::T
end

"https://farside.ph.utexas.edu/teaching/plasma/Plasma/img2921.png"
function shock_adiabatic_eq(r, θ, V₁², VA₁², VS₁², γ=5 / 3)
   cos2 = cos(θ)^2
   sin2 = sin(θ)^2

   term1 = ((V₁² - r * VA₁²)^2) * (((γ + 1) - r * (γ - 1)) * V₁² * cos2 - 2 * r * VS₁²)
   term2 = -r * sin2 * V₁² * VA₁² * (((γ + r * (2 - γ)) * V₁²) - ((γ + 1) - r * (γ - 1)) * r * VA₁²)

   return term1 + term2
end

"Calculates the compression ratio of the perpendicular shock."
function get_compression_perp_shock(M₁², β₁, γ)
   a = 2 * (2 - γ)
   b = γ * (2 * (1 + β₁) + (γ - 1) * β₁ * M₁²)
   c = -γ * (γ + 1) * β₁ * M₁²

   return (-b + √(b^2 - 4*a*c)) / (2*a)
end

"Reference: https://henry2004y.github.io/KeyNotes/contents/shock.html#oblique-shock"
function transform_to_dht_frame(V, B, n̂, doCheck::Bool=true)
   V_HT = V[1] / B[1] .* B
   VHT = V .- V_HT
   B_HT = B
   θ = acos(dot(V_HT, n̂) / norm(V_HT))

   if doCheck
      Emotional = -cross(V_HT, B_HT)
      @info "VHT = $VHT"
      @info "Check: motional E-field in HT frame: $Emotional"
   end

   return θ, V_HT, B_HT, VHT
end

"""
    get_downstream_oblique_shock(n₁, V₁, B₁, p₁, θ, VA₁², VS₁², n̂, γ) -> (r, p₂, n₂, V₂, B₂)

Given the upstream states, calculate the downstream states of the shock using
Rankine-Hugoniot jump conditions. `r` is the compression ratio. We assume proton here.
# Arguments
- `n₁`: Upstream number density.
- `V₁`: Upstream velocity vector.
- `B₁`: Upstream magnetic field vector.
- `p₁`: Upstream total thermal pressure.
- `θ`: Angle between the upstream plasma flow and the shock normal.
- `VA₁²`: Upstream Alfven speed squared.
- `VS₁²`: Upstream sonic speed squared.
- `n̂`: Shock normal.
- `γ`: adiabatic index.
# Reference
[Oblique MHD Shocks](https://farside.ph.utexas.edu/teaching/plasma/Plasma/node105.html)
"""
function get_downstream_oblique_shock(n₁, V₁, B₁, p₁, θ, VA₁², VS₁², n̂, γ)
   # Obtain normal and tangential components
   Vn1 = dot(V₁, n̂)
   Vt1_v = V₁ - Vn1 * n̂
   Vt1 = norm(Vt1_v)
   Bn1 = dot(B₁, n̂)
   Bt1_v = B₁ - Bn1 * n̂
   Bt1 = norm(Bt1_v)
   if dot(Vt1_v, Bt1_v) < 0
      Bt1 = -Bt1
   end

   V₁² = sum(V₁ .^ 2)

   # Compression ratio
   f(x, p) = shock_adiabatic_eq(x, p...)
   uspan = (1.0, (γ + 1) / (γ - 1))
   p = (θ, V₁², VA₁², VS₁²)
   prob_int = IntervalNonlinearProblem(f, uspan, p)
   r = solve(prob_int).u

   # Downstream calculations
   tmp = (V₁² - VA₁²) / (V₁² - r * VA₁²)
   Vn2 = Vn1 / r
   Vt2 = Vt1 * tmp

   R = 1 + γ * V₁² * (r - 1) / (VS₁² * r) * (
      cos(θ)^2 -
      0.5 * r * VA₁² * sin(θ)^2 * ((r + 1) * V₁² - 2 * r * VA₁²) / (V₁² - r * VA₁²)^2)
   p₂ = p₁ * R
   n₂ = n₁ * r

   Bn2 = Bn1
   Bt2 = Bt1 * r * tmp

   t̂ = V₁ - Vn1 * n̂
   t_mag = norm(t̂)
   if t_mag > 1e-18
      t̂ = t̂ / t_mag
   else
      t̂ = zeros(3)
   end

   B₂ = Bn2 * n̂ + Bt2 * t̂
   V₂ = Vn2 * n̂ + Vt2 * t̂

   @info "de Hoffman-Teller check:"
   @info "shock normal = $n̂"
   @info "V₁ x B₁ = $(cross(V₁, B₁))"
   @info "V₁ x B₂ = $(cross(V₂, B₂))"

   return r, p₂, n₂, V₂, B₂
end

"""
    get_downstream_perpendicular_shock(n₁, V₁, B₁, p₁, VA₁², VS₁², γ) -> (r, p₂, n₂, V₂, B₂)

# Reference
[Perpendicular MHD Shocks](https://farside.ph.utexas.edu/teaching/plasma/Plasma/node104.html)
"""
function get_downstream_perpendicular_shock(n₁, V₁, B₁, p₁, VA₁², VS₁², γ)
   V₁² = sum(V₁ .^ 2)
   M₁² = V₁² / VS₁²
   β₁ = 2 / γ * VS₁² / VA₁²
   r = get_compression_perp_shock(M₁², β₁, γ)
   n₂ = n₁ * r
   V₂ = V₁ / r
   B₂ = B₁ * r
   R = 1 + γ * M₁² * (1 - 1 / r) + (1 - r^2) / β₁
   p₂ = p₁ * R

   return r, p₂, n₂, V₂, B₂
end

"""
    get_downstream_parallel_shock(n₁, V₁, B₁, p₁, VS₁², γ) -> (r, p₂, n₂, V₂, B₂)

# Reference
[Parallel MHD Shocks](https://farside.ph.utexas.edu/teaching/plasma/Plasma/node103.html)
"""
function get_downstream_parallel_shock(n₁, V₁, B₁, p₁, VS₁², γ)
   V₁² = sum(V₁ .^ 2)
   M₁² = V₁² / VS₁²
   r = (γ + 1) * M₁² / (2 + (γ - 1) * M₁²)
   R = 1 + γ * M₁² * (1 - 1 / r)
   n₂ = n₁ * r
   V₂ = V₁ / r
   B₂ = B₁
   p₂ = p₁ * R

   return r, p₂, n₂, V₂, B₂
end


function show(n, V, B, T, p, n̂, location::String, pi2pe, γ)
   VS² = γ * p / (mᵢ * n)
   VA² = sum(B .^ 2) / (μ₀ * mᵢ * n)

   println("------------------------")
   println("$location state")
   println("------------------------")
   @info "VA [km/s]: $(sqrt(VA²) / 1e3)"
   @info "VS [km/s]: $(sqrt(VS²) / 1e3)"

   Vn² = dot(V, n̂)^2
   println("MA: ", sqrt(Vn² / VA²))

   Bn = dot(B, n̂)
   Bt = norm(B - Bn * n̂)

   Pram = n * mᵢ * Vn²
   Pb = 0.5 * sum(Bt .^ 2) / μ₀

   println("Ptotal [nPa]: ", (Pram + p + Pb) * 1e9)
   println("Pram [nPa]: ", Pram * 1e9)
   println("Pb [nPa]: ", Pb * 1e9)
   println("Pth [nPa]: ", p * 1e9)
   println("n [/cc]: ", n / 1e6)
   println("T [K]: ", T)
   println("Pi [nPa]: ", p * 1e9 / (1 + 1 / pi2pe))
   println("Pe [nPa]: ", p * 1e9 / (1 + pi2pe))
   println("V [km/s]: ", V / 1e3)
   println("B [nT]: ", B * 1e9)

   return
end


function show(state, n̂, location, pi2pe, γ)
   (; n, V, B, T, p) = state
   show(n, V, B, T, p, n̂, location, pi2pe, γ)
end


function check_pressure_balance(state1, state2, n̂)
   n₁, V₁, B₁, p₁ = state1.n, state1.V, state1.B, state1.p
   n₂, V₂, B₂, p₂ = state2.n, state2.V, state2.B, state2.p
   
   Bn1 = dot(B₁, n̂)
   Bt1 = norm(B₁ - Bn1 * n̂)
   Vn1² = dot(V₁, n̂)^2
   Pram1 = n₁ * mᵢ * Vn1²
   Pb1 = 0.5 * sum(Bt1 .^ 2) / μ₀

   Bn2 = dot(B₂, n̂)
   Bt2 = norm(B₂ - Bn2 * n̂)
   Vn2² = dot(V₂, n̂)^2
   Pram2 = n₂ * mᵢ * Vn2²
   Pb2 = 0.5 * sum(Bt2 .^ 2) / μ₀

   if abs((Pram1 + p₁ + Pb1) / (Pram2 + p₂ + Pb2) - 1) < 1e-6
      @info "Shock equilibrium verified!"
   else
      @warn "Shock is not at an equilibrium state!"
   end
end

function get_RH_downstream(MA₁, Btotal1, θ₁, T₁, n₁, pi2pe=1.0, γ=5/3)
   V₁ = [-MA₁ * Btotal1 / sqrt(μ₀ * mᵢ * n₁), 0, 0]

   B₁ = [Btotal1 * cos(θ₁), 0.0, Btotal1 * sin(θ₁)]
   n̂ = [1.0, 0.0, 0.0] |> normalize
   Vsh = 0.0

   V₁ .-= Vsh * n̂

   p₁ = n₁ * kB * T₁
   VS₁² = γ * p₁ / (mᵢ * n₁)
   VA₁² = sum(B₁ .^ 2) / (μ₀ * mᵢ * n₁)

   state1 = State(n₁, V₁, B₁, T₁, p₁)

   show(state1, n̂, "Upstream", pi2pe, γ)

   if θ₁ == 0.5 * π  # Perpendicular shock
      r, p₂, n₂, V₂, B₂ =
         get_downstream_perpendicular_shock(n₁, V₁, B₁, p₁, VA₁², VS₁², γ)
   elseif θ₁ == 0.0  # Parallel shock
      r, p₂, n₂, V₂, B₂ = get_downstream_parallel_shock(n₁, V₁, B₁, p₁, VS₁², γ)
   else  # Oblique shock
      θ, V₁HT, B₁HT, VHT = transform_to_dht_frame(V₁, B₁, n̂, false)
      r, p₂, n₂, V₂HT, B₂ =
         get_downstream_oblique_shock(n₁, V₁HT, B₁HT, p₁, θ, VA₁², VS₁², n̂, γ)

      V₂ = V₂HT + VHT + Vsh * n̂
   end

   T₂ = p₂ / (n₂ * kB)
   state2 = State(n₂, V₂, B₂, T₂, p₂)
   show(state2, n̂, "Downstream", pi2pe, γ)

   @info "Shock normal direction: $n̂"
   @info "Cone angle [degree]: $(θ₁ * 180 / π)"
   @info "Plasma compression ratio: $r"
   @info "Magnetic compression ratio: $(norm(B₂) / norm(B₁))"
   println("------------------------")
   check_pressure_balance(state1, state2, n̂)

   return
end


######

# Input upstream parameters
MA₁ = 5.0  # Mach number
Btotal1 = 5e-9  # [T]
θ₁ = 90 * π / 180   # In radians
T₁ = 1e5  # [K]
n₁ = 1e6  # [/m³]

get_RH_downstream(MA₁, Btotal1, θ₁, T₁, n₁)