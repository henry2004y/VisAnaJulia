# Shock calculations based on conversations.

using Roots

"""
    oblique_shock(Vx₁, Vy₁, Bx₁, By₁, T₁, rho₁; mass=1, ϵ = 0.1) ->
    (Vx₂, Vy₂, Bx₂, By₂, T₂, rho₂, P₂, X)

Calculates the Rankine-Hugoniot jump conditions on the other side of the shock, where `P₂`
is pressure and `X` is the compression ratio. `mass` is the atomic mass of plasma, and `ϵ`
is the control factor for determining the uniqueness of the solution.
# Arguments
- `Vx₁`: velocity component parallel to the shock normal direction.
- `Vy₁`: Velocity component perpendicular to the shock normal direction.
- `Bx₁`: Magnetic field component parallel to the shock normal direction.
- `By1`: Magnetic field component perpendicular to the shock normal direction.
- `T₁`: Temperature.
- `rho₁`: Density.
"""
function oblique_shock(Vx₁, Vy₁, Bx₁, By₁, T₁, rho₁)

   γ = 5 / 3 # adiabatic index
   m = mp * mass

   # Calculate other variables
   θ = acos(Bx₁ / √(Bx₁^2 + By₁^2))   # field angle w.r.t. normal dir
   vA₁ = √( Bx₁^2 / (m * rho₁ * μ₀) ) # Alfven speed
   V₁ = √( Vx₁^2 + Vy₁^2 )            # velocity
   P₁ = rho₁ * kB * T₁                # pressure
   vs₁ = √( γ * kB * T₁ / m )         # thermal speed

   # Function for solving the compression ratio `x`
   f(x) =
      (V₁^2 - x * vA₁^2)^2 * (x * vs₁^2 + 0.5V₁^2 * cos(θ)^2 * (x * (γ - 1) - (γ + 1))) +
      0.5vA₁^2 * V₁^2 * sin(θ)^2 * x *
      ((γ + x * (2 - γ)) * V₁^2 - x * vA₁^2 * ((γ + 1) - x * (γ - 1)))

   solutions = Float64[]
   # Find the first root of f(x)=0 starting at x=0.1
   root = find_zero(f, 0.1)

   push!(solutions, root)

   # Test the uniqueness of solution
   ϵ = 0.1 # control factor

   test_values = range(0.1, 10.0, step=0.1)
   for value in test_values
      root = find_zero(f, value)
      if !(root - ϵ ≤ result ≤ root + ϵ) && (ϵ ≥ solver(result) ≥ -ϵ)
         push!(solutions, root)
      end
   end

   length(solutions) > 1 || error("MORE THAN ONE SOLUTION: $solutions")

   X = solutions[1]
   rho₂ = rho₁ * X
   Vx₂ = Vx₁ / X
   Vy₂ = Vy₁ * (V₁^2 - vA₁^2) / (V₁^2 - X * vA₁^2)

   V₂ = √( Vx₂^2 + Vy₂^2 )

   Bx₂ = Bx₁
   By₂ = By₁ * (V₁^2 - vA₁^2) * X / (V₁^2 - X * vA₁^2)
   P₂ = P₁ * (X + (γ - 1) * X * V₁^2 * (1 - V₂^2 / V₁^2) / (2.0 * vs₁^2))
   T₂ = P₂ / (rho₂ * kB)

   return Vx₂, Vy₂, Bx₂, By₂, T₂, rho₂, P₂, X
end