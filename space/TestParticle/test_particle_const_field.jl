# Particle motion in constant EM field.

using Plots
pyplot()
using DifferentialEquations
using LinearAlgebra # cross product

const qₑ = -1.60217662e-19
const mₑ = 9.10938356e-31
const qᵢ = 1.60217662e-19
const mᵢ = 1.673557546e-27

function derivatives!(dy,y,p,t)
    q, m, E, B = p
    dy[1:3] = y[4:6]
    dy[4:6] = q/m*(E + y[4:6] × B)
end

x0 = [0.0, 0.0, 0.0]
u0 = [1.0, 0.0, 0.0]

y0 = [x0..., u0...] # initial position, velocity

E = [0.0,1.0,0.0]*1e-9  #[V/m]
B = [0.0,0.0,10.0]*1e-9 #[T]

# numerical issues with electron scale. Try to normalize the calculation!
p = [qₑ, mₑ, E, B] # charge, mass, E, B
#p = [qᵢ, mᵢ, E, B] # charge, mass, E, B
#tspan = (0.0,10.0)
tspan = (0.0,0.1)
prob = ODEProblem(derivatives!,y0,tspan,p)

sol = solve(prob)

#xyzt = plot(sol, lw=1.5)

#xy = plot(sol, vars=(1,2), xlabel="x", ylabel="y")
#xz = plot(sol, vars=(1,3), xlabel="x", ylabel="z")
#yz = plot(sol, vars=(2,3), xlabel="y", ylabel="z")
xyz = plot(sol, vars=(1,2,3), xlabel="x", ylabel="y", zlabel="z")
#plot(plot(xyzt,xyz),plot(xy, xz, yz, layout=(1,3),w=1), layout=(2,1))

#plot(sol.t,sol.u[1])
