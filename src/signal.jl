# Generating solarwind wave fluctuations.

export VariableBackground, VariableBoundary
export generate_signal, save_signal

using DelimitedFiles

abstract type Variable end

struct VariableBackground{T} <: Variable
   n::T
   T::T
   Vx::T
   Vy::T
   Vz::T
   Bx::T
   By::T
   Bz::T
end

struct VariableBoundary{T} <: Variable
   t::AbstractVector{T}
   n::Vector{T}
   T::Vector{T}
   Vx::Vector{T}
   Vy::Vector{T}
   Vz::Vector{T}
   Bx::Vector{T}
   By::Vector{T}
   Bz::Vector{T}
end

"""
    generate_signal(varBG; fsignal=1.0, dv=1e3, fsample=5.0, tstart=0.0, tend=10.0,
       signaltype=:alfven, dir="xyz", modetype=:vlasiator)

Generate a fluctuation signal of frequency `fsignal` based on `varBG` from `tstart` to
`tend` with sampling frequency `fsample`. For an Alfenic `signaltype`, `dir` is used to
specify the perturbed direction in Cartesian coordinates, which can be any combination of
"x", "y" and "z".
Currently only support input perturbation in velocity `dv`.
In Vlasiator, the background magnetic field includes both a dipole and a constant. Therefore
in the solarwind input files magnetic field values should not contain the constant field.
"""
function generate_signal(varBG; fsignal=1.0, dv=1e3, fsample=5.0, tstart=0.0, tend=10.0,
   signaltype=:alfven, dir="xyz", modeltype=:vlasiator)

   n0, T0, Vx0, Vy0, Vz0, Bx0, By0, Bz0 = varBG.n, varBG.T, varBG.Vx, varBG.Vy, varBG.Vz,
      varBG.Bx, varBG.By, varBG.Bz
   
   B0 = hypot(Bx0, By0, Bz0)
   μ₀ = 4π*1e-7
   mp = 1.6726219*1e-27
   VA = B0 / √(μ₀ * n0 * mp) # background Alfven speed

   len = floor(Int, tend * fsample)
   t = range(tstart, tend, length=len)

   δV = @. dv * sin(2π * fsignal * t)

   δB = @. B0 / VA * δV

   n = fill(n0, len)
   T = fill(T0, len)
   Vx, Vy, Vz = fill(Vx0, len), fill(Vy0, len), fill(Vz0, len)
   if modeltype == :vlasiator
      Bx, By, Bz = zeros(len), zeros(len), zeros(len)
   else
      Bx, By, Bz = fill(Bx0, len), fill(By0, len), fill(Bz0, len)
   end

   if signaltype == :alfven
      if occursin("x", dir)
         Vx .+= δV
         Bx .+= δB
      end
      if occursin("y", dir)
         Vy .+= δV
         By .+= δB
      end
      if occursin("z", dir)
         Vz .+= δV
         Bz .+= δB
      end

   elseif signaltype == :fast
      @error "Not yet implemented!"
   end

   v = VariableBoundary(t, n, T, Vx, Vy, Vz, Bx, By, Bz)
end


"Save generated signal data `var` to `file`."
function save_signal(var::VariableBoundary; file="sw.dat", sigdigits=5)
   t, n, T = var.t, var.n, var.T
   Vx, Vy, Vz, Bx, By, Bz = var.Vx, var.Vy, var.Vz, var.Bx, var.By, var.Bz
   open(file, "w") do io
      write(io, "# time density temperature Vx Vy Vz Bx By Bz\n")
      writedlm(io, round.([t n T Vx Vy Vz Bx By Bz]; sigdigits))
   end
end