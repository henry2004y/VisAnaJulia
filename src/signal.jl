# Generating solarwind wave fluctuations.

export BackgroundVariable, BoundaryVariable
export generate_signal, save_signal

using DelimitedFiles, UnPack

abstract type Variable end

struct BackgroundVariable{T} <: Variable
   n::T
   T::T
   Vx::T
   Vy::T
   Vz::T
   Bx::T
   By::T
   Bz::T
end

struct BoundaryVariable{T} <: Variable
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
       signal=:alfven, dir="xyz", model=:vlasiator)

Generate a fluctuation signal as an upstream input file to plasma models.
Currently only support input perturbation in velocity `dv`.

# Optional Arguments
- `fsignal=1.0`: fluctuation signal frequency in [Hz].
- `fsample=5.0`: sampling frequency in [Hz].
- `tstart=0.0`: starting time in [s].
- `tend=10.0`: end time in [s].
- `dv=1e3`: velocity perturbation magnitude in [m/s].
- `dn=1`: mass density perturbation magnitude in [amu/cc].
- `signal=:alfven`: [:alfven, :fast].
- `dir="xyz"`: direction(s) for the applied Alfvénic perturbations in Cartesian coordinates,
which  can be any combination of "x", "y" and "z".
- `model=:vlasiator`: model specific settings. In Vlasiator, the background magnetic field
includes both a dipole and a constant. Therefore in the solarwind input files magnetic field
values should not contain the constant field.
"""
function generate_signal(varBG::BackgroundVariable; fsignal=1.0, fsample=5.0,
   tstart=0.0, tend=10.0, signal=:alfven, dir="xyz", dv=1e3, dn=1, model=:vlasiator)

   n0, T0, Vx0, Vy0, Vz0, Bx0, By0, Bz0 = varBG.n, varBG.T, varBG.Vx, varBG.Vy, varBG.Vz,
      varBG.Bx, varBG.By, varBG.Bz
   
   B0 = hypot(Bx0, By0, Bz0)
   μ₀ = 4π*1e-7
   mp = 1.6726219*1e-27
   VA = B0 / √(μ₀ * n0 * mp) # background Alfven speed

   len = floor(Int, (tend - tstart) * fsample) + 1
   t = range(tstart, tend, length=len)

   n = fill(n0, len)
   T = fill(T0, len)
   Vx, Vy, Vz = fill(Vx0, len), fill(Vy0, len), fill(Vz0, len)
   if model == :vlasiator
      Bx, By, Bz = zeros(len), zeros(len), zeros(len)
   else
      Bx, By, Bz = fill(Bx0, len), fill(By0, len), fill(Bz0, len)
   end

   if signal == :alfven
      δV = @. dv * sinpi(2 * fsignal * t)
      δB = @. B0 / VA * δV
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

   elseif signal == :fast
      δn = @. dn * 1e6 * sinpi(2 * fsignal * t)
      n .+= δn
   end

   v = BoundaryVariable(t, n, T, Vx, Vy, Vz, Bx, By, Bz)
end


"Save generated signal data `var` to `file`."
function save_signal(var::BoundaryVariable; file="sw.dat", sigdigits=5)
   @unpack t, n, T, Vx, Vy, Vz, Bx, By, Bz = var
   open(file, "w") do io
      write(io, "# time density temperature Vx Vy Vz Bx By Bz\n")
      writedlm(io, round.([t n T Vx Vy Vz Bx By Bz]; sigdigits))
   end
end