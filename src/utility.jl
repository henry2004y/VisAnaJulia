# Utility functions for space data analysis.

using ImageFiltering, OffsetArrays, FFTW, Statistics, LinearAlgebra

export mag2db, MVA, sma, spectrum

"Convert `x` from magnitude to decibels."
mag2db(x) = 10. * log10.(x)

"""
    MVA(Bx, By, Bz; verbose=false)

Perform minimum variance analysis to vector components defined in orthogonal
coordinates `Bx`, `By` and `Bz`.
If λ₁ ≥ λ₂ ≥ λ₃ are 3 eigenvalues of the constructed matrix M, then a good
indicator of nice fitting LMN coordinate system should have λ₂/λ₃ > 5. Set 
`verbose=true` to turn on the check.
"""
function MVA(Bx, By, Bz; verbose=false)

   B̄1 = mean(Bx)
   B̄2 = mean(By)
   B̄3 = mean(Bz)
   B̄11= mean(Bx.*Bx) - B̄1*B̄1
   B̄22= mean(By.*By) - B̄2*B̄2
   B̄33= mean(Bz.*Bz) - B̄3*B̄3
   B̄12= mean(Bx.*By) - B̄1*B̄2
   B̄23= mean(By.*Bz) - B̄2*B̄3
   B̄31= mean(Bz.*Bx) - B̄3*B̄1
   # Construct the matrix
   M = [B̄11 B̄12 B̄31; B̄12 B̄22 B̄23; B̄31 B̄23 B̄33]

   # Compute the eigen values and ratios (descending order)
   F = eigen(M, sortby = x -> -abs(x))

   if verbose
      println("Eigenvalues:", F.values)
      println("Eigenvectors:")
      println(F.vectors[:,1])
      println(F.vectors[:,2])
      println(F.vectors[:,3])
      r = F.values[2] / F.values[3]
      println("Ratio of intermediate variance to minimum variance = ",
         round(r, digits=3))
      if r ≥ 5
         println("Seems to be a proper MVA attempt!")
      else
         @warn "Take the MVA result with a grain of salt!"
      end
   end
   F
end

"""
    MFA(x, y, z, Bx, By, Bz; method="MAVG", nw=15)

Obtain the time series of rotation matrix from geocentric solar ecliptic (GSE)
to Mean field-aligned (MFA) reference frame, given the positions of the
spacecraft `x`, `y`, `z` and magnetic fields `Bx`, `By`, `Bz` in GSE.
We can choose `method` from moving average `"MAVG"` with `nw` points, or 
empirical mode decomposition `"EMD"`.

Reference: [Regi+, 2016] https://doi.org/10.4401/ag-7067
"""
function MFA(x::AbstractVector, y::AbstractVector, z::AbstractVector,
   Bx::AbstractVector, By::AbstractVector, Bz::AbstractVector;
   method="MAVG", nw=15)
   if !(length(x) == length(y) == length(z) ==
      length(Bx) == length(By) == length(Bz))
      throw(DimensionMismatch("requires same vector length for inputs!"))
   end

   μx, μy, μz = similar(x), similar(x), similar(x)
   ϕx, ϕy, ϕz = similar(x), similar(x), similar(x)
   νx, νy, νz = similar(x), similar(x), similar(x)
   Bx′, By′, Bz′ = similar(x), similar(x), similar(x)

   if method == "MAVG"
      B̄x = sma(Bx, nw)
      B̄y = sma(By, nw)
      B̄z = sma(Bz, nw)

      for i = axes(x)
         B̄mag = hypot(B̄x[i], B̄y[i], B̄z[i])
         μx[i] = B̄x[i] / B̄mag
         μy[i] = B̄y[i] / B̄mag
         μz[i] = B̄z[i] / B̄mag
         ϕx[i] = y[i]*μz[i] - z[i]*μy[i]
         ϕy[i] = z[i]*μx[i] - x[i]*μz[i]
         ϕz[i] = x[i]*μy[i] - y[i]*μx[i]
         rBmag = hypot(ϕx[i], ϕy[i], ϕz[i])
         ϕx[i] /= rBmag
         ϕy[i] /= rBmag
         ϕz[i] /= rBmag
         νx[i] = μy[i]*ϕz[i] - μz[i]*ϕy[i]
         νy[i] = μz[i]*ϕx[i] - μx[i]*ϕz[i]
         νz[i] = μx[i]*ϕy[i] - μy[i]*ϕx[i]

         Bx′[i] = μx[i]*Bx[i] + μy[i]*By[i] + μz[i]*Bz[i]
         By′[i] = ϕx[i]*Bx[i] + ϕy[i]*By[i] + ϕz[i]*Bz[i]
         Bz′[i] = νx[i]*Bx[i] + νy[i]*By[i] + νz[i]*Bz[i]
      end
   elseif method == "EMD"

   end
   Bx′, By′, Bz′
end

"""
    sma(x, n=(5,5))

Simple moving average of `x` with `n[1]` leading and `n[2]` trailing points.
Default boundary is "replicate", meaning that the border pixels extend beyond
the boundaries.
"""
function sma(x, n::Tuple{Int,Int}=(5,5))
   len = sum(n) + 1
   kernel = OffsetArray(fill(1/len, len), -n[1]:n[2])
   imfilter(x, kernel)
end

"""
    sma(x, n=5)

Simple moving box average of the vector data `x` with box length 'n'.
One-sided average on the left and right edge with replicate border.
"""
function sma(x::Vector, n::Int=5)
   iseven(n) && @warn "Even box length detected!"
   sma(x, ((n-1)÷2, (n-1)÷2))
end

"""
    spectrum(x, Fs)

Return the frequency and amplitude for the single-sided spectrum of vector `x`
given sample frequency of `Fs`.
"""
function spectrum(x::Vector, Fs)
   L = length(x)
   x̃ = fft(x)
   # Compute the two-sided spectrum amplitude P2.
   P₂ = @. abs(x̃/L)
   # Compute the single-sided spectrum amplitude P1.
   P₁ = P₂[1:floor(Int, L/2+1)]
   P₁[2:end-1] .= 2 .* P₁[2:end-1]

   f = Fs * (0:(L/2)) / L

   f, P₁
end