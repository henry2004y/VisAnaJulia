# Mean Field Aligned (MFA) coordinates transformation.

using EmpiricalModeDecomposition
using EmpiricalModeDecomposition: zerocrossing!
using Statistics: mean

export mfa, emd_average

"""
    mfa(x, y, z, Bx, By, Bz; method="MAVG", nw=15, fmax=-1.0, fsample=-1.0, threshold=0.1,
       verbose=false) -> Bμ, Bϕ, Bν

Obtain the time series of rotation matrix from geocentric solar ecliptic (GSE) to Mean
field-aligned (MFA) reference frame, given the positions of the spacecraft `x`, `y`, `z` and
magnetic fields `Bx`, `By`, `Bz` in GSE. 

# Optional Arguments
- `method::String="MAVG"`: approach to obtain background ambient field, from moving average
`"MAVG"` with `nw` points, to empirical mode decomposition `"EMD"` with cutoff frequency
`fmax`, sampling frequency `fsample` and power selection `threshold` ∈ (0,1).
"""
function mfa(x::AbstractVector, y::AbstractVector, z::AbstractVector,
   Bx::AbstractVector, By::AbstractVector, Bz::AbstractVector;
   method="MAVG", nw=15, fsample=-1.0, fmax=fsample/50 , threshold=0.1, verbose=false)

   if !(length(x) == length(y) == length(z) ==
      length(Bx) == length(By) == length(Bz))
      throw(DimensionMismatch("requires same vector length for inputs!"))
   end

   μx, μy, μz = similar(x), similar(x), similar(x)
   ϕx, ϕy, ϕz = similar(x), similar(x), similar(x)
   νx, νy, νz = similar(x), similar(x), similar(x)
   Bμ, Bϕ, Bν = similar(x), similar(x), similar(x)

   if method == "MAVG"
     B̄x = sma(Bx, nw)
     B̄y = sma(By, nw)
     B̄z = sma(Bz, nw)
   elseif method == "EMD"
     B̄x = emd_average(Bx, fmax, fsample, threshold; verbose)
     B̄y = emd_average(By, fmax, fsample, threshold; verbose)
     B̄z = emd_average(Bz, fmax, fsample, threshold; verbose)
   end

   for i in axes(x, 1)
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

      Bμ[i] = μx[i]*Bx[i] + μy[i]*By[i] + μz[i]*Bz[i]
      Bϕ[i] = ϕx[i]*Bx[i] + ϕy[i]*By[i] + ϕz[i]*Bz[i]
      Bν[i] = νx[i]*Bx[i] + νy[i]*By[i] + νz[i]*Bz[i]
   end

   if verbose
      @info "method = $method"
   end

   Bμ, Bϕ, Bν
end


"""
    emd_average(y, fmax, fsample, threshold=0.1; verbose=false)

Extracts low frequency component from a signal using Empirical Mode Decomposition (EMD).
The signal is reconstructed by using Intrinsic Mode Functions (IMF) which have an average
frequency lower than a specified cut-off frequency `fmax` and a percentage of the power of
the signal, attributable to frequencies lower than `fmax` and greater than `thresold`.

Reference:
The use of the empirical mode decomposition for the identification of mean field aligned
reference frames, M Regi, Alfredo Del Corpo, Marcello De Lauretis, 2016
https://doi.org/10.4401/ag-7067
"""
function emd_average(y, fmax, fsample, threshold=0.1; verbose=false)

   len = length(y)
   x = range(1, length=len, step=1/fsample)

   imfs = emd(y, x)  # intrinsic mode functions

   if length(imfs) == 1
      verbose && @warn "No oscillation detected!"
      return imfs[1]
   end

   res = imfs[end]
   imf = imfs[1:end-1]

   ## Zero-crossing
   # maximum number of zero-crossing is the number of zero_crossing of a sin wave of
   # frequency fmax in a sampling time interval T0 = len/fsample
   nmax = floor(Int, length(y)*(2*fmax)/fsample+1)  

   # number of zeros in each IMF
   crossings = [Int[] for _ = 1:length(imfs)-1] # initialization

   for i in 1:length(crossings)
      zerocrossing!(imf[i], crossings[i])
   end

   nzeros = [length(crossings[i]) for i in 1:length(crossings)]

   ## Average period and signal power
   tmin = floor(Int, fsample/fmax)    # miminum period [# of samples]
   tminhalf = floor(Int, fsample/(2*fmax)) # miminum half period [# of samples]
   #T = zeros(length(crossings))
   ind1_ = Vector{Bool}(undef, length(crossings))

   for i = 1:length(imf)
      # indices of zero crossings plus first and last entries
      ind_ = [1, crossings[i]..., len]
      # lengths of consecutive zero-crossing [# of samples]
      tp = diff(ind_)
      #T[i] = mean(2*tp) # average period of the mode [# of samples]

      p = zeros(length(tp))
      # Evaluating the IMF power in each segment
      for j = 1:length(tp)
         p[j] = sum( abs.( imf[i][ind_[j]:ind_[j+1]-1] ) ) / fsample
      end

      ps = p[tp .> tminhalf] # select powers linked to periods greater than tmin 

      # select IMF if the percentage of the selected power is greater than the thresold
      ind1_[i] = sum(ps) / sum(p) > threshold

      verbose && @info "IMF $i, threshold = $threshold, power = $(sum(ps) / sum(p))"
   end
   
   ## Reconstructing ambient signal
   ind2_ = nzeros .< nmax
   imf_bk = imf[ind1_ .| ind2_] # number of IMFs to add to residue to obtain ambient signal

   y_bk = res + sum(imf_bk)
      
   return y_bk
end