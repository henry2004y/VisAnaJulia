# Test of SpaceAnalysis

using SpaceAnalysis, Test

@testset "SpaceAnalysis.jl" begin
   @testset "Spectrum" begin
      Fs = 1000.     # Sampling frequency                    
      T = 1/Fs       # Sampling period       
      L = 1500       # Length of signal
      t = (0:L-1)*T  # Time vector
      # Form a signal containing a 50 Hz sinusoid of amplitude 0.7 and
      # a 120 Hz sinusoid of amplitude 1.
      S = @. 0.7*sin(2π*50*t) + sin(2π*120*t)
      # Corrupt the signal with zero-mean white noise with a variance of 4.
      X = S + 2*randn(length(t))

      f, P1 = spectrum(X, Fs)
      @test f[end] == 500. && maximum(P1) == 1.0303307012340501

      @test mag2db(f)[end] == -12.760865229408134
   end

   @testset "MVA" begin # minimum variance analysis
      include("MVA_galileo_G8.jl")
      F = MVA_analysis("data/Galileo_G8_flyby_MAG.dat.gz", 2656:2775)
      eigenRef = [2079.360,  78.142,   34.309]
      @test F.values ≈ eigenRef atol=1e-3
   end
   @testset "MFA" begin # mean field-aligned coordinate
      
   end
end
