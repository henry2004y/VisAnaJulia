# Test of SpaceAnalysis

using SpaceAnalysis, Test, Random, DelimitedFiles

rng = MersenneTwister(1234)

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
      X = S + 2*randn(rng, length(t))

      f, P1 = spectrum(X, Fs)
      @test f[end] == 500. && maximum(P1) == 0.9513531477953495

      @test mag2db(f)[end] == 26.989700043360187
   end

   @testset "MVA" begin # minimum variance analysis
      include("MVA_galileo_G8.jl")
      F = MVA_analysis("data/Galileo_G8_flyby_MAG.dat.gz", 2656:2775)
      eigenRef = [2079.360,  78.142,   34.309]
      @test F.values ≈ eigenRef atol=1e-3
   end

   @testset "MFA" begin # mean field-aligned coordinate
      
   end

   @testset "signal" begin # signal generation
      n0 = 1e6 # density, [amu/m^3]
      T0 = 1e6 # temperature, [K]
      V0 = [-700., 0.0, 0.0]*1e3 # velocity, [m/s] 
      B0 = [-3., 3., 0.]*1e-9 # magnetic field, [T]

      varBG = VariableBackground(n0, T0, V0..., B0...)

      s = generate_signal(varBG; fsignal=1, dv=1e5, fsample=2, tend=10, dir="y")

      save_signal(s)

      d = readdlm("sw.dat", '\t', Float64, header=true)

      @test d[1][10, end-1] == -4.569e-9
      rm("sw.dat")
   end
end
