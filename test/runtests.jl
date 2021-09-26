# Test of SpaceAnalysis

using SpaceAnalysis, Test, Random, DelimitedFiles

"Calculates the magnetic field from a dipole with magnetic moment `M` at `r`."
function dipole(rIn, M)
   x, y, z = rIn
   r = sqrt(x^2 + y^2 + z^2)
   Coef = 1e-7 / r^5 # μ₀/(4*π*r^5)

   B = [3*x^2-r^2 3*x*y     3*x*z;
        3*y*x     3*y^2-r^2 3*y*z;
        3*z*x     3*z*y     3*z^2-r^2
       ] * M * Coef
end

@testset "SpaceAnalysis.jl" begin
   @testset "Spectrum" begin
      Fs = 1000.     # Sampling frequency                    
      T = 1/Fs       # Sampling period       
      L = 1500       # Length of signal
      t = (0:L-1)*T  # Time vector
      # Form a signal containing a 50 Hz sinusoid of amplitude 0.7 and
      # a 120 Hz sinusoid of amplitude 1.
      S = @. 0.7*sinpi(2*50*t) + sinpi(2*120*t)
      # Corrupt the signal with zero-mean white noise with a variance of 4.
      X = S + 2*randn(MersenneTwister(1234), length(t))

      f, P1 = spectrum(X, Fs)
      @test f[end] == 500. && maximum(P1) == 0.9513531477953538

      @test mag2db(f)[end] == 26.989700043360187
   end

   @testset "MVA" begin # minimum variance analysis
      include("MVA_galileo_G8.jl")
      F = mva_analysis("data/Galileo_G8_flyby_MAG.dat.gz", 2656:2775)
      eigenRef = [2079.360,  78.142,   34.309]
      @test F.values ≈ eigenRef atol=1e-3
   end

   @testset "MFA" begin # mean field-aligned coordinate
      # Generate artificial signal
      f = 7e-3 # frequency, [Hz]
      h = 1    # sample period, [s]
      t = range(0.0, 286.0, step=h)
      y = @. sinpi(2 * f * t)
      ynoise = y + 0.1*randn(MersenneTwister(1234), length(t))

      fmax = 1e-2
      fsample = 1 / h
      ybk = emd_average(ynoise, fmax, fsample)

      @test ybk[100] ≈ -0.9506051118519016

      # Transfroming an ideal Earth dipole field along a straight line in x
      Rₑ = 6.378e6
      x = range(3, 20, step=0.1) * Rₑ
      y = ones(size(x))
      z = 2 .* y
      Bx, By, Bz = similar(z), similar(z), similar(z)
      BMoment_Earth = [0.0, 0.0, 7.94e22] # [V*s/(A*m)]
      for i in axes(x, 1)
         Bx[i], By[i], Bz[i] = dipole([x[i], y[i], z[i]], BMoment_Earth)
      end
      # Bμ ∥ B should be dominant.
      Bμ, Bϕ, Bν = mfa(x, y, z, Bx, By, Bz; method="MAVG", nw=15)

      @test Bμ[end] == 3.8254010321827695e-9

      # Introduce 2 oscillations with different freqs s.t. there will be multiple IMFs.
      Bx += 1e-9*sinpi.(2 * f * t[1:length(x)])
      Bx += 1e-10*cospi.(2 * 10f * t[1:length(x)])

      # Bμ ∥ B should be dominant.
      Bμ, Bϕ, Bν = mfa(x, y, z, Bx, By, Bz; method="EMD", fsample, fmax, verbose=true)

      @test Bμ[end] == 3.956605392869394e-9
   end

   @testset "signal" begin # signal generation
      # Alfvénic perturbation
      n0 = 1e6 # density, [amu/m^3]
      T0 = 1e6 # temperature, [K]
      V0 = [-700., 0.0, 0.0]*1e3 # velocity, [m/s] 
      B0 = [-3., 3., 0.]*1e-9 # magnetic field, [T]

      varBG = BackgroundVariable(n0, T0, V0..., B0...)

      s = generate_signal(varBG; fsignal=0.1, signal=:alfven, dv=1e5, fsample=2, tend=10,
         dir="y")

      save_signal(s)

      d = readdlm("sw.dat", '\t', Float64, header=true)

      @test d[1][10, end-1] == 1.4167e-9
      rm("sw.dat")

      # Density pulse
      n0 = 2e6 # density, [amu/m^3]
      T0 = 5e5 # temperature, [K]
      V0 = [-600., 0.0, 0.0]*1e3 # [m/s] 
      B0 = [5., 0., 0.]*1e-9 # [T]
      
      varBG = BackgroundVariable(n0, T0, V0..., B0...)
      
      fsignal = 1/150 # perturbation frequency, [Hz]
      signal  = :fast # perturbation type
      dv      = 5e4   # perturbation magnitude in velocity, [m/s]
      dn      = 1     # perturbation magnitude in density, [amu/cc]
      dir     = "y"   # perturbed direction(s)
      fsample = 2     # sampling frequency, [Hz]
      tstart  = 300   # [s]
      tend    = 450   # [s]
      
      s = generate_signal(varBG; fsignal, signal, dv, fsample, tstart, tend, dir)
      @test s.n[end] ≈ n0

   end
end
