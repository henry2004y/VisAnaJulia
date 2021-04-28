# Test of VisAna

ENV["MPLBACKEND"]="agg" # no GUI

using VisAna, Batsrus, PyPlot, Test

@info "VisAna is calling PyPlot for plotting."

@testset "VisAna.jl" begin
   @testset "MVA" begin
      # Minimum variance analysis
      include("../space/MVA.jl")
      F = MVA_analysis("data/Galileo_G8_flyby_MAG.dat.gz", 2656:2775)
      eigenRef = [2079.360,  78.142,   34.309]
      @test F.values ≈ eigenRef atol=1e-3
   end
end
