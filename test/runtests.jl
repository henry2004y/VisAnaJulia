# Test of SpaceAnalysis

ENV["MPLBACKEND"]="agg" # no GUI

using SpaceAnalysis, Test

@info "SpaceAnalysis is calling PyPlot for plotting."

@testset "SpaceAnalysis.jl" begin
   @testset "MVA" begin # minimum variance analysis
      include("MVA_galileo_G8.jl")
      F = MVA_analysis("data/Galileo_G8_flyby_MAG.dat.gz", 2656:2775)
      eigenRef = [2079.360,  78.142,   34.309]
      @test F.values ≈ eigenRef atol=1e-3
   end
   @testset "MFA" begin # mean field-aligned coordinate
      
   end
end
