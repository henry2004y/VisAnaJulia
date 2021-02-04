# Test of VisAna

ENV["MPLBACKEND"]="agg" # no GUI

using VisAna, Batsrus, PyPlot, Test

@info "VisAna is calling PyPlot for plotting."

@testset "VisAna.jl" begin
   @testset "plotting 1D ascii" begin
      filename = "1d__raw_2_t25.60000_n00000258.out"
      data = readdata(filename, dir="data", verbose=true)
      plotdata(data, "p", plotmode="line")
      line = gca().lines[1]
      @test line.get_xdata() ≈ data.x
      @test line.get_ydata() ≈ data.w[:,10]
   end

   @testset "plotting 2D structured binary" begin
      filename = "z=0_raw_1_t25.60000_n00000258.out"
      data = readdata(filename, dir="data")
      plotdata(data,"p bx;by",plotmode="contbar streamover")
      ax = gca()
      @test isa(ax, PyPlot.PyObject)
      contourf(data,"p")
      ax = gca()
      @test isa(ax, PyPlot.PyObject)
   end

   @testset "MVA" begin
      # Minimum variance analysis
      include("../space/MVA.jl")
      F = MVA_analysis("data/Galileo_G8_flyby_MAG.dat.gz", 2656:2775)
      eigenRef = [2079.360,  78.142,   34.309]
      @test F.values ≈ eigenRef atol=1e-3
   end
end
