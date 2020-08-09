# Test of VisAna

ENV["MPLBACKEND"]="agg" # no GUI

using VisAna, Batsrus, PyPlot, Test

@info("VisAna is calling PyPlot for plotting.")

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

@testset "2D field tracing" begin
   # Test for the bilinear interpolation function
   x, y = 0.1, 0.2
   Q00, Q01, Q10, Q11 = 3.0, 5.0, 40.0, 60.0
   sol1 = 7.460

   println("Testing bilin_reg")
   out = VisAna.bilin_reg(x, y, Q00, Q01, Q10, Q11)
   @test out ≈ sol1 atol=1e-5

   # Test cEuler 1
   nx, ny, maxstep = 841, 121, 10000

   xgrid = Vector{Float64}(undef,nx)   # grid x
   ygrid = Vector{Float64}(undef,nx)   # grid y
   ux = Vector{Float64}(undef,nx*ny)   # field x
   uy = Vector{Float64}(undef,nx*ny)   # field y
   xt = Vector{Float64}(undef,maxstep) # output x
   yt = Vector{Float64}(undef,maxstep) # output y
   ds = 1.0

   for i=1:nx
      xgrid[i] = -10.0+0.25*(i-1)
      ygrid[i] = xgrid[i]
   end

   for i=1:nx, j=1:ny
      ux[(i-1)*ny+j] = xgrid[i]
      uy[(i-1)*ny+j] = -1.0*ygrid[j]
   end

   @time npoints = VisAna.Euler!(nx, ny, maxstep, ds, 1.0, 10.0, xgrid, ygrid,
 			 ux, uy, xt, yt)

   #= This will cause error if libtrace.so not in path!
   @time npoints = ccall((:cEuler,"libtrace.so"),Cint,
             (Cint,Cint,Cint,Float64,Float64,Float64,Ptr{Cdouble},Ptr{Cdouble},
             Ptr{Cdouble},Ptr{Cdouble},Ptr{Cdouble},Ptr{Cdouble}),
             nx, ny, maxstep, ds, 1.0, 10.0, xgrid, ygrid, ux, uy, xt, yt)
   =#
   @test npoints == 800
   println("Npoints = ", npoints)
   println("Grid goes from ", round(xgrid[1],digits=2), " to ", round(xgrid[nx],digits=2))
   println("Our trace starts at ", round(xt[1],digits=2), " ", round(yt[1],digits=2))
   println("...and ends at ", round(xt[npoints],digits=2), " ",round(yt[npoints],digits=2))

   @time npoints = VisAna.RK4!(nx, ny, maxstep, ds, 1.0, 10.0, xgrid, ygrid,
 		 ux, uy, xt, yt)

   #= This will cause error if libtrace.so not in path!
   @time npoints = ccall((:cRk4,"libtrace.so"),Cint,
             (Cint,Cint,Cint,Float64,Float64,Float64,Ptr{Cdouble},Ptr{Cdouble},
             Ptr{Cdouble},Ptr{Cdouble},Ptr{Cdouble},Ptr{Cdouble}),
             nx, ny, maxstep, ds, 1.0, 10.0, xgrid, ygrid, ux, uy, xt, yt)
   =#
   @test npoints == 799
   println("Npoints = ", npoints)
   println("Grid goes from ", round(xgrid[1],digits=2), " to ", round(xgrid[nx],digits=2))
   println("Our trace starts at ", round(xt[1],digits=2), " ", round(yt[1],digits=2))
   println("...and ends at ", round(xt[npoints],digits=2), " ",round(yt[npoints],digits=2))

   @test VisAna.test_trace_asymptote()    # single precision
   @test VisAna.test_trace_asymptote(true)# double precision

   @test VisAna.test_dipole()            # dipole field plotting in 2D
   @test VisAna.test_trace_dipole()      # dipole tracing in 2D

   # field tracing using dipole+background uniform field in BATSRUS
   filename = "y=0_var_1_t00000000_n00000000.out"
   data = readdata(filename, dir="data")

   streamplot(data, "bx;bz")

   bx = data.w[:,:,5]
   bz = data.w[:,:,7]
   x  = data.x[:,1,1]
   z  = data.x[1,:,2]

   #xstart = -120.0:10.0:10.
   #zstart = fill(10.0, size(xstart))

   seeds = select_seeds(x,z)
   xs, zs = seeds[1,end], seeds[2,end]
   # Forward
   x1, z1 = trace2d_eul(bx, bz, xs, zs, x, z, ds=0.1, gridType="ndgrid", direction="forward")
   # Backward
   x2, z2 = trace2d_rk4(bx, bz, xs, zs, x, z, ds=0.1, gridType="ndgrid", direction="backward")
   @test length(x1) == 54
   @test length(x2) == 25
end

@testset "MVA" begin
   @info("Minimum variance analysis test.")
   include("../space/MVA.jl")
   F = MVA_analysis("data/Galileo_G8_flyby_MAG.dat.gz", 2656:2775)
   eigenRef = [2079.360,  78.142,   34.309]
   @test F.values ≈ eigenRef atol=1e-3
end
