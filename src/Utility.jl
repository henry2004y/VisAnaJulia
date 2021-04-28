"Import satellite data. Require one line header presented."
function read_satellite_data(fname)
   f = readdlm(fname, ',', Float32, '\n'; header=true)

   header = f[2][1:end-1]
   data   = f[1]

   satelliteNo = unique(data[:,1]) # number of static satellites

   if minimum(satelliteNo) > 0.1
      @error "Some thing is wrong with the input data!"
   end

   return header, data, satelliteNo
end

"Convert magnitude to decibels."
function mag2db(y)
   ydb = 10. * log10.(y)
end


"""
    sma(x, n=100)

Return the moving box average of the vector data `x` with box length 'n'.
One-sided average on the left and right edge.
"""
function sma(x::Vector, n=100)
   nx = length(x)
   x̄ = zeros(eltype(x),length(x))

   # left points
   for i = 1:n
      x̄[i] = mean(x[1:(i+n)])
   end

   # middle points
   for i = (n+1):(nx-n)
      x̄[i] = mean(x[(i-n):(i+n)])
   end

   # right points
   for i = (nx-n+1):nx
      x̄[i] = mean(x[(i-n):end])
   end

   return x̄
end


"""
    sma(x, n)

Return the moving box average of the 3D array data `x` with box length 'n'.
One-sided average on the left and right edge.
This is so confusing:
1. what about arbitrary dimensions?
2. how to choose which dimension to do smoothing?
Probably going to change in the next version, with meta-programming!
"""
function sma(x::Array{T,3}, n=100) where T <: AbstractFloat
   nx = size(x)[3]
   if nx < n
      println("fewer snapshots than nsma, changing nS to nx-1")
      n = nx-1
   end
   x̄ = zeros(eltype(x),size(x))

   # left points
   for i = 1:n
      x̄[:,:,i] = mean(x[:,:,1:(i+n)]; dims=3)
   end

   # middle points
   for i = n+1:nx-n
      x̄[:,:,i] = mean(x[:,:,(i-n):(i+n)]; dims=3)
   end

   # right points
   for i = nx-n+1:nx
      x̄[:,:,i] = mean(x[:,:,(i-n):end]; dims=3)
   end

   return x̄
end


"""
    ema(x, n=5)

Return the exponentially moving box average of the vector data `x` with box
length 'n'.
Check https://github.com/JuliaQuant/MarketTechnicals.jl/blob/master/src/movingaverages.jl
"""
function ema(x::Vector, n=5)
   @warn "To be implemented!"
end
