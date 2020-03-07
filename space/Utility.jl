"Import satellite data. Require one line header presented."
function read_data(fname)
   f = readdlm(fname, ',', Float32, '\n'; header=true)

   header = f[2][1:end-1]
   data   = f[1]

   satelliteNo = unique(data[:,1]) # number of static satellites

   return header, data, satelliteNo
end


"""
	smooth(x, n=100)

Return the moving box average of the vector data `x` with box length 'n'.
One-sided average on the left and right edge.
"""
function smooth(x::Vector, n=100)
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
	smooth(x, n)

Return the moving box average of the array data `x` with box length 'n'.
One-sided average on the left and right edge.
"""
function smooth(x::Array{T,3}, n=100) where T <: AbstractFloat
   nx = size(x)[3]
   if nx < n
      println("fewer snapshots than nSmooth, changing nS to nx-1")
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
