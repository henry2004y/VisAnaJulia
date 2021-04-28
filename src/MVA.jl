# Minimum Variance Analysis (MVA)

using Statistics, LinearAlgebra

export MVA

"""
    MVA(Bx, By, Bz; verbose=true)

If the λ1 > λ2 > λ3 are 3 eigenvalues of the matrix M, then a good indicator
of nice fitting LMN coordinate system should have λ2/λ3 > 5.
"""
function MVA(Bx, By, Bz; verbose=false)

   B̄1 = mean(Bx)
   B̄2 = mean(By)
   B̄3 = mean(Bz)
   B̄11= mean(Bx.*Bx) - B̄1*B̄1
   B̄22= mean(By.*By) - B̄2*B̄2
   B̄33= mean(Bz.*Bz) - B̄3*B̄3
   B̄12= mean(Bx.*By) - B̄1*B̄2
   B̄23= mean(By.*Bz) - B̄2*B̄3
   B̄31= mean(Bz.*Bx) - B̄3*B̄1
   # Construct the matrix
   M = [B̄11 B̄12 B̄31; B̄12 B̄22 B̄23; B̄31 B̄23 B̄33]

   # Compute the eigen values and ratios (descending order)
   F = eigen(M, sortby = x -> -abs(x))

   if verbose
      println("Eigenvalues:",F.values)
      println("Eigenvectors:")
      println(F.vectors[:,1])
      println(F.vectors[:,2])
      println(F.vectors[:,3])
      println("Ratio of intermediate variance to minimum variance = ",
         round(F.values[2]/F.values[3],digits=3))
   end
   F
end


