# Wave polarization analysis.
#
# Hongyang Zhou, hyzhou@umich.edu

using FFTW: fft
using LinearAlgebra: svd, Diagonal

const N = 10 # points-1 for DFT in time-frequency decomposition
#const c = 3e8 # speed of light
const c = 1.0 # normalized speed of light
const Fs = 1.0 # sampling rate, [Hz]

t = range(0, 10, length=N+1)
bx = @. sin(t)
by = @. sin(t + π/2)
#by = @. sin(t)
bz = zeros(length(t))

ex = @. -cos(t)
#ey = @. cos(t + π/2)
ey = @. cos(t)
#ex = zeros(length(t))
#ey = zeros(length(t))
ez = zeros(length(t))

f = Fs * (0:(length(t)/2)) / length(t)


## B-field SVD

b̃x = fft(bx)
b̃y = fft(by)
b̃z = fft(bz)

b̃ = (b̃x, b̃y, b̃z)
#=
L = length(x)
x̃ = fft(x)
# Compute the two-sided spectrum amplitude P2.
P₂ = @. abs(x̃/L)
# Compute the single-sided spectrum amplitude P1.
P₁ = P₂[1:floor(Int, L/2+1)]
P₁[2:end-1] .= 2 .* P₁[2:end-1]

=#


for it = eachindex(t)
   q = zeros(ComplexF64, 3, 3)
   for j in axes(q,2), i in axes(q,1)
      q[i,j] = b̃[i][it] * conj(b̃[j][it])
   end
   a11, a41 = reim(q[1,1])
   a21, a51 = reim(q[1,2])
   a31, a61 = reim(q[1,3])
   a12, a42 = reim(q[2,1])
   a22, a52 = reim(q[2,2])
   a32, a62 = reim(q[2,3])
   a13, a43 = reim(q[3,1])
   a23, a53 = reim(q[3,2])
   a33, a63 = reim(q[3,3])

   A = [
      a11 a12 a13;
      a21 a22 a23;
      a31 a32 a33;
      a41 a42 a43;
      a51 a52 a53;
      a61 a62 a63]

   # F.U * Diagonal(F.S) * F.Vt
   A_svd = svd(A)

   # I don't know about this
   #θk = atan(√(F.Vt[1]^2 ))

   # planarity of polarization
   F = 1 - √(A_svd.S[3] / A_svd.S[1])

   # ellipticity
   Lp = A_svd.S[2] / A_svd.S[1]

   @info "frequency comp $it"
   @info "Lp = $Lp"

end

#=
a11 = sum(@. real(b̃x * conj(b̃x)) )
a12 = sum(@. real(b̃y * conj(b̃x)) )
a13 = sum(@. real(b̃z * conj(b̃x)) )
a21 = sum(@. real(b̃x * conj(b̃y)) )
a22 = sum(@. real(b̃y * conj(b̃y)) )
a23 = sum(@. real(b̃z * conj(b̃y)) )
a31 = sum(@. real(b̃x * conj(b̃z)) )
a32 = sum(@. real(b̃y * conj(b̃z)) )
a33 = sum(@. real(b̃z * conj(b̃z)) )

a41 = sum(zeros(eltype(bx), N) )
a42 = sum(@. imag(b̃y * conj(b̃x)) )
a43 = sum(@. imag(b̃z * conj(b̃x)) )
a51 = sum(@. imag(b̃x * conj(b̃y)) )
a52 = sum(zeros(eltype(bx), N) )
a53 = sum(@. imag(b̃z * conj(b̃y)) )
a61 = sum(@. imag(b̃x * conj(b̃z)) )
a62 = sum(@. imag(b̃y * conj(b̃z)) )
a63 = sum(zeros(eltype(bx), N) )

A = [a11 a12 a13; a21 a22 a23; a31 a32 a33; a41 a42 a43; a51 a52 a53; a61 a62 a63]

# F.U * Diagonal(F.S) * F.Vt
A_svd = svd(A)

# I don't know about this
#θk = atan(√(F.Vt[1]^2 ))

# planarity of polarization
F = 1 - √(A_svd.S[3] / A_svd.S[1])

# ellipticity
Lp = A_svd.S[2] / A_svd.S[1]

@info "B-field analysis"
@info "Lp = $Lp"
=#
#=
## EM field SVD

b̃x = fft(bx)
b̃y = fft(by)
b̃z = fft(bz)

ẽx = fft(ex)
ẽy = fft(ey)
ẽz = fft(ez)

ζ = (c.*b̃x, c.*b̃y, c.*b̃z, ẽx, ẽy, ẽz)

# loop over each frequency
for it = eachindex(t)
   q = zeros(ComplexF64, 6, 6)
   for j in axes(q,2), i in axes(q,1)
      q[i,j] = ζ[i][it] * conj(ζ[j][it])
   end

   Aₑ = [
      0              real(q[6,1])  -real(q[5,1])
      0              real(q[6,2])  -real(q[5,2])
      0              real(q[6,3])  -real(q[5,3])
      0              real(q[6,4])  -real(q[5,4])
      0              real(q[6,5])  -real(q[5,5])
      0              real(q[6,6])  -real(q[5,6])
      -real(q[6,1])  0              real(q[4,1])
      -real(q[6,2])  0              real(q[4,2])
      -real(q[6,3])  0              real(q[4,3])
      -real(q[6,4])  0              real(q[4,4])
      -real(q[6,5])  0              real(q[4,5])
      -real(q[6,6])  0              real(q[4,6])
      real(q[5,1])   -real(q[4,1])  0
      real(q[5,2])   -real(q[4,2])  0
      real(q[5,3])   -real(q[4,3])  0
      real(q[5,4])   -real(q[4,4])  0
      real(q[5,5])   -real(q[4,5])  0
      real(q[5,6])   -real(q[4,6])  0
      0              imag(q[6,1])  -imag(q[5,1])
      0              imag(q[6,2])  -imag(q[5,2])
      0              imag(q[6,3])  -imag(q[5,3])
      0              imag(q[6,4])  -imag(q[5,4])
      0              imag(q[6,5])  -imag(q[5,5])
      0              imag(q[6,6])  -imag(q[5,6])
      -imag(q[6,1])  0              imag(q[4,1])
      -imag(q[6,2])  0              imag(q[4,2])
      -imag(q[6,3])  0              imag(q[4,3])
      -imag(q[6,4])  0              imag(q[4,4])
      -imag(q[6,5])  0              imag(q[4,5])
      -imag(q[6,6])  0              imag(q[4,6])
      imag(q[5,1])   -imag(q[4,1])  0
      imag(q[5,2])   -imag(q[4,2])  0
      imag(q[5,3])   -imag(q[4,3])  0
      imag(q[5,4])   -imag(q[4,4])  0
      imag(q[5,5])   -imag(q[4,5])  0
      imag(q[5,6])   -imag(q[4,6])  0
   ]
   
   rhs = [
      real(q[1,1]), real(q[1,2]), real(q[1,3]), real(q[1,4]), real(q[1,5]), real(q[1,6]),
      real(q[2,1]), real(q[2,2]), real(q[2,3]), real(q[2,4]), real(q[2,5]), real(q[2,6]),
      real(q[3,1]), real(q[3,2]), real(q[3,3]), real(q[3,4]), real(q[3,5]), real(q[3,6]),
      0,            imag(q[1,2]), imag(q[1,3]), imag(q[1,4]), imag(q[1,5]), imag(q[1,6]),
      imag(q[2,1]), 0           , imag(q[2,3]), imag(q[2,4]), imag(q[2,5]), imag(q[2,6]),
      imag(q[3,1]), imag(q[3,2]), 0           , imag(q[3,4]), imag(q[3,5]), imag(q[3,6]),
   ]

   # A_svd.U * Diagonal(A_svd.S) * A_svd.Vt
   A_svd = svd(Aₑ)

   n = transpose(A_svd.Vt) * Diagonal(inv.(A_svd.S)) * transpose(A_svd.U) * rhs

   # planarity of polarization
   F = 1 - √(A_svd.S[3] / A_svd.S[1])

   # ellipticity
   Lp = A_svd.S[2] / A_svd.S[1]

   @info "it = $it"
   @info "Lp = $Lp"
   @info "n = $n"
end
=#




