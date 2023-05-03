ϵ₀ = 8.85e-12
mₑ = 9.109e-31
e = 1.6e-19
n = 1e9  # [m^-3]
B = 1e-7 # [T]
ωₑ = √(n/(mₑ*ϵ₀))*e
Ωₑ = e*B/mₑ
ν = ωₑ / (32π) * 10 / exp(10)

σ₀ = n*e^2/(mₑ*ν)
σₚ = ν^2 / (ν^2 + Ωₑ^2) * σ₀
σₕ = ν^2 / (ν^2 + Ωₑ^2) * σ₀