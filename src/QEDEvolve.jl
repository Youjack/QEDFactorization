#= Part of QEDFactorization.jl                                                                    =#
module  QEDEvolve

using SpecialFunctions: gamma, beta
using QuadGK
using StaticArrays

export evolve_αEM_from_μ₀²_to_μ²

export P̂ll, P̂s, ℙ̂s
export ΔP̂ll, ΔP̂s, Δℙ̂s

export fll_toy, f̂ll₀_toy
export fll_nlo, f̂ll₀
export fγl_nlo, f̂γl₀
export gll_nlo, ĝll₀
export gγl_nlo, ĝγl₀
export Dll_nlo, D̂ll₀
export Dlγ_nlo, D̂lγ₀

export evolve_nonsinglet_from_μ₀²_to_μ², evolve_singlet_from_μ₀²_to_μ²
export mellin_inv_integrand, mellin_inv
export integrate_distribution

#= ExtraMath ======================================================================================#

include("ExtraMath.jl")

const γE = -digamma(1)

#= αEM ============================================================================================#

# `d/dlnμ² αEM = - β₀ αEM²`
β₀(nf) = - nf / 3π
"""
    evolve_αEM_from_μ₀²_to_μ²(αEM, μ₀², μ², nf)

Evolve αEM from `μ₀²` to `μ²`. \\
Threshold crossing is not considered.
"""
function evolve_αEM_from_μ₀²_to_μ²(αEM::Real, μ₀²::Real, μ²::Real, nf::Int)::Float64
    inv( inv(αEM) + β₀(nf) * log( μ² / μ₀² ) )
end

#= Splitting kernels ==============================================================================#

# unpolarized
P̂ll(N) = 3/2 - 1/N - 1/(1+N) - 2( digamma(N) - digamma(1) )
P̂γl(N) = (2 + N + N^2)/((N - 1) * N * (N + 1))
P̂lγ(N) = (2 + N + N^2)/(N * (N + 1) * (N + 2))
P̂γγ(N) = - 2/3
P̂s(N) = @SMatrix [ P̂ll(N) 2P̂lγ(N); P̂γl(N) P̂γγ(N) ] # space-like singlet
ℙ̂s(N) = @SMatrix [ P̂ll(N) 2P̂γl(N); P̂lγ(N) P̂γγ(N) ] # time-like  singlet

# helicity
ΔP̂ll(N) = P̂ll(N)
ΔP̂γl(N) = (2 + N)/(N * (N + 1))
ΔP̂lγ(N) = (-1 + N)/(N * (N + 1))
ΔP̂γγ(N) = P̂γγ(N)
ΔP̂s(N) = @SMatrix [ ΔP̂ll(N) 2ΔP̂lγ(N); ΔP̂γl(N) ΔP̂γγ(N) ] # space-like singlet
Δℙ̂s(N) = @SMatrix [ ΔP̂ll(N) 2ΔP̂γl(N); ΔP̂lγ(N) ΔP̂γγ(N) ] # time-like  singlet

#= QED parton distributions =======================================================================#

"Fixed order toy `l→l` distribution."
fll_toy(a, b, x) =
    x^a * (1-x)^b / beta(1+a, 1+b)
"Mellin moment of toy `l→l` distribution at initial scale."
f̂ll₀_toy(a, b, N) =
    beta(a+N, 1+b) / beta(1+a, 1+b)

"Fixed order NLO `l→l` distribution."
fll_nlo(αEM, x, m², μ²) = αEM/2π *
    (1+x^2)/(1-x) * ( log(μ²/m²) - 2log(1-x) - 1 )
"Mellin moment of NLO `l→l` distribution at `μ∼m`."
f̂ll₀(αEM, N) = 1 + αEM/2π * (
    -P̂ll(N) +
    #-- from Mathematica --
    -1/6*(12 + N*(1 + N)*(36 + 12*γE^2*N*(1 + N) + (-21 + 2*pi^2)*N*(1 + N) + 
    12*γE*(1 + 2*N)) - 12*N*(1 + N)*(-(digamma(N)*(1 + 2*N*(1 + γE + γE*N) + 
    N*(1 + N)*digamma(N))) + N*(1 + N)*trigamma(N)))/(N^2*(1 + N)^2)
    #----------------------
)

"Fixed order NLO `l→γ` distribution."
fγl_nlo(αEM, x, m², μ²) = αEM/2π *
    (1+(1-x)^2)/x * ( log(μ²/m²) - 2log(x) - 1 )
"Mellin moment of NLO `l→γ` distribution at `μ∼m`."
f̂γl₀(αEM, N) = αEM/2π * (
    -P̂γl(N) +
    #-- from Mathematica --
    4/(-1 + N)^2 - 4/N^2 + 2/(1 + N)^2
    #----------------------
)

"Fixed order NLO `l→l` helicity distribution."
gll_nlo(αEM, x, m², μ²) = fll_nlo(αEM, x, m², μ²)
"Mellin moment of NLO `l→l` helicity distribution at `μ∼m`."
ĝll₀(αEM, N) = f̂ll₀(αEM, N)

"Fixed order NLO `l→γ` helicity distribution."
gγl_nlo(αEM, x, m², μ²) = αEM/2π *
    (2-x) * ( log(μ²/m²) - 2log(x) - (1-x)/(2-x) )
"Mellin moment of NLO `l→γ` helicity distribution at `μ∼m`."
ĝγl₀(αEM, N) = αEM/2π * (
    #-- from Mathematica --
    (4 + N * (7 + N))/(N^2 * (1 + N)^2)
    #----------------------
)

"Fixed order NLO `l→l` fragmentation function."
Dll_nlo(αEM, x, m², μ²) = fll_nlo(αEM, x, m², μ²)
"Mellin moment of NLO `l→l` fragmentation function at `μ∼m`."
D̂ll₀(αEM, N) = f̂ll₀(αEM, N)

"Fixed order NLO `γ→l` fragmentation function."
Dlγ_nlo(αEM, x, m², μ²) = αEM/2π *
    (x^2+(1-x)^2) * log(μ²/m²)
"Mellin moment of NLO `γ→l` fragmentation function at `μ∼m`."
D̂lγ₀(αEM, N) = 0.0

#= DGLAP evolutions ===============================================================================#

"""
    evolve_nonsinglet_from_μ₀²_to_μ²(P̂::Function, αEM, f̂::Function, μ₀², μ², nf)::Function

Evolve `l→(l-l̄)` distribution from `μ₀²` to `μ²`. \\
Threshold crossing is not considered.
"""
function evolve_nonsinglet_from_μ₀²_to_μ²(P̂::Function, αEM::Real, f̂::Function,
        μ₀²::Real, μ²::Real, nf::Int)::Function
    N -> let
        if iszero(nf) return f̂(N) end
        αEM_ratio = αEM / evolve_αEM_from_μ₀²_to_μ²(αEM, μ₀², μ², nf)
        return αEM_ratio^( P̂(N) /( 2π * β₀(nf) ) ) * f̂(N)
    end
end
"""
    evolve_singlet_from_μ₀²_to_μ²(P̂::Function, αEM, f̂::Function, μ₀², μ², nf)::Function

Evolve `l→(l+l̄)` and `l→γ` distributions from `μ₀²` to `μ²`. \\
`P̂` should give a 2×2 matrix and `f̂` should give a 2D vector. \\
Threshold crossing is not considered.
"""
function evolve_singlet_from_μ₀²_to_μ²(_P̂::Function, αEM::Real, f̂::Function,
        μ₀²::Real, μ²::Real, nf::Int)::Function
    N -> let
        if iszero(nf) return f̂(N) end
        αEM_ratio = αEM / evolve_αEM_from_μ₀²_to_μ²(αEM, μ₀², μ², nf)
        I = [ 1 0; 0 1 ]
        P̂ = _P̂(N)
        λ₊ = 1/2 * ( P̂[1,1] + P̂[2,2] + sqrt( (P̂[1,1] - P̂[2,2])^2 + 4 * P̂[1,2] * P̂[2,1] ) )
        λ₋ = 1/2 * ( P̂[1,1] + P̂[2,2] - sqrt( (P̂[1,1] - P̂[2,2])^2 + 4 * P̂[1,2] * P̂[2,1] ) )
        U = (P̂ - λ₋ * I)/(λ₊ - λ₋) * αEM_ratio^( λ₊ /( 2π * β₀(nf) ) ) +
            (P̂ - λ₊ * I)/(λ₋ - λ₊) * αEM_ratio^( λ₋ /( 2π * β₀(nf) ) )
        return U * f̂(N)
    end
end

#= Mellin transformation ==========================================================================#

const _rtol = 1e-4

"""
    mellin_inv_integrand(f̂::Function, x)::Function

The integrand in Mellin inversion formula of `f̂`.
"""
function mellin_inv_integrand(f̂::Function, x::Real)::Function
    r -> let c = 1.9, ϕ = 3π/4
        N = c + r * exp(im * ϕ)
        1/π * inv(x)^( c + r * cos(ϕ) ) * imag(
            exp( im * ( ϕ - r * sin(ϕ) * log(x) ) ) * f̂(N) )
    end
end
"""
    mellin_inv(f̂::Function)::Function

Mellin inversion of `f̂`.
"""
function mellin_inv(f̂::Function)::Function
    x -> let
        I, E = quadgk(mellin_inv_integrand(f̂, x), 0,Inf)
        if E > _rtol
            @warn "mellin_inv: relative err = $E"
        end
        return I
    end
end

"""
    integrate_distribution(f̂::Function, xmin)::Float64

Integrate distribution `f` from `xmin` to `1`, using its Mellin moment `f̂`.
"""
function integrate_distribution(f̂::Function, xmin::Real)::Float64
    xmin * mellin_inv( N -> f̂(N) /(N - 1) )(xmin)
end

end # module
