#= Part of QEDFactorization.jl                                                                    =#
module  QEDEvolve

using SpecialFunctions: gamma, beta
using QuadGK

export evolve_αEM_from_μ₀²_to_μ²
export fll_toy, f̂ll₀_toy
export fll_nlo, f̂ll₀
export fγl_nlo, f̂γl₀
export Dll_nlo, D̂ll₀
export Dγl_nlo, D̂γl₀
export evolve_nonsinglet_f̂_from_μ₀²_to_μ², evolve_singlet_f̂_from_μ₀²_to_μ²
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

"Mellin moment of `l→l` splitting kernel."
P̂ll(s) = 3/2 - 1/s - 1/(1+s) - 2( digamma(s) - digamma(1) )
"Mellin moment of `l→γ` splitting kernel."
P̂γl(s) = (2 + s + s^2)/((s - 1) * s * (s + 1))
"Mellin moment of `γ→l` splitting kernel."
P̂lγ(s) = (2 + s + s^2)/(s * (s + 1) * (s + 2))
"Mellin moment of `γ→γ` splitting kernel."
P̂γγ(s) = - 2/3

#= QED parton distributions =======================================================================#

"Fixed order toy `l→l` distribution."
fll_toy(a, b, x) =
    x^a * (1-x)^b / beta(1+a, 1+b)
"Mellin moment of toy `l→l` distribution at initial scale."
f̂ll₀_toy(a, b, s) =
    beta(a+s, 1+b) / beta(1+a, 1+b)

"Fixed order NLO `l→l` distribution."
fll_nlo(αEM, x, μ₀², μ²) = αEM/2π *
    (1+x^2)/(1-x) * ( log(μ²/μ₀²) - 2log(1-x) - 1 )
"Mellin moment of NLO `l→l` distribution at `μ∼μ₀`."
f̂ll₀(αEM, s) = 1 + αEM/2π * (
    -P̂ll(s) +
    #-- from Mathematica --
    -1/6*(12 + s*(1 + s)*(36 + 12*γE^2*s*(1 + s) + (-21 + 2*pi^2)*s*(1 + s) + 
    12*γE*(1 + 2*s)) - 12*s*(1 + s)*(-(digamma(s)*(1 + 2*s*(1 + γE + γE*s) + 
    s*(1 + s)*digamma(s))) + s*(1 + s)*trigamma(s)))/(s^2*(1 + s)^2)
    #----------------------
)

"Fixed order NLO `l→γ` distribution."
fγl_nlo(αEM, x, μ₀², μ²) = αEM/2π *
    (1+(1-x)^2)/x * ( log(μ²/μ₀²) - 2log(x) - 1 )
"Mellin moment of NLO `l→γ` distribution at `μ∼μ₀`."
f̂γl₀(αEM, s) = αEM/2π * (
    -P̂γl(s) +
    #-- from Mathematica --
    4/(-1 + s)^2 - 4/s^2 + 2/(1 + s)^2
    #----------------------
)

#= QED fragmentation functions ====================================================================#

"Fixed order NLO `l→l` fragmentation function."
Dll_nlo(αEM, z, μ₀², μ²) = αEM/2π *
    (1+z^2)/(1-z) * ( log(μ²/μ₀²) + 2log(z) - 2log(1-z) - 1 )
"Mellin moment of NLO `l→l` fragmentation function at `μ∼μ₀`."
D̂ll₀(αEM, s) = f̂ll₀(αEM, s) + αEM/2π * (
    #-- from Mathematica --
    2*(s^(-2) + (1 + s)^(-2) - 2*trigamma(s))
    #----------------------
)

"Fixed order NLO `l→γ` fragmentation function."
Dγl_nlo(αEM, z, μ₀², μ²) = αEM/2π *
    (1+(1-z)^2)/z * ( log(μ²/μ₀²) - 1 )
"Mellin moment of NLO `l→γ` fragmentation function at `μ∼μ₀`."
D̂γl₀(αEM, s) = αEM/2π *
    -P̂γl(s)

# DGLAP evolutions ---------------------------------------------------------------------------------

"""
    evolve_nonsinglet_f̂_from_μ₀²_to_μ²(αEM, f̂::Function, μ₀², μ², nf)::Function

Evolve `l→(l-l̄)` distribution from `μ₀²` to `μ²`. \\
Threshold crossing is not considered.
"""
function evolve_nonsinglet_f̂_from_μ₀²_to_μ²(αEM::Real, f̂::Function,
        μ₀²::Real, μ²::Real, nf::Int)::Function
    s -> let
        if iszero(nf) return f̂(s) end
        αEM_ratio = αEM / evolve_αEM_from_μ₀²_to_μ²(αEM, μ₀², μ², nf)
        return αEM_ratio^( P̂ll(s) /( 2π * β₀(nf) ) ) * f̂(s)
    end
end
"""
    evolve_singlet_f̂_from_μ₀²_to_μ²(αEM, f̂::Function, μ₀², μ², nf)::Function

Evolve `l→(l+l̄)` and `l→γ` distributions from `μ₀²` to `μ²`. \\
`f̂` should be `Vector`-valued. \\
Threshold crossing is not considered.
"""
function evolve_singlet_f̂_from_μ₀²_to_μ²(αEM::Real, f̂::Function,
        μ₀²::Real, μ²::Real, nf::Int)::Function
    s -> let
        if iszero(nf) return f̂(s) end
        αEM_ratio = αEM / evolve_αEM_from_μ₀²_to_μ²(αEM, μ₀², μ², nf)
        I = [ 1 0; 0 1 ]
        P = [ P̂ll(s) 2P̂lγ(s); P̂γl(s) P̂γγ(s) ]
        λ₊ = 1/2 * ( P[1,1] + P[2,2] + sqrt( (P[1,1] - P[2,2])^2 + 4 * P[1,2] * P[2,1] ) )
        λ₋ = 1/2 * ( P[1,1] + P[2,2] - sqrt( (P[1,1] - P[2,2])^2 + 4 * P[1,2] * P[2,1] ) )
        U = (P - λ₋ * I)/(λ₊ - λ₋) * αEM_ratio^( λ₊ /( 2π * β₀(nf) ) ) +
            (P - λ₊ * I)/(λ₋ - λ₊) * αEM_ratio^( λ₋ /( 2π * β₀(nf) ) )
        return U * f̂(s)
    end
end

"""
    mellin_inv_integrand(f̂::Function, x)::Function

The integrand in Mellin inversion formula of `f̂`.
"""
function mellin_inv_integrand(f̂::Function, x::Real)::Function
    c = 1.9; ϕ = 3π/4
    r -> let
        s = c + r * exp(im * ϕ)
        1/π * inv(x)^( c + r * cos(ϕ) ) * imag(
            exp( im * ( ϕ - r * sin(ϕ) * log(x) ) ) * f̂(s) )
    end
end
"""
    mellin_inv(f̂::Function)::Function

Mellin inversion of `f̂`.
"""
function mellin_inv(f̂::Function)::Function
    x -> quadgk(mellin_inv_integrand(f̂, x), 0,Inf)[1]
end

"""
    integrate_distribution(f̂::Function, xmin)::Float64

Integrate distribution `f` from `xmin` to `1`, using its Mellin moment `f̂`.
"""
function integrate_distribution(f̂::Function, xmin::Real)::Float64
    xmin * mellin_inv( s -> f̂(s) /(s - 1) )(xmin)
end

end # module
