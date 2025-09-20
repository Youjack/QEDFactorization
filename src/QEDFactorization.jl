module QEDFactorization

using cLHAPDF

using QuadGK: quadgk
using HCubature: hcubature
using Cuba: vegas, suave, divonne

export get_qed_data
export get_αEM
export get_qed_density, get_qed_density_integral

export qed_conv

#= QEDData ========================================================================================#

"""
    get_qed_data(setname_nmem::String)::Ptr{Cvoid}
"""
function get_qed_data(setname_nmem::String)::Ptr{Cvoid}
    mkPDF(setname_nmem)
end

"""
    get_αEM(qed_data::Ptr{Cvoid}, μ²)::Float64

Get αEM at scale `μ²`.
"""
function get_αEM(qed_data::Ptr{Cvoid}, μ²::Real)::Float64
    alphasQ2(qed_data, μ²)
end

"""
    get_qed_density(qed_data::Ptr{Cvoid}, l, μ²)::Function

Get QED DF or FF at scale `μ²`, as `f(x)`. \\
`x` is not constrained to be in `(0,1)`.
"""
function get_qed_density(qed_data::Ptr{Cvoid}, l::Integer, μ²::Real)::Function
    x -> 1/x * xfxQ2(qed_data, l, 1-x, μ²)
end

"""
    get_qed_density_integral(qed_data::Ptr{Cvoid}, l, μ²)::Function

Get integral (from `x` to `1`) of QED DF or FF at scale `μ²`. \\
`x` is not constrained to be in `(0,1)`.
"""
function get_qed_density_integral(qed_data::Ptr{Cvoid}, l::Integer, μ²::Real)::Function
    x -> 1 - xfxQ2(qed_data, 100l, 1-x, μ²)
end

#= qed_conv =======================================================================================#

const INT_CUBA = 1
const INT_MC   = 2

"""
    qed_conv(f::Function, If::Function, D::Function, ID::Function, H::Function,
        Θ::Function, ξmin, ζmin; algor=INT_HCUBA, rtol=_rtol, seed=0)::Float64

Calculate `∬ dξ dζ f(ξ) D(ζ) Θ(ξ,ζ) H(ξ,ζ)` using the subtraction trick, where
`Θ` is a boolean-valued step function that constrains `(ξ,ζ)` to be in physical region.
- Cross-section usually vanishes at the lower bound of `(ξ,ζ)` and `INT_CUBA=1` could be used.
- For calculations of probability densities of kinematic variables, `INT_MC=2` could be used.
"""
function qed_conv(f::Function, If::Function, D::Function, ID::Function, H::Function,
        Θ::Function, ξmin, ζmin; algor=INT_CUBA, rtol=_rtol, seed=0)::Float64
    integrand2D(ξ,ζ) = isone(ξ) || isone(ζ) ? 0 :
        f(ξ) * D(ζ) * (
            + (Θ(ξ,ζ) ? H(ξ,ζ) : 0)
            - (Θ(ξ,1) ? H(ξ,1) : 0) - (Θ(1,ζ) ? H(1,ζ) : 0)
            + (Θ(1,1) ? H(1,1) : 0) )
    integral2D =
        if algor == INT_CUBA
            hcubature(X -> integrand2D(X[1],X[2]), (ξmin,ζmin),(1.,1.), rtol=rtol)[1]
        elseif algor == INT_MC
            vegas((X,F) -> F[1] =
                (1-ξmin)*(1-ζmin) * integrand2D(ξmin+(1-ξmin)*X[1],ζmin+(1-ζmin)*X[2]),
                2, rtol=rtol, seed=seed)[1][1]
        else
            throw(ArgumentError("qed_conv: algorithm $algor not supported."))
        end
    return integral2D +
        quadgk(ξ -> f(ξ) *( (Θ(ξ,1) ? H(ξ,1) : 0) - (Θ(1,1) ? H(1,1) : 0) ), ξmin,1., rtol=rtol)[1] * ID(ζmin) +
        quadgk(ζ -> D(ζ) *( (Θ(1,ζ) ? H(1,ζ) : 0) - (Θ(1,1) ? H(1,1) : 0) ), ζmin,1., rtol=rtol)[1] * If(ξmin) +
        (Θ(1,1) ? H(1,1) : 0) * If(ξmin) * ID(ζmin)
end

#= QEDEvolve ======================================================================================#

include("QEDEvolve.jl")

end # module
