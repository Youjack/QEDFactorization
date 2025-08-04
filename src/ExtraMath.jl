#= Part of QEDEvolve.jl                                                                           =#
#= Special functions copied from `SpecialFunctions.jl`                                            =#

const ComplexOrReal{T} = Union{T,Complex{T}}

digamma(x::Number) = _digamma(float(x))
function _digamma(z::ComplexOrReal{Float64})
    # Based on eq. (12), without looking at the accompanying source
    # code, of: K. S. Kölbig, "Programs for computing the logarithm of
    # the gamma function, and the digamma function, for complex
    # argument," Computer Phys. Commun.  vol. 4, pp. 221–226 (1972).
    x = real(z)
    if x <= 0 # reflection formula
        ### ψ = -π * _cotpi(z)
        ψ = -π * cot(π * z)
        z = 1 - z
        x = real(z)
    else
        ψ = zero(z)
    end
    X = 8
    if x < X
        # shift using recurrence formula
        n = X - floor(Int,x)
        for ν = 1:n-1
            ψ -= inv(z + ν)
        end
        ψ -= inv(z)
        z += n
    end
    t = inv(z)
    ψ += log(z) - 0.5*t
    t *= t # 1/z^2
    # the coefficients here are Float64(bernoulli[2:9] .// (2*(1:8)))
    ### ψ -= t * @evalpoly(t,0.08333333333333333,-0.008333333333333333,0.003968253968253968,-0.004166666666666667,0.007575757575757576,-0.021092796092796094,0.08333333333333333,-0.4432598039215686)
    ψ -= (
        t^1 *  0.08333333333333333  +
        t^2 * -0.008333333333333333 +
        t^3 *  0.003968253968253968 +
        t^4 * -0.004166666666666667 +
        t^5 *  0.007575757575757576 +
        t^6 * -0.021092796092796094 +
        t^7 *  0.08333333333333333  +
        t^8 * -0.4432598039215686    )
end

trigamma(x::Number) = _trigamma(float(x))
function _trigamma(z::ComplexOrReal{Float64})
    # via the derivative of the Kölbig digamma formulation
    x = real(z)
    if x <= 0 # reflection formula
        ### return (π / sinpi(z))^2 - trigamma(1 - z)
        return (π / sin(π * z))^2 - trigamma(1 - z)
    end
    ψ = zero(z)
    N = 10
    if x < N
        # shift using recurrence formula
        n = N - floor(Int,x)
        ψ += inv(z)^2
        for ν = 1:n-1
            ψ += inv(z + ν)^2
        end
        z += n
    end
    t = inv(z)
    w = t * t # 1/z^2
    ψ += t + 0.5*w
    # the coefficients here are Float64(bernoulli[2:9])
    ### ψ += t*w * @evalpoly(w,0.16666666666666666,-0.03333333333333333,0.023809523809523808,-0.03333333333333333,0.07575757575757576,-0.2531135531135531,1.1666666666666667,-7.092156862745098)
    ψ += t * (
        w^1 *  0.16666666666666666  +
        w^2 * -0.03333333333333333  +
        w^3 *  0.023809523809523808 +
        w^4 * -0.03333333333333333  +
        w^5 *  0.07575757575757576  +
        w^6 * -0.2531135531135531   +
        w^7 *  1.1666666666666667   +
        w^8 * -7.092156862745098     )
end
