# This file contains general functions used several places

# This is the largest possible 32-bit number
global EPS32 = 21474836480

# This is the largest possible 64-bit number
global EPS64 = 18446744073709551615

# Defines parameters used in both bequest and no bequest model
# FIXME: Do we really need all this parameters? Perhaps clean it up a bit
Base.@kwdef mutable struct Param
    T:: Float32 = 75 #Got a large error "AssertionError: isfinite(phi_d) && isfinite(gphi)..." when set to 100
    # Number of equidistant intervals in time
    N:: Int64 = 300
    Δt:: Float32 = T / N
    t:: Array{Float32} = collect(range(0, T, N + 1))

    # Spatial dimensions
    dim = 1

    # Common model parameters
    μ:: Float32 = 0.01 # Deterministic income growth rate
    σ:: Float32 = 0.1 # Stochastic income volatility (not used here)
    r:: Float32 = 0.02 # Return to savings = interest rate (currently fixed)
    ρ :: Float32 = 0.05 # Discount rate
    γ:: Float32 = 0.5 # Running utility exponent
    ψ:: Float32 = (r-ρ)/γ # Useful pre-calculated quantity


    # Bequest model extra parameters
    ε:: Float32 = 0.5 # Terminal utility exponent AVI: WE SHOULD CHANGE THIS NOTATION NOT TO CLASH WITH APPROXIMATION EPSILON
    λ:: Float32 = 15 # Terminal utility coefficient

    # Mock parameters
    β:: Float32 = 1000 # Terminal target
    κ:: Float32 = 5 # Running cost coefficient


    # Number of Monte-Carlo
    n:: Int64 = 50
    s:: Int64 = 1
    #x::Array{Float32, 0} = zeros()


    y0::Float32 = 1 # Initial income value
    y::Array{Float32} =y0 * exp.(μ*t)  # Pre-calculate income process

    a0::Float32 = 100 # Initial wealth
end

function crra(c, gamma)
    if c < 0
        return -Inf
    end
    if gamma == 1
        log(c)
    else
        c.^(1 - gamma) / (1 - gamma)
    end
end

function crra_talk(c, gamma)
    if c < 0
        10000*c
    else
        c^(1 - gamma) / (1 - gamma)
    end
end

function exp_util(c, a)
    (1 - exp(-a* c)) / a
end

function violation_penalty(x, eps)
    1000*eps.*x - eps
end

function smooth_crra(c, eps, gamma)
    if c < 0
        violation_penalty(c, eps)
    else
        crra(c, gamma)
    end
end

function mollified_crra(c, gamma, epsilon, p)
    if c >= epsilon
        c^(1 - gamma) / (1 - gamma)
    else
        -(c^2)/(2 * epsilon^p) + (1 / epsilon^gamma + 1 / epsilon^(p - 1)) * c - 1 / (2 * epsilon^(p -2)) + epsilon^(1 - gamma)*gamma / (1 - gamma)
    end
end