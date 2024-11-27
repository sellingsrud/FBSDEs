# This file contains general functions used several places
using Parameters

# This is the largest possible 32-bit number
global EPS32 = 21474836480

# This is the largest possible 64-bit number
global EPS64 = 18446744073709551615

# Defines parameters used in both bequest and no bequest model
# FIXME: Do we really need all this parameters? Perhaps clean it up a bit
# TODO: Check whether parameters are admissible, possible to do this in the constructor
@with_kw mutable struct Param
    t0::Float32 = 0.0
    T:: Float32 = 75 #Got a large error "AssertionError: isfinite(phi_d) && isfinite(gphi)..." when set to 100 Effective parameters are different when we change the time horizon.
    N:: Int32 = 300 # Number of equidistant intervals in time
    Δt:: Float32 = T / N #SHOULD REVERSE THESE: FIX DELTA t TOLERANCE FIRST
    t:: Array{Float32} = collect(range(0, T, N + 1))

    # Spatial dimensions
    dim = 1 #This is the dimension of the output of the control variable

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

    # Mock parameters for toy OU problem
    β:: Float32 = 1000 # Terminal target
    κ:: Float32 = 5 # Running cost coefficient


    # Number of Monte-Carlo
    n:: Int32 = 50
    s:: Int32 = 1
    #x::Array{Float32, 0} = zeros()


    y0::Float32 = 1 # Initial income value
    y::Array{Float32} =y0 * exp.(μ*t)  # Pre-calculate income process

    a0::Float32 = 100 # Initial wealth

    quad_eps::Float32 = 0.01 # Point to glue the quadtratic onto the CRRA
    quad_p::Float32 = 15 # Parameter deciding the steepness of the quadtratic glue
end

function crra(c, gamma) #c should probably be x since we apply these functions to consumpton and terminal wealth
    if c < 0
        return -Inf
    end
    if gamma == 1
        log(c)
    else
        c.^(1 - gamma) / (1 - gamma)
    end
end

function cara(c, a)
    (1 - exp(-a*c)) / a
end

function quad_interp(x, γ, ε, p)
    a = -2 * ε^p
    b = 1 / ε^γ + 1 / ε^(p - 1)
    c = - 1 / (2*ε^(p-2)) + ε^(1 - γ)*γ / (1 - γ)
    return a*x^2 + b*x + c
end

function crra_quad_interp(c, γ, ε, p)
  if c >= ε
      crra(c, γ)
  else
      quad_interp(c, γ, ε, p)
  end
end
