using Roots
using ForwardDiff
using Roots

include("utils.jl")


# Use automatic differentation to use gradient method instead of Bisection
#D(f) = x -> ForwardDiff.derivative(f, float(x))


function aT(x, p::Param)
    (;γ,ρ,r,ψ,μ,λ,ε,y0,a0,t,T)=p

    par0 = ρ-(1-γ)*r
    par = (γ/par0)*exp(-((r-ρ)/γ)*T)*(1-exp(-(par0/γ)*T))
    hT = (y0/(μ-r)) * (exp((μ-r)*T)-1)
    res =  a0 + hT - λ^(-1/γ)*par*x^(ε/γ)  - exp(-r*T)*x
end

function optimal_bequest(p::Param)
    (;γ,ρ,r,ψ,μ,λ,ε,y0,a0,t,T)=p

    f(x) = aT(x, p)

    par0 = ρ-(1-γ)*r
    factor = (γ/par0)*exp(-((r-ρ)/γ)*T)

    hT = (y0/(μ-r)) * (exp((μ-r)*T)-1)
    aT_0 = (a0 + hT) / exp(-r*T)

    aT_upper_bound = exp.(r*t) .*(a0 .+ hT)

    aT_sol = find_zero(x -> aT(x, p), (0, 10000), Roots.Bisection())
    

    ct = exp.(-((r-ρ)/γ)*(T.-t)) *λ^(-1/γ)*aT_sol^(ε/γ)
    yt = y0*exp.(μ*t)
    ht = (y0/(μ-r)) * (exp.((μ-r).*t).-1) # human capital from 0 to t
    cum_ct = factor*(1 .-exp.(-(par0/γ).*t))*λ^(-1/γ)*aT_sol^(ε/γ)
    at = exp.(r*t) .*(a0 .+ ht .- cum_ct)

    return ct,yt,at,aT
end

function optimal_nobequest(cp)
    (;γ,ρ,r,ψ,μ,y0,a0,t,T) = cp

    h0 = (y0/(μ-r)) * (exp((μ-r)*T)-1)
    c0 = (r-ψ)/(1-exp(-(r-ψ)*T)) *(a0+h0)
    ct = c0*exp.(-((ρ-r)/γ).*t) 
    yt = y0*exp.(μ*t)  # Income
    ht = (y0/(μ-r)) * (exp.((μ-r)*t).-1) # human capital from 0 to t
    cum_ct = (c0/(r-ψ)) * (1 .- exp.(-(r-ψ)*t))
    at = exp.(r*t) .*(a0 .+ ht - cum_ct)

    return ct,yt,at
end