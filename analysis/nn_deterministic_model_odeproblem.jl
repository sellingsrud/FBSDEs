# This file replaces our own EM scheme with a diff eq package

using Lux
using OrdinaryDiffEq
using DifferentialEquations
using SciMLSensitivity
using ComponentArrays
using Zygote
using Plots #; pgfplotsx()
using StableRNGs

using Optimization
using OptimizationOptimisers
using OptimizationOptimJL

using LaTeXStrings

include("../src/utils.jl")
include("../src/deterministic_analytical.jl")


# Set seed for reproducability
rng = StableRNG(1111)

u1(x, p::Param) = crra(x, p.γ)
u2(x, p::Param) = crra(x, p.ε)

function analytical_ct(aT, tt, p::Param)
    @unpack_Param p
    return exp.(-((r-ρ)/γ)*(T.-tt)) *λ^(-1/γ)*aT^(ε/γ)
end

function x!(dx, x, p, t)
    w, y, c = x
    dx[1] = w*r + y - c
    dx[2] = y*μ
    dx[3] = gradient(t -> analytical_ct(aT_sol, t, params), t) |> first
end

function diffusion!(dx, x, p, t)
    dx[1] = 0
    dx[2] = 0
    dx[3] = 0
end

function x_nn!(dx, x, p, t)
    w, y, c = x
    ps = ComponentArray(p, ax)
    control = ann([w, t], ps, st)[1][1]

    dx[1] = w*r + y - control
    dx[2] = y*μ
    dx[3] = control
end

function predict(θ)
    Array(solve(prob, Euler(), p = θ, dt = params.Δt))
end

function running_utility(p, ct)
  (;n,N,ρ,t,Δt,γ) = p

    println(length(ct))
  r_util = 0

  for i in 1:N
    u = 0

    u = crra_quad_interp(ct[i], γ, 0.01, 20)
    r_util = r_util .+ Δt * u * exp.(-ρ*t[i])
  end
  return r_util
end

# Terminal utility
function terminal_utility(p,wT)
  (; λ,ε,ρ,N,T) = p
  exp.(-ρ*T).*λ.*crra_quad_interp(wT, ε, 0.01, 20)
end

function total_utility(θ)
    x = predict(θ)
    w = x[1, :]
    c = x[3, :]

    total_utility = (terminal_utility(params,w[end]) .+ running_utility(params, c))[1]

    return -total_utility
end

cb = function (state, l; doplot = true)
    println(l)

    ps = ComponentArray(state.u, ax)

    if doplot
        # p = plot(solve(remake(prob, p = state.u), Tsit5(), saveat = 0.01),
            # ylim = (-6, 6), lw = 3)
        # plot!(p, ts, [first(first(ann([t], ps, st))) for t in ts], label = "u(t)", lw = 3)
        # display(p)
    end

    return false
end

losses = Float64[]
params = Param()
r = params.r
μ = params.μ
γ = params.γ
tspan = (params.t0, params.T)
ct, yt, at, aT_sol = optimal_bequest(params)
c0 = analytical_ct(aT_sol, Float32(0.0), params) |> Float32
x0 = [params.a0, params.y0, c0]

prob_det = ODEProblem(x!, x0, tspan, [])
res_analytical = Array(solve(prob_det, Vern9(), saveat=params.Δt))

prob_stoch = SDEProblem(x!, diffusion!, x0, tspan, [])
res_stoch = solve(prob_stoch, SRIW1(), dt = params.Δt, adaptive = false, abstol = 1e-10, reltol = 1e-10)

dims = 15
ann = Chain(
    Dense(2, dims, relu),
    Dense(dims, dims, relu),
    Dense(dims, 1)
)


ps, st = Lux.setup(rng, ann)
p = ComponentArray(ps)
θ, _ax = getdata(p), getaxes(p)
const ax = _ax


prob = ODEProblem(x_nn!, x0, tspan, θ)
# solve(prob, Vern9(), dt = params.Δt)
res= Array(solve(prob, Vern9(), saveat=params.Δt))

adtype = Optimization.AutoZygote()
optf = Optimization.OptimizationFunction((x, p) -> total_utility(x), adtype)

optprob = Optimization.OptimizationProblem(optf, θ)
res1 = Optimization.solve(
    optprob, OptimizationOptimisers.Adam(0.01), callback = cb, maxiters = 100)

optprob2 = Optimization.OptimizationProblem(optf, res1.u)
res2 = Optimization.solve(
    optprob2, OptimizationOptimJL.BFGS(), callback = cb, maxiters = 100)
