using Lux
using SciMLSensitivity
using ComponentArrays
using Zygote
using Plots
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

# Loss for constraint violation
# violation_penalty(x, eps) = eps.*x - eps

# Callback function prints results every nth epoch
function callback(p, l)
    push!(losses, l)
    if length(losses) % 10 == 0
        println("Current loss after $(length(losses)) iterations: $(losses[end])")
    end
    return false
end


function wealth(p, alpha, θ)
  (; r, Δt, n, N,y) = p
  w = p.a0 # initial wealth

  for i in 1:N
      w_new = w[i] .+ r * w[i] * Δt .+ y[i] *Δt .- alpha([i,w[i]], θ, st)[1][1]*Δt
      w = vcat(w,w_new')
  end

  return w
end

function running_utility(p, w, alpha, θ)
  (;n,N,ρ,t,Δt,γ) = p

  r_util = 0

  for i in 1:N
    u = 0

    ct = alpha([i,w[i]],θ,st)[1][1]

    if ct <=0
      u = exp.(ρ *t[i])*violation_penalty(ct, EPS64)/Δt
    else
      u = u1(ct, p)
    end

      r_util = r_util .+ Δt * u * exp.(-ρ*t[i])
  end
  return r_util
end

# Terminal utility
function terminal_utility(p,wT)
  (; λ,ε,ρ,N,T) = p
  if wT<=0
      exp.(ρ*T)*violation_penalty(wT, EPS64)/λ
  else
      return exp.(-ρ*T).*λ.*u2(wT, p)
  end
end

function total_utility(θ)
  w = wealth(p, alpha, θ)


  total_utility = (terminal_utility(p,w[end]) .+ running_utility(p, w, alpha,θ))[1]

  return -total_utility
end

function plot_nn_analytical(p,pp_star)
    #Compare NN and analytic optimal wealth and consumption functions
    ct_b,yt_b,at_b,aT_b = optimal_bequest(p)

    w = wealth(p,alpha,pp_star)

    optimal_control = alpha([1,w[1]],pp_star, st)[1][1]

    for i in 1:p.N
        optimal_control = vcat(optimal_control,alpha([i,w[i]],pp_star,st)[1][1])
    end

    p1 = plot(p.t, w, xlabel="t", label=L"w_t NN",linestyle=:solid,legend=:topright)  # Wealth post-train no running
    p1 = plot!(p.t, at_b, xlabel="t",label=L"w_t Analytic",legend=:topright) # Wealth pre-train

    p2 = plot(p.t, optimal_control, xlabel="t", label=L"c_t NN",linestyle=:dash,legend=:bottomleft)
    p2 = plot!(p.t, ct_b, xlabel="t", label=L"c_t Analytic",linestyle=:dash,legend=:bottomleft)

    pg = plot(p1,p2, size=(800, 400),linewidth=2, margins = 1.25Plots.cm)
    savefig(pg, "nn_deterministic_solution.png")

    println("w_T NN = $(w[end]) and w_T Analytic = $(aT_b)")

end

#TODO make a structure with Neural Network (hyper-)parameters
#TODO make a structure to hold the results of this train function
#TODO pass alpha, p, st etc to the loss, now they're global
losses = Float64[]
p = Param()

dims = 10
alpha = Chain(
    Dense(2, dims, relu),
    Dense(dims, dims, relu),
    Dense(dims, dims, relu),
    Dense(dims, 1)
)

# alpha_pretrain = deepcopy(alpha)

pp, st = Lux.setup(rng, alpha)
# pp_pretrain, st_pretrain = Lux.setup(rng, alpha_pretrain)

adtype = Optimization.AutoZygote()
optf = Optimization.OptimizationFunction((x, pp) -> total_utility(x), adtype)
optprob = Optimization.OptimizationProblem(optf, ComponentVector{Float64}(pp))

@info "Training with 0.1 learning"
res1 = Optimization.solve(optprob, ADAM(0.1), callback = callback, maxiters = 50)

@info "Training with 0.01 learning"
optprob2 = Optimization.OptimizationProblem(optf, res1.u)
res2 = Optimization.solve(optprob2, ADAM(0.01), callback = callback, maxiters = 250)

@info "Training with 0.001 learning"
optprob3 = Optimization.OptimizationProblem(optf, res2.u)
res3 = Optimization.solve(optprob3, ADAM(0.001), callback = callback, maxiters = 1000)

@info "Training with 0.0001 learning"
optprob4 = Optimization.OptimizationProblem(optf, res3.u)
res4 = Optimization.solve(optprob4, ADAM(0.0001), callback = callback, maxiters = 1000)

@info "Training with LBFGS"
optprob5 = Optimization.OptimizationProblem(optf, res4.u)
res5 = Optimization.solve(optprob5, Optim.BFGS(), callback = callback, maxiters = 1500)

# Optimal weights
pp_star = res5.u


plot_nn_analytical(p, pp_star)

# Use determinstic solution as initial values for the stochastic one. How much speed gain?
