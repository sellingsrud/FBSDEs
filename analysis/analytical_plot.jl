using Plots
using LaTeXStrings

include("../src/deterministic_analytical.jl")

function plot_analytical(p)
    ct_b, yt_b, at_b, aT_b = optimal_bequest(p)

    p1 = plot(p.t, at_b, xlabel=L"$t$",label=L"$w_t$",linestyle=:dash, legend=:topright) # Wealth

    p2 = plot(p.t, ct_b, xlabel=L"$t$", label=L"$c_t$",linestyle=:dash,legend=:topright) # Consumption

    pg = plot(p1,p2, size=(1200, 800),linewidth=2, margins = 1.25Plots.cm)
    savefig("plots/analytical_solution_deterministic_model.png")
end

#Include comment of which parameters are admissible to enter in Param() argument below

plot_analytical(Param(;λ=5))


