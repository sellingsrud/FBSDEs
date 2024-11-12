using Plots
using LaTeXStrings

include("../src/deterministic_analytical.jl")

function plot_analytical(p)
    ct_b, yt_b, at_b, aT_b = optimal_bequest(p)
    ct_nb, yt_nb, at_nb = optimal_nobequest(p)


    p1 = plot(p.t, at_b, xlabel=L"$t$",label=L"$w_t$", legend=:topright) # Wealth
    p1 = plot!(p1, p.t, at_nb, label = L"$w_t^{nb}$", linestyle=:dash)

    p2 = plot(p.t, ct_b, xlabel=L"$t$", label=L"$c_t$",linestyle=:dash,legend=:topright) # Consumption
    p2 = plot!(p2, p.t, ct_nb, label = L"$c_t^{nb}$", linestyle=:dash)

    pg = plot(p1,p2, size=(1200, 800),linewidth=2, margins = 1.25Plots.cm)
    savefig("analytical_solution_deterministic_model.png")
end

plot_analytical(Param(;γ=0.5))
