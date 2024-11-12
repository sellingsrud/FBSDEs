using Plots #; pgfplotsx()
using LaTeXStrings

include("../src/deterministic_analytical.jl")

function plot_analytical(p)
    ct_b, yt_b, at_b, aT_b = optimal_bequest(p)
    ct_nb, yt_nb, at_nb = optimal_nobequest(p)


    p1 = plot(p.t, at_b, xlabel=L"$t$",label=L"$w_t$ - bequest", legend=:topright, color=:blue) # Wealth
    p1 = plot!(p1, p.t, at_nb, label = L"$w_t$ - no bequest", linestyle=:solid, color=:orange)

    p2 = plot(p.t, ct_b, xlabel=L"$t$", label=L"$c_t$ - bequest",linestyle=:solid,legend=:topright, color=:blue) # Consumption
    p2 = plot!(p2, p.t, ct_nb, label = L"$c_t$ - no bequest", linestyle=:solid, color=:orange)

    pg = plot(p1,p2, size=(1200, 800),linewidth=2, margins = 1.25Plots.cm)
    pg2 = plot(p1,p2, size=(1200, 800),linewidth=2, margins = 1.25Plots.cm)


    #pg2 = plot(p1,p2, size=(1200, 800),linewidth=2, margins = 1.25Plots.cm, tex_output_standalone = true)


    
    savefig(pg,"plots/analytical_solution_deterministic_model.png")
    #savefig(pg2,"plots/analytical_solution_deterministic_model.tex")
end

#Include comment of which parameters are admissible to enter in Param() argument below

plot_analytical(Param())


