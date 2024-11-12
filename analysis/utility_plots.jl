using Plots
using LaTeXStrings

include("../src/utils.jl")


p = Param()

viol_rng = -1.0e-3:1.0e-5:1.0e-4
c_range = 0:0.01:5
c_range_smooth = -1.0e-16:1.0e-18:1.0e-16

plot_dict = (
    xlabel = L"$c$",
    label = L"$u(c)$",
    legend = :bottomright,
    size = (1200, 800),
    linewidth = 2,
    margins = 2Plots.cm
)

plot_crra = plot(c_range, crra.(c_range, p.γ); plot_dict...)
savefig(plot_crra, "../plots/crra_plot.png")

plot_exp = plot(c_range, exp_util.(c_range, 2); plot_dict...)
savefig(plot_exp, "exp_plot.png")

plot_smooth = plot(c_range_smooth, smooth_crra.(c_range_smooth, EPS64, p.γ); plot_dict...)
savefig(plot_smooth, "crra_plot_penalty.png")

plot_violation = plot(viol_rng, violation_penalty.(viol_rng, EPS64); plot_dict...)
savefig(plot_violation, "violation_penalty_plot.png")

plot_smooth_2 = plot(viol_rng, smooth_crra.(viol_rng, EPS64, p.γ); plot_dict...)
savefig(plot_smooth_2, "crra_plot_penalty_2.png")


function plot_analytical(p)
    ct_b, yt_b, at_b, aT_b = optimal_bequest(p)
end
