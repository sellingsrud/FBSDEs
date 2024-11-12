using Plots #; pgfplotsx()
using LaTeXStrings

include("../src/utils.jl")


p = Param()

viol_rng = -0.1:0.0001:0.5
c_neg_range = -0.1:0.01:5
c_range = -1.0:0.01:5
c_range_smooth = -1.0e-16:1.0e-18:1.0e-16

plot_dict = (
    xlabel = L"$c$",
    legend = :bottomright,
    size = (1200, 800),
    linewidth = 2,
    margins = 2Plots.cm,
    ylims = (-5, 5),
    xlim = (-1, 5)
)

plot_crra = plot(c_range, crra_talk.(c_range, p.γ), label=L"$u^{crra}(c)$"; plot_dict...)
savefig(plot_crra, "plots/crra_plot.png")

plot_exp = plot(c_range, exp_util.(c_range, 2), label=L"$u^{exp}(c)$"; plot_dict...)
savefig(plot_exp, "plots/exp_plot.png")

plot_moll = plot(c_range, mollified_crra.(c_range, p.γ, 0.01, 1), label=L"$u^{crra\_quad}(c)$"; plot_dict...)
savefig(plot_moll, "plots/crra_quad_interp.png")