using CSV, DataFrames, Plots, Printf

data = CSV.read("heatmap_results_0.2.csv", DataFrame)

phis = sort(unique(data.phi))
dts  = sort(unique(data.dt))

qmat = [begin
    vals = data.q[(data.phi .== ϕ) .& (data.dt .== δ)]
    isempty(vals) ? NaN : vals[1]
end for ϕ in phis, δ in dts]

xidx = 1:length(dts)
yidx = 1:length(phis)
xlabels = [@sprintf("%.2f", x) for x in dts]
ylabels = [@sprintf("%.1e", y) for y in phis]

plt = heatmap(xidx, yidx, qmat;
        xlabel="# days between dilution",
        ylabel="Virus adsorption rate",
        c=:viridis,
        colorbar_title="Quiescence rate",
        xticks=(xidx, xlabels),
        yticks=(yidx, ylabels),
        xlims=(0.5, length(dts)+0.5),
        ylims=(0.5, length(phis)+0.5),
        aspect_ratio=:equal,
        title="Quiescence rate",
        size=(800,800))

display("image/png", plt)
