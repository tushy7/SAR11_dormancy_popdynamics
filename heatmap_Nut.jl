using CSV, DataFrames, Plots, Printf

phi = 0.0

df = CSV.read("qstar_dilution_map_0.0.csv", DataFrame)

effs = sort(unique(df.dilEffect))
dts  = sort(unique(df.dt))

qmat = [begin
    vals = df.q[(df.dilEffect .== e) .& (df.dt .== d)]
    isempty(vals) ? NaN : vals[1]
end for e in effs, d in dts]

xidx = 1:length(dts)
yidx = 1:length(effs)
xlabels = [@sprintf("%.2f", d) for d in dts]
ylabels = [@sprintf("%.2f", e) for e in effs]

plt = heatmap(xidx, yidx, qmat;
        xlabel="Δt (days between dilutions)",
        ylabel="dilution effect",
        c=:viridis, colorbar_title="q* (day⁻¹)",
        xticks=(xidx, xlabels),
        yticks=(yidx, ylabels),
        xlims=(0.5, length(xidx)+0.5),
        ylims=(0.5, length(yidx)+0.5),
        aspect_ratio=:equal,
        title=@sprintf("q* vs dt and dilEffect (ϕ=%.1e, umax=%.2f)", phi, umax),
        size=(950,820))

display("image/png", plt)
