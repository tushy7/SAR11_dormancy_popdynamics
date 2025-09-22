using Plots

function plotSolution(sol::NTuple{6, Vector{Float64}}, withDormancy::Bool)
    A, Q, _, V, I, t = sol  # `_` ignores the nutrients 

    plot(t, A .+ Q .+ I,
         label="Total Cells", linewidth=2, yscale=:log10, ylim=(1e-10, 1e5),
         color=:black, grid=true, legend=:outertop)
    plot!(t, A, label="Active Cells", linestyle=:dash, linewidth=2, color=:red)
    plot!(t, V, label="Viruses", linestyle=:dot, linewidth=2, color=:purple)

    titleStr = withDormancy ? "SAR11 Dynamics with Dormancy" : "SAR11 Dynamics without Dormancy"
    title!(titleStr)
    xlabel!("Time (Days)")
    ylabel!("Cell Density (pmol cellular C/mL)")
    display(plot!())
end

