#timeSeriesPlot.jl
using DifferentialEquations, Plots
include("parameters.jl")
include("growthSimulator.jl")
include("objecFunctions.jl")

phi = 1.0e-7
dt  = 5.0
T   = 500.0

#optimizes locally
opt_any = optimizeDormancy(phi, dt)
opt = length(opt_any) == 2 ? opt_any : [defaultOptParameters[1], opt_any]

A, Q, N, V, IQ, t = simulateGrowth(opt; phi=phi, withDormancy=true,
                                   dilution_interval=dt, total_time=T)

#plot
p = plot(t, A, label="Active cells")
plot!(p, t, V, label="Viruses")

#shade
t1, t2 = maximum(t) - 300, maximum(t)
plot!(p, [t1, t2], seriestype=:vspan, alpha=0.15, label="mask window")

xlabel!(p, "Time (days)"); ylabel!(p, "Concentration")
display(p)
