using Plots, Printf
include("parameters.jl")
include("simFunctions.jl")

phi        = 1.0e-7
dt         = 5.0
dilEffect  = 0.1
diluteDormant = false
total_days = 500.0
umax       = defaultOptParameters[1]
q          = defaultOptParameters[2]
withDormancy = true
savepath   = nothing

q = 0.2

A, Q, N, V, IQ, M, V2, t = simulateGrowth([umax, q];
    phi=phi,
    withDormancy=withDormancy,
    dilution_interval=dt,
    total_time=total_days,
    run_partial_last=true,
    dilute_dormant=diluteDormant,
    dilEffect=dilEffect)

total_cells = A .+ Q .+ IQ


p1 = plot(t, A;  lw=2, label="Active (A)")
plot!(p1, t, Q;  lw=2, label="Dormant (Q)")
plot!(p1, t, IQ; lw=2, label="Infected (IQ)")
plot!(p1, t, M;  lw=2, label="Menace (M)")
plot!(p1, t, total_cells; lw=2, ls=:dash, label="A + Q + IQ")
xlabel!(p1, "Time (days)")
ylabel!(p1, "Cells (model units)")
title!(p1, @sprintf("Cells  (ϕ=%.1e, Δt=%.2f d, q=%.3g, dilEff=%.2f, %s)",
                    phi, dt, q, dilEffect, diluteDormant ? "dilute A+Q" : "dilute A only"))

p2 = plot(t, V; lw=2, label="Viruses (V)")
#plot!(p2, t, V2;  lw=2, ls=:dash, label="M Viruses")
xlabel!(p2, "Time (days)")
ylabel!(p2, "Viruses (model units)")
title!(p2, "Viruses")

pulses = 0.0:dt:total_days
for τ in pulses
    vline!(p1, [τ]; lw=1, ls=:dot, c=:gray, label=false)
    vline!(p2, [τ]; lw=1, ls=:dot, c=:gray, label=false)
end

plt = plot(p1, p2; layout=(2,1), size=(900,700), legend=:topright)

if savepath !== nothing
    savefig(plt, savepath)
    @info "Saved figure to $savepath"
end

display(plt)
