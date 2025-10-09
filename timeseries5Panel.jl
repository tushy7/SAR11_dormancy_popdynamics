using Plots, Printf, Statistics
include("simFunctions.jl")

phi             = 0.0
dt              = 3.0
umax            = 0.69
q               = 0.0
total_days      = 1000.0
dilute_dormant  = false
dilution_effect = 0.06
log_y           = false

A, Q, N, V, I, C, t = simulateGrowth([umax, q];
    phi=phi,
    withDormancy=true,
    dilution_interval=dt,
    total_time=total_days,
    run_partial_last=true,
    dilute_dormant=dilute_dormant,
    dilEffect=dilution_effect)

Total = A .+ Q .+ I

safe(v) = log_y ? map(x -> x > 0 ? x : NaN, v) : v
ys = log_y ? :log10 : :identity

title_main = @sprintf("ϕ=%.1e, Δt=%.2f d, q=%.3g, umax=%.3g, %s, dilEffect=%.2f",
                      phi, dt, q, umax, (dilute_dormant ? "dilute A+Q" : "dilute A only"),
                      dilution_effect)

lay = @layout [ a{0.32h}; b; c; d; e; f ]

p_all = plot(t, safe(A), lw=2, label="Active", yscale=ys, legend=:topright,
             title=title_main, ylabel="cells", xlabel="")
plot!(p_all, t, safe(Q),     lw=2, label="Dormant")
plot!(p_all, t, safe(V),     lw=2, label="Viruses")
plot!(p_all, t, safe(I),     lw=2, label="Infected")
plot!(p_all, t, safe(Total), lw=2, label="Total")

p_A = plot(t, safe(A), lw=2, yscale=ys, label=false, title="Active",   ylabel="cells",      xlabel="")
p_Q = plot(t, safe(Q), lw=2, yscale=ys, label=false, title="Dormant",  ylabel="cells",      xlabel="")
p_V = plot(t, safe(V), lw=2, yscale=ys, label=false, title="Viruses",  ylabel="particles",  xlabel="")
p_I = plot(t, safe(N), lw=2, yscale=ys, label=false, title="Nutrients", ylabel="cells",      xlabel="")
p_T = plot(t, safe(Total), lw=2, yscale=ys, label=false, title="Total Cells", ylabel="cells", xlabel="time (days)")

plt = plot(p_all, p_A, p_Q, p_V, p_I, p_T; layout=lay, size=(950, 1250))
display("image/png", plt)
