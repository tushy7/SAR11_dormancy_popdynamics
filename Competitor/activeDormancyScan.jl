using Plots, Printf, Statistics
include("simFunctions.jl")
include("optimizeDormancy.jl")

phi             = 1e-9
dt              = 7.0
umax            = 0.69
dilute_dormant  = false
total_days      = 1000.0
eval_window     = 300.0
nq              = 201
dilution_effect = 0.25 
if !isempty(ARGS)
    try
        dilution_effect = parse(Float64, ARGS[1])
    catch
        @warn "Could not parse ARGS[1] as Float64; using default dilution_effect=$(dilution_effect)"
    end
end

tail_mean(q) = begin
    A, _, _, _, _, _, t = simulateGrowth([umax, q];
        phi=phi, withDormancy=true,
        dilution_interval=dt, total_time=total_days,
        run_partial_last=true,
        dilute_dormant=dilute_dormant,
        dilEffect=dilution_effect)
    t_end = t[end]
    keep  = t .>= (t_end - eval_window)
    vals  = collect(A[keep])
    isempty(vals) ? NaN : mean(vals)
end

qs  = range(0.0, umax; length=nq)
mA  = [tail_mean(q) for q in qs]

qopt, fmin, _ = optimizeDormancy(phi, dt;
    growth=umax, bounds=(0.0, umax),
    total_days=total_days, eval_window=eval_window,
    dilute_dormant=dilute_dormant,
    dilution_effect=dilution_effect)   # <- pass through

mode_str = dilute_dormant ? "dilute A+Q" : "dilute A only"
plt = plot(qs, mA; lw=3, c=:royalblue, label="Mean Active (last $(Int(eval_window)) d)",
           xlabel="Dormancy rate q (day⁻¹)",
           ylabel=@sprintf("Mean Active over last %.0f d", eval_window))
vline!(plt, [qopt]; c=:black, lw=2, ls=:dash, label=@sprintf("q* = %.3g", qopt))
title!(plt, @sprintf("Active vs Dormancy  (ϕ=%.1e, Δt=%.2f d, umax=%.3g, %s, dilEffect=%.2f)",
                     phi, dt, umax, mode_str, dilution_effect))
display(plt)
