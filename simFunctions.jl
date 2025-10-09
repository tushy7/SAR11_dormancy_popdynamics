using DifferentialEquations, Plots
include("opt_inf_qui.jl")
include("parameters.jl")

using OrdinaryDiffEq, DiffEqCallbacks

function simulateGrowth(optParameters; phi::Float64,
        withDormancy::Bool=true,
        dilution_interval::Float64=default_dilution_interval,
        total_time::Float64=default_total_time,
        run_partial_last::Bool=true,
        dilute_dormant::Bool = true,
        dilEffect::Float64=dilEffect,
        customModelParams::Union{Nothing,Dict}=nothing)

    beta = defaultBeta
    dormancyRate = withDormancy ? optParameters[2] : 0.0
    initialState = [initialCellCount, initialCellCount, nutrientReset, initialVirusLoad, 0.0,
        includeCompetitors ? initialCompetitorCount : 0.0]

    num_full = floor(Int, total_time / dilution_interval)
    rem_time = total_time - num_full * dilution_interval

    tspan  = (0.0, dilution_interval)
    mp = isnothing(customModelParams) ? modelParams : customModelParams
    params = ([optParameters[1], dormancyRate], phi, beta, mp, withDormancy)
    
    has_negative(u, t, integrator) = any(u .< 0.0)

    # Positivity guard (prevents negative excursions blowing up other terms)
    function clamp_negatives!(integrator)
        u = integrator.u
        @inbounds for i in eachindex(u)
            if u[i] < 0.0
                u[i] = 0.0
            end
        end
    end
    cb = DiscreteCallback(has_negative, clamp_negatives!)
    
    # Stiff algorithm + headroom + dt floor
    alg  = KenCarp4()    # TRBDF2() also works well here
    rtol = 1e-4
    atol = 1e-8
    mxi  = 2_000_000
    dmin = 1e-12

    save_step = min(0.1, dilution_interval/100)

    # First cycle
    prob = ODEProblem(opt_inf_qui_model!, initialState, tspan, params)
    sol  = solve(prob, alg; saveat=0:save_step:dilution_interval,
                 reltol=rtol, abstol=atol, maxiters=mxi, dtmin=dmin, callback=cb)

    actDil = copy(sol[1,:]); quiDil = copy(sol[2,:]); nutDil = copy(sol[3,:])
    virDil = copy(sol[4,:]); infDil = copy(sol[5,:]); menDil = copy(sol[6,:]); time  = copy(sol.t)

    # Full cycles
    for i in 1:(num_full - 1)
        act_mult = dilEffect
        men_mult = dilEffect
        qui_mult = dilute_dormant ? dilEffect : 1.0
        dilState = [actDil[end]*act_mult, quiDil[end]*qui_mult, nutrientReset,
                    virDil[end]*dilEffect, infDil[end]*dilEffect, menDil[end]*men_mult]
        total_main = dilState[1] + dilState[2]  # A + Q + infected-dormant

        if total_main < extinctionThreshold
            # wipe the main organism entirely (no resuscitation later)
            dilState[1] = 0.0   # Active
            dilState[2] = 0.0   # Quiescent
        end
        if i > 2 && dilState[6] < extinctionThreshold; dilState[6] = 0.0; end

        prob  = ODEProblem(opt_inf_qui_model!, dilState, tspan, params)
        newSol = solve(prob, alg; saveat=0:save_step:dilution_interval,
                       reltol=rtol, abstol=atol, maxiters=mxi, dtmin=dmin, callback=cb)

        append!(actDil, newSol[1,:]); append!(quiDil, newSol[2,:]); append!(nutDil, newSol[3,:])
        append!(virDil, newSol[4,:]); append!(infDil, newSol[5,:]); append!(menDil, newSol[6,:]); append!(time, newSol.t .+ time[end])

        checkExtinct!(actDil, extinctionThreshold); checkExtinct!(quiDil, extinctionThreshold); checkExtinct!(infDil, extinctionThreshold); checkExtinct!(menDil, extinctionThreshold)
    end

    # Partial last cycle (if any)
    if run_partial_last && rem_time > 1e-9
        act_mult = dilEffect
        men_mult = dilEffect
        qui_mult = dilute_dormant ? dilEffect : 1.0
        dilState = [actDil[end]*act_mult, quiDil[end]*qui_mult, nutrientReset,
                    virDil[end]*dilEffect, infDil[end]*dilEffect, menDil[end]*men_mult]
        
        total_main = dilState[1] + dilState[2]  # A + Q + infected-dormant

        if total_main < extinctionThreshold
            # wipe the main organism entirely (no resuscitation later)
            dilState[1] = 0.0   # Active
            dilState[2] = 0.0   # Quiescent
        end

        if dilState[6] < extinctionThreshold; dilState[6] = 0.0; end

        prob_last = ODEProblem(opt_inf_qui_model!, dilState, (0.0, rem_time), params)
        newSol = solve(prob_last, alg; saveat=0:save_step:rem_time,
                       reltol=rtol, abstol=atol, maxiters=mxi, dtmin=dmin, callback=cb)

        append!(actDil, newSol[1,:]); append!(quiDil, newSol[2,:]); append!(nutDil, newSol[3,:])
        append!(virDil, newSol[4,:]); append!(infDil, newSol[5,:]); append!(menDil, newSol[6,:]); append!(time, newSol.t .+ time[end])
    end

    return (actDil, quiDil, nutDil, virDil, infDil, menDil, time)
end


function checkExtinct!(trajectory::AbstractVector, threshold::Real)
    indices = findall(x -> x <= threshold, trajectory)
    if !isempty(indices)
        firstIdx = minimum(indices) 
        trajectory[firstIdx:end] .= zero(eltype(trajectory))
    end
end



function scanDormancyOld(phi, dt; qRange=range(0.01, 10.0, length=100))
    results = []
    for q in qRange
        optParams = [0.69, q]  # fixed max growth rate for now
        sol = simulateGrowth(optParams; phi=phi, withDormancy=true)
        active = sol[1]
        timePerCycle = round(Int, dt / 0.1) + 1
        endingCycles = 100 * timePerCycle
        if length(active) >= endingCycles && mean(active[end-endingCycles+1:end]) > 0
            avgActive = mean(active[end-endingCycles+1:end])
        else
            avgActive = 0.0  # or NaN
        end
        push!(results, (q, avgActive))
    end
    
    q = [result[1] for result in results if result[1] <= 2.5]
    avgActive = [result[2] for result in results if result[1] <= 2.5]

    p = plot(q, avgActive, linestyle=:dash, linewidth=2, color=:red)
    xlabel!("Quiescent Rate")
    ylabel!("Average Active Cell")
    title!("phi = $phi, dt = $dt")
    display(p)
end
using CSV, DataFrames
using Plots

"""
    scanDormancyFromCSV(csv_path, phi_idx, dt_idx;
                        qrange=range(0.0, 3.0; length=121),
                        growth=0.69,
                        total_days=500.0,
                        eval_window=300.0)

Picks φ and Δt by **index** from the unique sorted values in `csv_path`,
runs your existing `scanDormancy(φ, Δt, ...)`, and draws a vertical line
at the optimized q from the CSV.
"""
function scanDormancyFromCSV(csv_path::String, phi_idx::Int, dt_idx::Int;
    qrange=range(0.0, 3.0; length=501),
    growth=0.69,
    total_days=500.0,
    eval_window=300.0,
)
    df   = CSV.read(csv_path, DataFrame)
    phis = sort(collect(unique(df.phi)))
    dts  = sort(collect(unique(df.dt)))

    φ  = phis[phi_idx]
    Δt = dts[dt_idx]

    # run your existing scanner (the one you already have)
    qs, mtail, L = scanDormancy(φ, Δt;
        qrange=qrange, growth=growth,
        total_days=total_days, eval_window=eval_window)

    # overlay vertical line for the CSV's optimal q
    row = first(filter(r -> r.phi == φ && r.dt == Δt, df))
    qopt = row.q
    vline!([qopt]; lw=2, ls=:dot, label="opt q=$(round(qopt, digits=3))")
    display(current())

    return (φ, Δt, qopt)
end
