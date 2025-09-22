using DifferentialEquations, Plots
include("opt_inf_qui.jl")
include("parameters.jl")

using OrdinaryDiffEq, DiffEqCallbacks

function simulateGrowth(optParameters; phi::Float64,
        withDormancy::Bool=true,
        dilution_interval::Float64=default_dilution_interval,
        total_time::Float64=default_total_time,
        run_partial_last::Bool=true)

    beta = defaultBeta
    dormancyRate = withDormancy ? optParameters[2] : 0.0
    initialState = [initialCellCount, initialCellCount, nutrientReset, initialVirusLoad, 0.0]

    num_full = floor(Int, total_time / dilution_interval)
    rem_time = total_time - num_full * dilution_interval

    tspan  = (0.0, dilution_interval)
    params = ([optParameters[1], dormancyRate], phi, beta, modelParams, withDormancy)

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
    virDil = copy(sol[4,:]); infDil = copy(sol[5,:]); time  = copy(sol.t)

    # Full cycles
    for i in 1:(num_full - 1)
        dilState = [actDil[end]*dilEffect, quiDil[end]*dilEffect, nutrientReset,
                    virDil[end]*dilEffect, infDil[end]*dilEffect]
        if i > 2 && dilState[1] < extinctionThreshold; dilState[1] = 0; end

        prob  = ODEProblem(opt_inf_qui_model!, dilState, tspan, params)
        newSol = solve(prob, alg; saveat=0:save_step:dilution_interval,
                       reltol=rtol, abstol=atol, maxiters=mxi, dtmin=dmin, callback=cb)

        append!(actDil, newSol[1,:]); append!(quiDil, newSol[2,:]); append!(nutDil, newSol[3,:])
        append!(virDil, newSol[4,:]); append!(infDil, newSol[5,:]); append!(time, newSol.t .+ time[end])

        checkExtinct!(actDil, extinctionThreshold); checkExtinct!(quiDil, extinctionThreshold); checkExtinct!(infDil, extinctionThreshold)
    end

    # Partial last cycle (if any)
    if run_partial_last && rem_time > 1e-9
        dilState = [actDil[end]*dilEffect, quiDil[end]*dilEffect, nutrientReset,
                    virDil[end]*dilEffect, infDil[end]*dilEffect]
        prob_last = ODEProblem(opt_inf_qui_model!, dilState, (0.0, rem_time), params)
        newSol = solve(prob_last, alg; saveat=0:save_step:rem_time,
                       reltol=rtol, abstol=atol, maxiters=mxi, dtmin=dmin, callback=cb)

        append!(actDil, newSol[1,:]); append!(quiDil, newSol[2,:]); append!(nutDil, newSol[3,:])
        append!(virDil, newSol[4,:]); append!(infDil, newSol[5,:]); append!(time, newSol.t .+ time[end])
    end

    return (actDil, quiDil, nutDil, virDil, infDil, time)
end


function checkExtinct!(trajectory::AbstractVector, threshold::Real)
    indices = findall(x -> x <= threshold, trajectory)
    if !isempty(indices)
        firstIdx = minimum(indices) 
        trajectory[firstIdx:end] .= zero(eltype(trajectory))
    end
end



function scanDormancy(phi, dt; qRange=range(0.01, 10.0, length=100))
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
