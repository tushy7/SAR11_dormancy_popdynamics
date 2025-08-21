using DifferentialEquations, Plots
include("opt_inf_qui.jl")
include("parameters.jl")

function simulateGrowth(optParameters; phi::Float64, withDormancy::Bool=true)
    # Parameters
    beta = defaultBeta

    # Set dormancy rate manually
    dormancyRate = 0.0
    if withDormancy
        dormancyRate = optParameters[2]
    end

    # Initial state: [Active, Quiescent, Nutrients, Viruses, Infected Quiescent]
    initialState = [initialCellCount, initialCellCount, nutrientReset, initialVirusLoad, 0.0]
    tspan = (0.0, dt)

    # Params: (optParameters, phi, beta, modelParams, dormancyFlag)
    effectiveOptParams = [optParameters[1], dormancyRate]
    params = (effectiveOptParams, phi, beta, modelParams, withDormancy)

    # Solve
    prob = ODEProblem(opt_inf_qui_model!, initialState, tspan, params)
    sol = solve(prob, Rodas5(); saveat=0:0.1:dt, reltol=1e-6, abstol=1e-6) #saveat creates equal time divides between 0 to dt 

    actDil = copy(sol[1, :])
    quiDil = copy(sol[2, :])
    nutDil = copy(sol[3, :])
    virDil = copy(sol[4, :])
    infDil = copy(sol[5, :])
    time = copy(sol.t)

    for i in 1:(Dilutions - 1)
        dilState = [actDil[end]*dilEffect, quiDil[end]*dilEffect, nutrientReset, 
                    virDil[end]*dilEffect, infDil[end]*dilEffect]

        if i > 2 && dilState[1] < 5.41e-4 
            dilState[1] = 0
        end 

        prob = ODEProblem(opt_inf_qui_model!, dilState, tspan, params)
        newSol = solve(prob, Rodas5(); saveat=0:0.1:dt, reltol=1e-6, abstol=1e-6, maxiters=1000)

        append!(actDil, newSol[1,:])
        append!(quiDil, newSol[2,:])
        append!(nutDil, newSol[3,:])
        append!(virDil, newSol[4,:])
        append!(infDil, newSol[5,:])
        append!(time, newSol.t .+ time[end])

        checkExtinct!(actDil, 5.41e-4)
        checkExtinct!(quiDil, 5.41e-4)
        checkExtinct!(infDil, 5.41e-4)
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
