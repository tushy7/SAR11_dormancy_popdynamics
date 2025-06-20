using DifferentialEquations, Plots
include("opt_inf_qui.jl")
include("parameters.jl")

function simulateGrowth(optParameters; phi::Float64, simulationTime::Float64, withDormancy::Bool=true)
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
    sol = solve(prob, Rodas5(); reltol=1e-6, abstol=1e-6)

    # Dilute results
    newSol = dilute(sol, tspan, params)

    return newSol
end

function dilute(sol::ODESolution, tspan, params)
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
        newSol = solve(prob, Rodas5(); reltol=1e-6, abstol=1e-6)

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

function checkExtinct!(trajectory::Vector{Float64}, threshold::Float64)
    i = findfirst(x -> x <= threshold, trajectory)
    if !isnothing(i)
        trajectory[i:end] .= 0.0
    end
end

