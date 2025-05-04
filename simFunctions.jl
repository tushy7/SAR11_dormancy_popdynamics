using DifferentialEquations, Plots
include("opt_inf_qui.jl")
include("parameters.jl")

function simulateGrowth(; withDormancy::Bool)
    # Parameters
    optParameters = copy(defaultOptParameters)
    phi = defaultPhi
    beta = defaultBeta
    simulationTime = 200.0

    # Set dormancy rate manually
    dormancyRate = 0.0
    if withDormancy
        dormancyRate = optParameters[2]
    end

    # Initial state: [Active, Quiescent, Nutrients, Viruses, Infected Quiescent]
    initialState = [initialCellCount, initialCellCount, nutrientReset, initialVirusLoad, 0.0]
    tspan = (0.0, simulationTime)

    # Params: (optParameters, phi, beta, modelParams, dormancyFlag)
    effectiveOptParams = [optParameters[1], dormancyRate]
    params = (effectiveOptParams, phi, beta, modelParams, withDormancy)

    # Solve
    prob = ODEProblem(opt_inf_qui_model!, initialState, tspan, params)
    sol = solve(prob, Tsit5(); reltol=1e-6, abstol=1e-6)

    # Extract results
    t = sol.t
    A = sol[1, :]
    Q = sol[2, :]
    V = sol[4, :]
    I = sol[5, :]

    return sol
end
