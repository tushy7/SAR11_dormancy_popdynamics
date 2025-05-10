include("parameters.jl")
include("simFunctions.jl")
include("plotting.jl")

function regGrowthSim(eventType::Int)
    # Default parameters
    optParameters = copy(defaultOptParameters)
    phi = defaultPhi
    simTime = simulationTime  # from parameters.jl

    if eventType == 0  # with dormancy
        sol = simulateGrowth(optParameters; phi=phi, simulationTime=simTime, withDormancy=true)
        plotSolution(sol, true)
    else  # without dormancy
        optParameters[2] = 0.0  # set dormancy rate to 0
        sol = simulateGrowth(optParameters; phi=phi, simulationTime=simTime, withDormancy=false)
        plotSolution(sol, false)
    end
end
