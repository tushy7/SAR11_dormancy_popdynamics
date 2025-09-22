include("parameters.jl")
include("growthSimulator.jl")
include("plotResults.jl")

function regGrowthSim(eventType::Int)
    # Default parameters
    optParameters = copy(defaultOptParameters)
    phi = defaultPhi

    if eventType == 0  # with dormancy
        sol = simulateGrowth(optParameters; phi=phi, withDormancy=true)
        plotSolution(sol, true)
    else  # without dormancy
        optParameters[2] = 0.0  # set dormancy rate to 0
        sol = simulateGrowth(optParameters; phi=phi, withDormancy=false)
        plotSolution(sol, false)
    end
end
