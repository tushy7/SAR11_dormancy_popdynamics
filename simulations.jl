include("parameters.jl")
include("simFunctions.jl")
include("plotting.jl")

function regGrowthSim(eventType::Int)
    if eventType == 0  # with dormancy
        sol = simulateGrowth(withDormancy=true)
        plotSolution(sol, true)
    else  # without dormancy
        sol = simulateGrowth(withDormancy=false)
        plotSolution(sol, false)
    end
end
