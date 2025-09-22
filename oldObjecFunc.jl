using Statistics
include("parameters.jl")
include("growthSimulator.jl")
include("plotResults.jl")

function objective(optParams, phi, dt) 
    #run simulateGrowth for each input parameters   
    sol = simulateGrowth(optParams; phi=phi)
    active = sol[1] #save only active output 
    
   # Make sure there's enough cycles to analyze, then check if population mean is above 0 to ensure they survived
    if length(active) >= 3001 && mean(active[end-2999:end]) > 0
        avgActive = mean(active[end-2999:end])
        return 1/avgActive
    else
        return 1e10 #really big number = bad
    end
end

#actual numerical optimization 
using Optimization, OptimizationOptimJL 
#optimize over [umax (max growth rate), q (quiescence rate)] for each phi and dt
function optimizeDormancy(phi, dt)
    obj(p, _) = objective(p, phi, dt)
    
    optf = OptimizationFunction(obj, Optimization.AutoForwardDiff())
    prob = OptimizationProblem(optf, [0.66, 0.1]; lb=[0.1, 0.01], ub=[2.0, 5.0])
    sol = solve(prob, LBFGS(); maxiters=100, maxtime=500.0)
    
    
    optParams = sol.u # Get the optimal parameters
    
    solGrowth = simulateGrowth(optParams; phi=phi) # Calculate the  active population
    active = solGrowth[1]
    
    if length(active) >= 3001 && mean(active[end-2999:end]) > 0
        avgActive = mean(active[end-2999:end])
    else
        avgActive = 0.0
    end
    
    return sol.u[2], avgActive  # return both dormancy rate and active population
end

#define range of interest for phi and dt, and loop through them
#these are your x and y axes in the heatmap
using CSV, DataFrames

function heatMap(phiRange::Tuple{Float64, Float64}, dtRange::Tuple{Float64, Float64}; outputFile::String="heatMapData.csv")
    phis = range(phiRange[1], phiRange[2], length=8)
    dts = range(dtRange[1], dtRange[2], length=14)
    
    # Create DataFrame to store all three values
    results = DataFrame(phi=Float64[], dt=Float64[], q=Float64[], avgActive=Float64[])
    
    for phi in phis
        for dt in dts
            q_opt, avg_active = optimizeDormancy(phi, dt)  # Now returns both values
            push!(results, (phi, dt, q_opt, avg_active))
        end
    end
    
    CSV.write(outputFile, results)
    println("Done")
    return results 
end