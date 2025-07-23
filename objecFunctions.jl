using Statistics
include("parameters.jl")
include("simFunctions.jl")
include("juliaPlotting.jl")

function objective(optParams, phi, dt)
	#run simulateGrowth for each input parameters
   	sol = simulateGrowth(optParams; phi=phi)
    active = sol[1] #save only active output

    # Calculate average over last 3 cycles 
    timePerCycle = round(Int, dt / 0.1) + 1
    endingCycles = 3 * timePerCycle
    
    if length(active) >= endingCycles && mean(active[end-endingCycles+1:end]) > 0
        avgActive = mean(active[end-endingCycles+1:end])
        return 1/avgActive
    else
        return 1e10 #really big number = bad
    end
end


#actual numerical optimization
using Optimization, OptimizationOptimJL

#optimize over [umax (max growth rate), q (quiescence rate)] for each phi and dt
function optimizeDormancy(phi, dt)

    obj(p, _) = objective(p, phi, dt) # two-arg function, Julia requires "_" 
    optf = OptimizationFunction(obj, Optimization.AutoForwardDiff()) #telling Optimizer what/how to optimize 
    prob = OptimizationProblem(optf, [0.66, 0.1]; lb=[0.1, 0.01], ub=[2.0, 1.0]) #like writing down the intial question with initial guesses etc.
    #[0.66, 0.1] = initial guess for q, lb = lower bounds for umax and q, ub = upper bounds.
    #these can be messed with

    sol = solve(prob, LBFGS(); maxiters=100, maxtime=500.0) #can start max iters lower for now if it’s slow

    return sol.u[2]  #return optimal dormancy rate q
end


#define range of interest for phi and dt, and loop through them
#these are your x and y axes in the heatmap
using CSV, DataFrames

function heatMap(phiRange::Tuple{Float64, Float64}, dtRange::Tuple{Float64, Float64}; outputFile::String="heatMapData.csv")
    phis = range(phiRange[1], phiRange[2], length=8) #section virus count range into 8 subsets
    dts = range(dtRange[1], dtRange[2], length=14) #section time range into 14 subsets

    # create empty DataFrame to store (phi, dt, optimal_quiescence)
    results = DataFrame(phi=Float64[], dt=Float64[], q=Float64[])

    # loop through every combo of phi and dt
    for phi in phis
        for dt in dts
            q_opt = optimizeDormancy(phi, dt)
            push!(results, (phi, dt, q_opt))
        end
    end

    CSV.write(outputFile, results)
    println("Done")
end
