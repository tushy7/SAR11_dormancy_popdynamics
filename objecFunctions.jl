include("parameters.jl")
include("simFunctions.jl")
include("plotting.jl")

#set up the objective function
function objective(optParams, phi, dt)
	#run simulateGrowth for each input parameters
   	sol = simulateGrowth(optParams; phi=phi, simulationTime=dt*Dilutions)
	#save only active output
    	active = sol[1, :] #or whatever the output actually is

    if active[end] > 0
        return active[end]
    else
        return 1e10 #really big number = bad
    end
end


#actual numerical optimization
using Optimization, OptimizationOptimJL

#optimize over [umax, q] for each phi and dt
function optimizeDormancy(phi, dt)

    obj = (p, _) -> objective(p, phi, dt)  # two-arg function
    optf = OptimizationFunction(obj, Optimization.AutoForwardDiff())
    prob = OptimizationProblem(optf, [0.66, 0.1]; lb=[0.0, 0.0], ub=[2.0, 3.0])
    #[0.66, 0.1] = initial guess, lb = lower bounds for phi, dt, ub = upper bounds.
    #these can be messed with

    sol = solve(prob, BFGS(); maxiters=60) #can start max iters lower for now if it’s slow

    return sol.u[2]  #return optimal dormancy rate q
end


#define range of interest for phi and dt, and loop through them
#these are your x and y axes in the heatmap
using CSV
using DataFrames


function heatMap(phiRange::Tuple{Float64, Float64}, dtRange::Tuple{Float64, Float64}; outputFile::String="heatMapData.csv")
    phis = range(phiRange[1], phiRange[2], length=8)
    dts = range(dtRange[1], dtRange[2], length=14)

    # Create an empty DataFrame to store (phi, dt, optimal_q)
    results = DataFrame(phi=Float64[], dt=Float64[], q=Float64[])

    # loop through every combo of phi and dt
    for phi in phis
        for dt in dts
            q_opt = optimizeDormancy(phi, dt)
            push!(results, (phi, dt, q_opt))
        end
    end

    # Write to CSV
    CSV.write(outputFile, results)
    println("Done")
end
