using CSV, DataFrames

function plotSolutions(sol, withDormancy::Bool)
    filename = ""
    if withDormancy
        filename = "growth_with_dormancy.csv"
    else
        filename = "growth_without_dormancy.csv"
    end

    # Create a DataFrame from the ODE solution
    df = DataFrame(time = sol.t)
    for i in 1:length(sol.u[1])
        df[!, variableName(i)] = sol[i, :]
    end

    CSV.write(filename, df)
    println("Saved solution to $filename")
end

function variableName(i)
    # Replace with your real variable names if known
    names = ["Active", "Quiescent", "Nutrients", "Viruses", "Infected_Quiescent"]
    if i <= length(names)
        return names[i]
    else
        return "var_$i"
    end
end
