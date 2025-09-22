# SAR11 Dormancy Population Dynamics

This repository contains a modular Julia implementation of a population dynamics model for SAR11 bacteria.

---

## Initial Steps

To run the simulations locally:

1. **Open the Julia REPL**
   - Use `Shift + Cmd + P` in VS Code
   - Select `Julia: Start REPL`

2. **Install required packages (once)**
   ```julia
   import Pkg
   Pkg.add("Brent")
   Pkg.add("CSV")
   Pkg.add("DifferentialEquations")
   Pkg.add("DiffEqCallbacks")
   Pkg.add("DataFrames")
   Pkg.add("Optimization")
   Pkg.add("OptimizationOptimJL")
   Pkg.add("optimize")
   Pkg.add("Options")
   Pkg.add("OrdinaryDiffEq")
   Pkg.add("Plots")
   Pkg.add("Statistics")
   Pkg.add("UnPack")
   

  
3. You must run `include("filename.jl")` to load functions from that file into your session.
  - Example:
    ```julia
    include("growthSimulator.jl")
    ```
  - Once included, all functions in the file are accessible globally within the REPL session.


## File Overview

| File              | Description |
|-------------------|-------------|
| `opt_inf_qui.jl`  | Defines the core **ODE model** describing SAR11 population dynamics. This includes how active cells, quiescent cells, nutrients, viruses, and infected quiescent cells change over time. |
| `parameters.jl`   | Stores all **default model parameters** including growth rates, mortality rates, viral infection parameters, nutrient levels, and simulation settings. |
| `plotResults.jl`     | Includes the `plotSolution` function to generate plots from the simulation results. It takes in the solution object (`sol`) and whether dormancy was used, and outputs a log-scaled population dynamics plot. |
| `growthSimulator.jl` | Contains the `simulateGrowth` function, which sets up and solves the ODE problem. It takes a `withDormancy` flag and returns the simulation result (`sol`) without plotting. |
| `runSimulator.jl`  | The **main runner script**. It defines `regGrowthSim(eventType)` which calls `simulateGrowth` and then passes the result to `plotSolution`. |
| `optimizeDormancy.jl` | Defines the **objective function**  which determines the optimal dormancy rate that maximizes active cell density for given virus adsorption rate (`phi`) and dilution interval (`dt`) |
| `main.jl` | Calls and runs optimizeDormancy.jl |
| `timeSeriesPlot.jl` | Plots SAR11 population against time |


