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
   Pkg.add("UnPack")
   Pkg.add("Plots")
   Pkg.add("DifferentialEquations")

  
3. You must run `include("filename.jl")` to load functions from that file into your session.
  - Example:
    ```julia
    include("simFunctions.jl")
    ```
  - Once included, all functions in the file are accessible globally within the REPL session.


## File Overview

| File              | Description |
|-------------------|-------------|
| `opt_inf_qui.jl`  | Defines the core **ODE model** describing SAR11 population dynamics. This includes how active cells, quiescent cells, nutrients, viruses, and infected quiescent cells change over time. |
| `parameters.jl`   | Stores all **default model parameters** including growth rates, mortality rates, viral infection parameters, nutrient levels, and simulation settings. |
| `plotting.jl`     | Includes the `plotSolution` function to generate plots from the simulation results. It takes in the solution object (`sol`) and whether dormancy was used, and outputs a log-scaled population dynamics plot. |
| `simFunctions.jl` | Contains the `simulateGrowth` function, which sets up and solves the ODE problem. It takes a `withDormancy` flag and returns the simulation result (`sol`) without plotting. |
| `simulations.jl`  | The **main runner script**. It defines `regGrowthSim(eventType)` which calls `simulateGrowth` and then passes the result to `plotSolution`. |
