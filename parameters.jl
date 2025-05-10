# Model parameter defaults for SAR11 viral/dormancy simulations

# Growth and dormancy parameters
defaultOptParameters = [0.66, 0.1]  # [maxGrowthRate, dormancyRate]

# Virus parameters
defaultPhi = 1e-6                  # viral adsorption rate (mL/day)
defaultBeta = 10.0                 # virus burst size (mL/day)
initialVirusLoad = 1.0727          # starting virus particles if virus is present (mL/day)

# Nutrient settings 
nutrientReset = 100000.0    # initial nutrient concentration

# Simulation settings
simulationTime = 200.0       # total duration (days)

# Cell concentrations
initialCellCount = 47.3606         # starting number for active & quiescent cells (pmol/mL)

#Used for time stamps
Dilutions = 30


# Internal model constants used in ODEs
modelParams = Dict(
    :uptakeAffinity => 2.67, #pmol*day/mL (Maybe should be mL/(pmol*day))
    :resuscitationRate => 0.2, # 1/day
    :outflowRate => 0.0, # 1/day
    :activeMortalityRate => 0.04, # 1/day
    :quiescentMortalityRate => 0.01, # 1/day
    :yieldCoefficient => 0.6, # not sure yet? 
    :nutrientInflowRate => 5.0, # 1/day
    :virusMortalityRate => 0.2, # 1/day
    :phiCell => 1 / 5.41e-3, # mL/day
    :phiVirus => 1 / 5.41e-3, # mL/day
    :burstScale => 1.0727 # infected_cell/pmol_cell_carbon 
)
