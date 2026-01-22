defaultOptParameters = [0.69, 0.1]  # [maxGrowthRate, dormancyRate]

# Virus parameters
defaultPhi = 1.0e-7              # viral adsorption rate (mL/day)
defaultBeta = 10.0                 # virus burst size (mL/day)
initialVirusLoad = 1.0727          # starting virus particles if virus is present (mL/day)

# Nutrient settings 
nutrientReset = 10000.0    # initial nutrient concentration

# Dilution Simulation settings
Range = 300
dt = 5
dilEffect = 0.12

# Cell concentrations
initialCellCount = 47.3606         # starting number for active & quiescent cells (pmol/mL)
initialCompetitorCount = initialCellCount   # starting number for competitor cells (pmol/mL)
includeCompetitors = true                   # optional toggle
#initialCompetitorCount = 0         # starting number for active & quiescent cells (pmol/mL)


extinctionThreshold = 5.41e-4  # cells (units: same as active/quiescent)


modelParams = Dict(
    :uptakeAffinity => 2.67, #pmol*day/mL (Maybe should be mL/(pmol*day))
    :resuscitationRate => 0.1, # 1/day
    :outflowRate => 0.0, # 1/day
    :activeMortalityRate => 0.2, # 1/day
    :quiescentMortalityRate => 0.02, # 1/day
    :sporeMortalityRate => 0.02, # 1/day
    :yieldCoefficient => 0.6, # not sure yet? 
    :nutrientInflowRate => 5.0, # 1/day
    :virusMortalityRate => 0.2, # 1/day
    :phiCell => 1 / 5.41e-3, # mL/day
    :phiVirus => 1 / 5.41e-3, # mL/day
    :burstScale => 1.0727 # infected_cell/pmol_cell_carbon 
)

defaultPhiMain = defaultPhi
defaultBetaMain = defaultBeta
initialVirusLoadMain = initialVirusLoad

defaultPhiMenace = defaultPhi
defaultBetaMenace = defaultBeta
initialVirusLoadMenace = initialVirusLoad
