using UnPack

function opt_inf_qui_model!(dstate, state, params, t)
    # Unpack parameters
    optParameters, phi, beta, modelParams, includeDormancy = params

    # Unpack model constants
    @unpack uptakeAffinity, resuscitationRate, outflowRate,
            activeMortalityRate, quiescentMortalityRate, yieldCoefficient,
            nutrientInflowRate, virusMortalityRate, phiCell, phiVirus, burstScale = modelParams

    # Unpack state variables
    active = state[1]
    quiescent = state[2]
    nutrients = state[3]
    viruses = state[4]
    infectedQ = state[5]

    # Optimization parameters
    maxGrowthRate = optParameters[1]
    dormancyRate = 0.0
    if includeDormancy
        dormancyRate = optParameters[2]
    end

    # Monod kinetics for nutrient-limited growth
    ks = maxGrowthRate / uptakeAffinity
    specificGrowthRate = maxGrowthRate * nutrients / (ks + nutrients)

    # ODEs
    dstate[1] = specificGrowthRate * active -
                dormancyRate * active +
                resuscitationRate * quiescent -
                (outflowRate + activeMortalityRate) * active -
                phi * phiCell * active * viruses

    dstate[2] = dormancyRate * active -
                resuscitationRate * quiescent -
                (outflowRate + quiescentMortalityRate) * quiescent -
                phi * phiCell * quiescent * viruses

    dstate[3] = nutrientInflowRate * outflowRate -
                outflowRate * nutrients -
                (1 / yieldCoefficient) * specificGrowthRate * active

    dstate[4] = beta * burstScale * phi * phiCell * active * viruses -
                phi * phiVirus * active * viruses -
                phi * phiVirus * quiescent * viruses -
                phi * phiVirus * infectedQ * viruses -
                (outflowRate + virusMortalityRate) * viruses

    dstate[5] = phi * phiCell * quiescent * viruses -
                quiescentMortalityRate * infectedQ
end
