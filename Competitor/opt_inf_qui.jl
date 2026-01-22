using UnPack

function opt_inf_qui_model!(dstate, state, params, t)
    optParameters, phi_main, beta_main, phi_menace, beta_menace, modelParams, includeDormancy, isSporeFormer = params

    @unpack uptakeAffinity, resuscitationRate, outflowRate,
        activeMortalityRate, quiescentMortalityRate, sporeMortalityRate, yieldCoefficient,
        nutrientInflowRate, virusMortalityRate, phiCell, phiVirus, burstScale = modelParams

    active       = state[1]
    quiescent    = state[2]
    nutrients    = state[3]
    virus_main   = state[4]
    infectedQ    = state[5]
    menace_raw   = state[6]
    virus_menace = state[7]

    menace = max(menace_raw, 0.0)  # menace can't go below 0


    maxGrowthRate = optParameters[1]
    dormancyRate  = (includeDormancy && length(optParameters) >= 2) ? optParameters[2] : 0.0

    ks = maxGrowthRate / uptakeAffinity
    #specificGrowthRate = maxGrowthRate * nutrients / (ks + nutrients)
    N  = max(nutrients, 0.0)
    specificGrowthRate = maxGrowthRate * N / (ks + N)

    if isSporeFormer
        qMort = sporeMortalityRate
        viralSusceptibility = 0.0  # spores are impervious to new infections
    else
        qMort = quiescentMortalityRate
        viralSusceptibility = 1.0  # normal infection rate
    end

    # ODEs
    #Active
    dstate[1] = specificGrowthRate * active -
                dormancyRate * active +
                resuscitationRate * quiescent -
                (outflowRate + activeMortalityRate) * active -
                phi_main * phiCell * active * virus_main

    #Quiescent
    dstate[2] = dormancyRate * active -
                resuscitationRate * quiescent -
                (outflowRate + qMort) * quiescent -
                viralSusceptibility * phi_main * phiCell * quiescent * virus_main

    #Nutrients
    dstate[3] = nutrientInflowRate * outflowRate -
                outflowRate * nutrients -
                (1 / yieldCoefficient) * specificGrowthRate * (active + menace)

    #Virus 1 dynamics (main host only)
    dstate[4] = beta_main * burstScale * phi_main * phiCell * active * virus_main -
                phi_main * phiVirus * (active + quiescent + infectedQ) * virus_main -
                (outflowRate + virusMortalityRate) * virus_main

    #Infected quiescent (virus 1 only)
    dstate[5] = viralSusceptibility * phi_main * phiCell * quiescent * virus_main -
                qMort * infectedQ
    #Menace
    if menace <= 0.0
        dstate[6] = 0.0  # absorbing boundary at 0
    else
        dstate[6] = specificGrowthRate * menace -
                    (outflowRate + activeMortalityRate) * menace -
                    phi_menace * phiCell * menace * virus_menace
    end

    #Virus 2 dynamics (menace host only)
    dstate[7] = beta_menace * burstScale * phi_menace * phiCell * menace * virus_menace -
                phi_menace * phiVirus * menace * virus_menace -
                (outflowRate + virusMortalityRate) * virus_menace
                        
end
