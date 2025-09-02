
  #####################
  #      EGGMASS      #
  #####################
function eggaging!(Anchovy, model)
    if !Anchovy.Dead
    Anchovy.Age += 1.0
    Anchovy.Age_years = floor(Anchovy.Age / 365.0)
    end
    return
end

function eggDEB!(Anchovy, model)
    if !Anchovy.Dead
        # Anchovy Volume
        V = Anchovy.L^3.0

        ## Initialize the variation in the state variables
        deltaV = 0.0
        deltaEggEn = 0.0
        deltaH = 0.0
        
        ## Energy fluxes

        pS = (Anchovy.pM_i * model.Tc_value) * V  
        # Mobilized energy
        pC = ((Anchovy.maternal_EggEn / V) * (model.Eg * (model.v_rate * model.Tc_value) * (V ^ (2/3)) + pS)/(Anchovy.K_i * (Anchovy.maternal_EggEn / V) + model.Eg))


        #Maturity maintenance
        pJ = model.k_J * Anchovy.H * model.Tc_value
        
     
        deltaEggEn = 0.0 - pC #
        
        if ((Anchovy.K_i * pC) < pS)
            model.dead_eggmass += 1.0
            Anchovy.Dead = true
            return
        end
        

        deltaH =  (( 1.0 - Anchovy.K_i) * pC - pJ)
        if (deltaH < 0.0 )
            deltaH = 0.0
        end
    
        deltaV = ((Anchovy.K_i * pC - pS) / model.Eg)
        if (deltaV < 0.0)
            deltaV = 0.0
        end
    
        Anchovy.En = Anchovy.En + deltaEggEn
        Anchovy.maternal_EggEn = Anchovy.maternal_EggEn + deltaEggEn
        Anchovy.H = Anchovy.H + deltaH 
        Anchovy.L = (V + deltaV)^(1/3) 
    end
    return
end

function egghatch!(Anchovy, model)
    if !Anchovy.Dead && (Anchovy.H >= model.Hb)
        Anchovy.type = :juvenile
        Anchovy.f_i = model.f 
        Anchovy.Lw = (Anchovy.L / model.del_M)
        Anchovy.Lb_i = Anchovy.L
        Anchovy.Age = model.Ap * (Anchovy.Lw * model.del_M) / model.Lp
        Anchovy.Age_years = floor(Anchovy.Age / 365.0)
        Anchovy.Nind = Float64(ceil((1 - model.M_egg) * Float64((Anchovy.Nind))))
        Anchovy.Nind0 = Anchovy.Nind
        
        Anchovy.s_M_i = if model.Hb >= Anchovy.H
            1.0
        elseif model.Hb < Anchovy.H < model.Hj
            Anchovy.Lw * model.del_M / Anchovy.Lb_i
        else
            model.s_M
        end

        
        Anchovy.pA = Anchovy.f_i * Anchovy.pAm_i * model.Tc_value* Anchovy.s_M_i * ((Anchovy.Lw * model.del_M)^2.0)
        Anchovy.Ww = (model.w * (model.d_V * ((Anchovy.Lw * model.del_M) ^ 3.0) + model.w_E / model.mu_E *(Anchovy.En + 0.0))) #R
        Anchovy.Scaled_En = Anchovy.En / ( Anchovy.Em_i * ((Anchovy.Lw * model.del_M)^3.0))
        Anchovy.t_puberty = Anchovy.Age
        model.dead_eggmass += 1.0                                              
        return
    end
    return
end


    #####################
    #      JUVENILE 
    #####################
    

function juvedie!(Anchovy, model)

    # Initialize deaths
    natural_deaths = 0.0
    total_deaths = 0.0
    fishing_deaths = 0.0
      
    ee = Anchovy.Scaled_En > 1.0 ? 1.0 : Anchovy.Scaled_En
    ee = ee < 0.0 ? 0.0 : ee

    Lm_i = Anchovy.K_i .* Anchovy.pAm_i .* Anchovy.s_M_i ./ Anchovy.pM_i 
    g_i = model.Eg / (Anchovy.K_i * Anchovy.Em_i) 

    if model.thinning
        if  ee >= Anchovy.L / Lm_i / Anchovy.s_M_i
            r = model.Tc_value * Anchovy.s_M_i * model.v_rate * ( ee / Anchovy.L - 1 / (Lm_i * Anchovy.s_M_i)) / ( ee + g_i)
        else
            r = model.Tc_value * Anchovy.s_M_i * model.v_rate * ( ee / Anchovy.L - 1 / (Lm_i * Anchovy.s_M_i)) / ( ee + model.kap_G * g_i)
        end
        thin = 2/3 * r
        thin = thin / 15.0 # adjust thinning rate to be more realistic

        #if negative not applied                         
        thin < 0.0 ? thin = 0.0 : thin = thin
    else
            ## If thinning is not applied, set thin to 0.0
        thin = 0.0
    end
                

                                
    if !Anchovy.Dead

        # 1st case: sardine too small to be fished
        if Anchovy.Lw < 10.0 || model.MF0_value == 0.0
            # only natural mortality
            natural_deaths = Float64(rand(abmrng(model),Binomial(Int64(Anchovy.Nind), 1 - exp(-(model.M0 + thin)))))
            Anchovy.Nind -= natural_deaths
            model.deadJ_nat += natural_deaths
            model.natJ_biom += natural_deaths * Anchovy.Ww

            # keep track of the age
            if floor(Anchovy.Age / 365.0) == 0.0
                model.deadJ_nat0 += natural_deaths
                model.natJ_biom0 += natural_deaths * Anchovy.Ww
            elseif floor(Anchovy.Age / 365.0) == 1.0
                model.deadJ_nat1 += natural_deaths
                model.natJ_biom1 += natural_deaths * Anchovy.Ww
            end
        end
                                
        # juveniles that are big enough to be fished
        if Anchovy.Lw > 10.0 && !(model.MF0_value == 0.0)

            Mf = model.MF0_value / 365.0

            M = model.M0 + thin + Mf
            # natural mortality
            Mn = model.M0 + thin

            if Mf == 0.0
                total_deaths = Float64(rand(abmrng(model),Binomial(Int64(Anchovy.Nind), 1 - exp(-Mn))))
                natural_deaths = total_deaths
                fishing_deaths = 0.0
            else
                total_deaths = Float64(rand(abmrng(model),Binomial(Int64(Anchovy.Nind), 1 - exp(-M))))
                natural_deaths = Float64(rand(abmrng(model),Binomial(total_deaths, 1 - exp(-(Mn / M)))))
                fishing_deaths = total_deaths - natural_deaths
            end

            Anchovy.Nind -= total_deaths

            # total deaths
            model.fishedW += fishing_deaths * Anchovy.Ww
            model.fished += fishing_deaths
            model.deadJ_nat += natural_deaths
            model.natJ_biom += natural_deaths * Anchovy.Ww
                                
            # keep track of the ages
            if floor(Anchovy.Age / 365.0) == 0.0
                model.fished0 += fishing_deaths
                model.fished0_biom += fishing_deaths * Anchovy.Ww
                model.deadJ_nat0 += natural_deaths
                model.natJ_biom0 += natural_deaths * Anchovy.Ww
            elseif floor(Anchovy.Age / 365.0) == 1.0
                model.fished1 += fishing_deaths
                model.fished1_biom += fishing_deaths * Anchovy.Ww
                model.deadJ_nat1 += natural_deaths
                model.natJ_biom1 += natural_deaths * Anchovy.Ww
            end
        end
    end

    # if less than a certain threshold of ind, superindividual dies
    if Anchovy.Nind <= Anchovy.Nind0 * model.death_threshold
        Anchovy.Dead = true
        model.deadJ_nat += Anchovy.Nind
        Anchovy.death_type = :decline
    end
    return
end


                                
function juveDEB!(Anchovy, model)

    if !Anchovy.Dead

        Anchovy.f_i = model.f 


        #initialize the state variables before the fluxes
        Vdyn = (Anchovy.Lw * model.del_M) ^ 3.0
        Endyn = Anchovy.En
        Hdyn = Anchovy.H
        Rdyn = Anchovy.R

        p_M_T = Anchovy.pM_i * model.Tc_value 

        #initialize the variation in the state variables
        deltaV = 0.0
        deltaEn  = 0.0
        deltaH = 0.0
        deltaR = 0.0

        v_T = model.v_rate * model.Tc_value

        # Energy fluxes
        pA = (Anchovy.f_i * Anchovy.pAm_i* model.Tc_value * Anchovy.s_M_i * (Vdyn ^ (2/3)))
        pS = p_M_T * Vdyn
        pC = ((Endyn/Vdyn) * (model.Eg * v_T * Anchovy.s_M_i * (Vdyn ^ (2/3)) + pS) / (Anchovy.K_i * (Endyn/ Vdyn) + model.Eg))
        pJ = model.k_J * Hdyn * model.Tc_value
        deltaEn = (pA - pC) * model.DEB_timing

        # die due to starvation
        if ((Anchovy.K_i * pC) < pS)
            model.deadJ_starved += Anchovy.Nind
            model.starvedJ_biom += Anchovy.Nind * Anchovy.Ww

            #keep track of the ages
            if (floor(Anchovy.Age / 365.0 ) == 0.0)
                model.deadJ_starved0 += Anchovy.Nind
                model.starvedJ_biom0 += Anchovy.Nind * Anchovy.Ww
            elseif (floor(Anchovy.Age / 365.0 ) == 1.0)
                model.deadJ_starved1 += Anchovy.Nind
                model.starvedJ_biom1 += Anchovy.Nind * Anchovy.Ww
            end

            Anchovy.Dead = true
            Anchovy.death_type = :starved
            return
        end


        deltaV = ((Anchovy.K_i * pC - pS) / model.Eg) * model.DEB_timing
        if (deltaV < 0.0) 
        deltaV = 0.0
        end

        # maturing energy
        deltaH = (((1.0 - Anchovy.K_i) * pC - pJ) * model.DEB_timing)
        if deltaH < 0.0
            deltaH = 0.0
        end

        # update state variables
        Anchovy.En = Endyn + deltaEn
        V = Vdyn + deltaV
        Anchovy.Lw = (V ^ (1/3)) / model.del_M
        Anchovy.H = Hdyn + deltaH
        Anchovy.R = Rdyn + deltaR
        Anchovy.Ww = (model.w *(model.d_V * V + model.w_E/ model.mu_E * (Anchovy.En + Anchovy.R)))
        Anchovy.Scaled_En = Anchovy.En / (Anchovy.Em_i * (( Anchovy.Lw * model.del_M)^3.0))
        Anchovy.L = Anchovy.Lw * model.del_M
        Anchovy.pA = Anchovy.f_i * Anchovy.pAm_i * model.Tc_value * Anchovy.s_M_i * ((Anchovy.Lw * model.del_M)^2.0)
  
        # adjust acceleration factor
        if !Anchovy.metamorph
           if Anchovy.H <= model.Hb
                Anchovy.s_M_i = 1.0
            elseif model.Hb < Anchovy.H < model.Hj
                Anchovy.s_M_i = (Anchovy.Lw * model.del_M) / Anchovy.Lb_i
            elseif Anchovy.H >= model.Hj
                Anchovy.Lj_i = Anchovy.Lw * model.del_M
                Anchovy.s_M_i = Anchovy.Lj_i / Anchovy.Lb_i
                Anchovy.metamorph = true
            end
        end

        Anchovy.pA = Anchovy.f_i * Anchovy.pAm_i * model.Tc_value * Anchovy.s_M_i * ((Anchovy.Lw * model.del_M)^2.0)
        Anchovy.CI = 100 * Anchovy.Ww / (Anchovy.Lw^3)
    end
return
end

function juvemature!(Anchovy, model)
    if !Anchovy.Dead && (Anchovy.H >= Anchovy.Hp_i)
         Anchovy.type = :adult
         Anchovy.R = 0.0
         Anchovy.pA = Anchovy.f_i * Anchovy.pAm_i * model.Tc_value * Anchovy.s_M_i * ((Anchovy.Lw * model.del_M)^2.0) #perchè non alla 2/3?
    end
    return
end

function juveaging!(Anchovy, model)
    if !Anchovy.Dead
    Anchovy.Age += 1.0
    Anchovy.Age_years = floor(Anchovy.Age / 365.0)
    Anchovy.t_puberty += 1.0
    end
return
end
                                  #####################
                                  #      ADULT 
                                  #####################


  function adult_step!(Anchovy, model)
    adultdie!(Anchovy, model)
    adultDEB!(Anchovy, model)
    adultaging!(Anchovy, model)
    adultspawn!(Anchovy, model)
end


function adultdie!(Anchovy, model)

    # Initialize deaths
    natural_deaths = 0.0
    total_deaths = 0.0
    fishing_deaths = 0.0

    if !Anchovy.Dead

        # set the new AGE DEPENDENT MORTALITIES -- If Mf is not 0, it is added to M
        if floor(Anchovy.Age / 365.0) == 0.0
            Mn = model.M0 
            Mf_agent = model.MF0_value / 365.0
        elseif floor(Anchovy.Age / 365.0) == 1.0
            Mn = model.M1 
            Mf_agent = model.MF1_value / 365.0
        elseif floor(Anchovy.Age / 365.0) == 2.0
            Mn = model.M2
            Mf_agent = model.MF2_value / 365.0
        elseif floor(Anchovy.Age / 365.0) == 3.0
            Mn = model.M3
            Mf_agent = model.MF3_value / 365.0
        else
            Mn = model.M4
            Mf_agent = model.MF4_value / 365.0
        end

        Mf = Mf_agent
        M = Mn + Mf

        if Mf != 0.0
            # Calculate the total number of deaths
            total_deaths = Float64(rand(abmrng(model),Binomial(Int64(Anchovy.Nind), 1 - exp(-M))))

            # Calculate the number of deaths due to natural causes
            natural_deaths = Float64(rand(abmrng(model),Binomial(total_deaths, 1 - exp(-(Mn / M)))))

            # The number of deaths due to fishing is the total deaths minus the natural deaths
            fishing_deaths = total_deaths - natural_deaths
        else
            total_deaths = Float64(rand(abmrng(model),Binomial(Int64(Anchovy.Nind), 1 - exp(-Mn))))
            natural_deaths = total_deaths
        end

        # Update Anchovy.Nind
        Anchovy.Nind -= total_deaths

        # Update total mortality events
        model.fished += fishing_deaths
        model.fishedW += fishing_deaths * Anchovy.Ww
        model.deadA_nat += natural_deaths
        model.natA_biom += natural_deaths * Anchovy.Ww

        # differentiate mortality with age
        if floor(Anchovy.Age / 365.0) == 0.0
            model.fished0 += fishing_deaths
            model.fished0_biom += fishing_deaths * Anchovy.Ww
            model.deadA_nat0 += natural_deaths
            model.natA_biom0 += natural_deaths * Anchovy.Ww
        elseif floor(Anchovy.Age / 365.0) == 1.0
            model.fished1 += fishing_deaths
            model.fished1_biom += fishing_deaths * Anchovy.Ww
            model.deadA_nat1 += natural_deaths
            model.natA_biom1 += natural_deaths * Anchovy.Ww
        elseif floor(Anchovy.Age / 365.0) == 2.0
            model.fished2 += fishing_deaths
            model.fished2_biom += fishing_deaths * Anchovy.Ww
            model.deadA_nat2 += natural_deaths
            model.natA_biom2 += natural_deaths * Anchovy.Ww
        elseif floor(Anchovy.Age / 365.0) == 3.0
            model.fished3 += fishing_deaths
            model.fished3_biom += fishing_deaths * Anchovy.Ww
            model.deadA_nat3 += natural_deaths
            model.natA_biom3 += natural_deaths * Anchovy.Ww
        else
            model.fished4more += fishing_deaths
            model.fished4more_biom += fishing_deaths * Anchovy.Ww
            model.deadA_nat4more += natural_deaths
            model.natA_biom4more += natural_deaths * Anchovy.Ww
        end
    end

    if Anchovy.Nind <= Anchovy.Nind0 * model.death_threshold
        Anchovy.Dead = true
        model.deadA_nat += Anchovy.Nind
        Anchovy.death_type = :decline
    end
    return
end


function adultDEB!(Anchovy, model)

if !Anchovy.Dead
    Anchovy.f_i = model.f
    Vdyn = (Anchovy.Lw * model.del_M) ^ 3.0
    Endyn = Anchovy.En
    Hdyn = Anchovy.Hp_i
    Rdyn = Anchovy.R

    p_M_T = Anchovy.pM_i * model.Tc_value 
    
    deltaV = 0.0
    deltaEn  = 0.0
    deltaH = 0.0
    deltaR = 0.0
    
    # Energy fluxes
    
    pA = (Anchovy.f_i * Anchovy.pAm_i * model.Tc_value * Anchovy.s_M_i * (Vdyn ^ (2/3)))
    pS = p_M_T * Vdyn
    pC = ((Endyn/Vdyn) * (model.Eg * (model.v_rate * model.Tc_value) * Anchovy.s_M_i * (Vdyn ^ (2/3)) + pS) / (Anchovy.K_i * (Endyn/ Vdyn) + model.Eg))
    pJ = model.k_J * Hdyn  * model.Tc_value # should not take into account the temperature?
    deltaEn = (pA - pC) * model.DEB_timing
    
    deltaV = ((Anchovy.K_i * pC - pS) / model.Eg) * model.DEB_timing #pG
    if (deltaV < 0.0) 
        deltaV = 0.0
    end
    
    #starvation
    if ((Anchovy.K_i * pC) < pS)
        if (Rdyn < ((pS - (Anchovy.K_i * pC)) * model.DEB_timing))
            model.deadA_starved += Anchovy.Nind
            model.starvedA_biom += Anchovy.Nind * Anchovy.Ww

            #keep track of the ages
            if (floor(Anchovy.Age / 365.0 ) == 0.0)
                model.deadA_starved0 += Anchovy.Nind
                model.starvedA_biom0 += Anchovy.Nind * Anchovy.Ww
            elseif (floor(Anchovy.Age / 365.0 ) == 1.0)
                model.deadA_starved1 += Anchovy.Nind
                model.starvedA_biom1 += Anchovy.Nind * Anchovy.Ww
            elseif (floor(Anchovy.Age / 365.0 ) == 2.0)
                model.deadA_starved2 += Anchovy.Nind
                model.starvedA_biom2 += Anchovy.Nind * Anchovy.Ww
            elseif (floor(Anchovy.Age / 365.0 ) == 3.0)
                model.deadA_starved3 += Anchovy.Nind
                model.starvedA_biom3 += Anchovy.Nind * Anchovy.Ww
            else
                model.deadA_starved4more += Anchovy.Nind
                model.starvedA_biom4more += Anchovy.Nind * Anchovy.Ww
            end

            Anchovy.Dead = true
            Anchovy.death_type = :starved
            return
        else
            #take energy from repro reserve in case of starvation
            Rdyn = (Rdyn - (pS - (Anchovy.K_i * pC)) * model.DEB_timing)
        end
    end

    #maturing energy
    deltaR = (((1- Anchovy.K_i)* pC - pJ)* model.DEB_timing)

    if (deltaR < 0.0)
        deltaR = 0.0
    end
    
    Anchovy.En = Endyn + deltaEn
    V = Vdyn + deltaV
    Anchovy.Lw = (V ^ (1/3)) / model.del_M
    Anchovy.H = Hdyn + deltaH
    Anchovy.R = Rdyn + deltaR
    Anchovy.Ww = (model.w *(model.d_V * V + model.w_E/ model.mu_E * (Anchovy.En + Anchovy.R)))
    Anchovy.Wg = (model.w * (model.w_E / model.mu_E) * Anchovy.R)
    Anchovy.L = Anchovy.Lw * model.del_M
    Anchovy.Scaled_En= Anchovy.En / (Anchovy.Em_i * (( Anchovy.Lw * model.del_M)^3.0))
    Anchovy.pA = Anchovy.f_i * Anchovy.pAm_i * model.Tc_value * Anchovy.s_M_i * ((Anchovy.Lw * model.del_M)^2.0)
    Anchovy.CI = 100 * Anchovy.Ww / (Anchovy.Lw^3)
    Anchovy.GSI = (model.w * (model.w_E / model.mu_E) * Anchovy.R) / Anchovy.Ww * 100
    Anchovy.dryGSI = ((model.w_E / model.mu_E) * Anchovy.R) / (Anchovy.Ww)* 100
    Anchovy.Scaled_En = Anchovy.En / (Anchovy.Em_i * (( Anchovy.Lw * model.del_M)^3.0))

end
return
end

function adultaging!(Anchovy, model) 
    if !Anchovy.Dead 
        Anchovy.Age += 1.0
        Anchovy.Age_years = floor(Anchovy.Age / 365.0)
    end
    return
end

function adultspawn!(Anchovy, model)

    Anchovy.reproduction = :nonspawner
    Anchovy.superind_Neggs = 0.0
    reprostart = model.repro_start + rand(abmrng(model), -14:14) 
    reproend = model.repro_end

#1st condition to reproduce not being dead
if  ((reprostart <= model.day_of_the_year) && (model.day_of_the_year <= reproend))

            Neggs_value_single = Float64((model.fecundity + randn(abmrng(model)) * 50) * (Anchovy.Ww - Anchovy.Wg))
            superind_Neggs_value = Neggs_value_single * ceil((Anchovy.Nind/2.0)) 
            
            # Then determine the energy content of the eggs from maternal effects
            Anchovy.maternal_EggEn = (Float64(((model.E0_max - model.E0_min) / (1.0 - model.ep_min)) * (Anchovy.Scaled_En - model.ep_min)) + model.E0_min) + randn() * 0.1 * model.E0_min #NOISE
            
            spawned_en = Neggs_value_single *  Anchovy.maternal_EggEn 

            if (spawned_en <= Anchovy.R * (model.KappaR + (randn() * 0.01 * model.KappaR))) # NOISE
                Anchovy.superind_Neggs = superind_Neggs_value
                Anchovy.reproduction = :spawner
                Anchovy.R = Float64(Anchovy.R - spawned_en) 
                Anchovy.spawned += 1.0 
            else
                Anchovy.superind_Neggs = 0.0
                Anchovy.reproduction = :nonspawner
            end
    end
        return
end
