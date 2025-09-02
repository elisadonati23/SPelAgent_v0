
  #####################
  #      EGGMASS      #
  #####################
function eggaging!(Sardine, model)
    if !Sardine.Dead
    Sardine.Age += 1.0
    Sardine.Age_years = floor(Sardine.Age / 365.0)
    end
    return
end

function eggDEB!(Sardine, model)
    if !Sardine.Dead
        # Sardine Volume
        V = Sardine.L^3.0

        ## Initialize the variation in the state variables
        deltaV = 0.0
        deltaEggEn = 0.0
        deltaH = 0.0
        
        ## Energy fluxes

        pS = (Sardine.pM_i * model.Tc_value) * V  
        # Mobilized energy
        pC = ((Sardine.maternal_EggEn / V) * (model.Eg * (model.v_rate * model.Tc_value) * (V ^ (2/3)) + pS)/(Sardine.K_i * (Sardine.maternal_EggEn / V) + model.Eg))


        #Maturity maintenance
        pJ = model.k_J * Sardine.H * model.Tc_value
        
     
        deltaEggEn = 0.0 - pC #
        
        if ((Sardine.K_i * pC) < pS)
            model.dead_eggmass += 1.0
            Sardine.Dead = true
            return
        end
        

        deltaH =  (( 1.0 - Sardine.K_i) * pC - pJ)
        if (deltaH < 0.0 )
            deltaH = 0.0
        end
    
        deltaV = ((Sardine.K_i * pC - pS) / model.Eg)
        if (deltaV < 0.0)
            deltaV = 0.0
        end
    
        Sardine.En = Sardine.En + deltaEggEn
        Sardine.maternal_EggEn = Sardine.maternal_EggEn + deltaEggEn
        Sardine.H = Sardine.H + deltaH 
        Sardine.L = (V + deltaV)^(1/3) 
    end
    return
end

function egghatch!(Sardine, model)
    if !Sardine.Dead && (Sardine.H >= model.Hb)
        Sardine.type = :juvenile
        Sardine.f_i = model.f 
        Sardine.Lw = (Sardine.L / model.del_M)
        Sardine.Lb_i = Sardine.L
        Sardine.Age = model.Ap * (Sardine.Lw * model.del_M) / model.Lp
        Sardine.Age_years = floor(Sardine.Age / 365.0)
        Sardine.Nind = Float64(ceil((1 - model.M_egg) * Float64((Sardine.Nind))))
        Sardine.Nind0 = Sardine.Nind
        
        Sardine.s_M_i = if model.Hb >= Sardine.H
            1.0
        elseif model.Hb < Sardine.H < model.Hj
            Sardine.Lw * model.del_M / Sardine.Lb_i
        else
            model.s_M
        end

        
        Sardine.pA = Sardine.f_i * Sardine.pAm_i * model.Tc_value* Sardine.s_M_i * ((Sardine.Lw * model.del_M)^2.0)
        Sardine.Ww = (model.w * (model.d_V * ((Sardine.Lw * model.del_M) ^ 3.0) + model.w_E / model.mu_E *(Sardine.En + 0.0))) #R
        Sardine.Scaled_En = Sardine.En / ( Sardine.Em_i * ((Sardine.Lw * model.del_M)^3.0))
        Sardine.t_puberty = Sardine.Age
        model.dead_eggmass += 1.0                                              
        return
    end
    return
end


    #####################
    #      JUVENILE 
    #####################
    

function juvedie!(Sardine, model)

    # Initialize deaths
    natural_deaths = 0.0
    total_deaths = 0.0
    fishing_deaths = 0.0
  
    ee = Sardine.Scaled_En > 1.0 ? 1.0 : Sardine.Scaled_En
    ee = ee < 0.0 ? 0.0 : ee

    Lm_i = Sardine.K_i .* Sardine.pAm_i .* Sardine.s_M_i ./ Sardine.pM_i 
    g_i = model.Eg / (Sardine.K_i * Sardine.Em_i) 

    if model.thinning
        if  ee >= Sardine.L / Lm_i / Sardine.s_M_i
            r = model.Tc_value * Sardine.s_M_i * model.v_rate * ( ee / Sardine.L - 1 / (Lm_i * Sardine.s_M_i)) / ( ee + g_i)
        else
            r = model.Tc_value * Sardine.s_M_i * model.v_rate * ( ee / Sardine.L - 1 / (Lm_i * Sardine.s_M_i)) / ( ee + model.kap_G * g_i)
        end
        thin = 2/3 * r
        thin = thin / 5.0 # adjust thinning rate to be more realistic

        #if negative not applied                         
        thin < 0.0 ? thin = 0.0 : thin = thin
    else
            ## If thinning is not applied, set thin to 0.0
        thin = 0.0
    end
                

                                
    if !Sardine.Dead

        # 1st case: sardine too small to be fished
        if Sardine.Lw < 10.0 || model.MF0_value == 0.0
            # only natural mortality
            natural_deaths = Float64(rand(abmrng(model),Binomial(Int64(Sardine.Nind), 1 - exp(-(model.M0 + thin)))))
            Sardine.Nind -= natural_deaths
            model.deadJ_nat += natural_deaths
            model.natJ_biom += natural_deaths * Sardine.Ww

            # keep track of the age
            if floor(Sardine.Age / 365.0) == 0.0
                model.deadJ_nat0 += natural_deaths
                model.natJ_biom0 += natural_deaths * Sardine.Ww
            elseif floor(Sardine.Age / 365.0) == 1.0
                model.deadJ_nat1 += natural_deaths
                model.natJ_biom1 += natural_deaths * Sardine.Ww
            end
        end
                                
        # juveniles that are big enough to be fished
        if Sardine.Lw > 10.0 && !(model.MF0_value == 0.0)

            Mf = model.MF0_value / 365.0

            M = model.M0 + thin + Mf
            # natural mortality
            Mn = model.M0 + thin

            if Mf == 0.0
                total_deaths = Float64(rand(abmrng(model),Binomial(Int64(Sardine.Nind), 1 - exp(-Mn))))
                natural_deaths = total_deaths
                fishing_deaths = 0.0
            else
                total_deaths = Float64(rand(abmrng(model),Binomial(Int64(Sardine.Nind), 1 - exp(-M))))
                natural_deaths = Float64(rand(abmrng(model),Binomial(total_deaths, 1 - exp(-(Mn / M)))))
                fishing_deaths = total_deaths - natural_deaths
            end

            Sardine.Nind -= total_deaths

            # total deaths
            model.fishedW += fishing_deaths * Sardine.Ww
            model.fished += fishing_deaths
            model.deadJ_nat += natural_deaths
            model.natJ_biom += natural_deaths * Sardine.Ww
                                
            # keep track of the ages
            if floor(Sardine.Age / 365.0) == 0.0
                model.fished0 += fishing_deaths
                model.fished0_biom += fishing_deaths * Sardine.Ww
                model.deadJ_nat0 += natural_deaths
                model.natJ_biom0 += natural_deaths * Sardine.Ww
            elseif floor(Sardine.Age / 365.0) == 1.0
                model.fished1 += fishing_deaths
                model.fished1_biom += fishing_deaths * Sardine.Ww
                model.deadJ_nat1 += natural_deaths
                model.natJ_biom1 += natural_deaths * Sardine.Ww
            end
        end
    end

    # if less than a certain threshold of ind, superindividual dies
    if Sardine.Nind <= Sardine.Nind0 * model.death_threshold
        Sardine.Dead = true
        model.deadJ_nat += Sardine.Nind
        Sardine.death_type = :decline
    end
    return
end


                                
function juveDEB!(Sardine, model)

    if !Sardine.Dead

        Sardine.f_i = model.f 


        #initialize the state variables before the fluxes
        Vdyn = (Sardine.Lw * model.del_M) ^ 3.0
        Endyn = Sardine.En
        Hdyn = Sardine.H
        Rdyn = Sardine.R

        p_M_T = Sardine.pM_i * model.Tc_value 

        #initialize the variation in the state variables
        deltaV = 0.0
        deltaEn  = 0.0
        deltaH = 0.0
        deltaR = 0.0

        v_T = model.v_rate * model.Tc_value

        # Energy fluxes
        pA = (Sardine.f_i * Sardine.pAm_i* model.Tc_value * Sardine.s_M_i * (Vdyn ^ (2/3)))
        pS = p_M_T * Vdyn
        pC = ((Endyn/Vdyn) * (model.Eg * v_T * Sardine.s_M_i * (Vdyn ^ (2/3)) + pS) / (Sardine.K_i * (Endyn/ Vdyn) + model.Eg))
        pJ = model.k_J * Hdyn * model.Tc_value
        deltaEn = (pA - pC) * model.DEB_timing

        # die due to starvation
        if ((Sardine.K_i * pC) < pS)
            model.deadJ_starved += Sardine.Nind
            model.starvedJ_biom += Sardine.Nind * Sardine.Ww

            #keep track of the ages
            if (floor(Sardine.Age / 365.0 ) == 0.0)
                model.deadJ_starved0 += Sardine.Nind
                model.starvedJ_biom0 += Sardine.Nind * Sardine.Ww
            elseif (floor(Sardine.Age / 365.0 ) == 1.0)
                model.deadJ_starved1 += Sardine.Nind
                model.starvedJ_biom1 += Sardine.Nind * Sardine.Ww
            end

            Sardine.Dead = true
            Sardine.death_type = :starved
            return
        end


        deltaV = ((Sardine.K_i * pC - pS) / model.Eg) * model.DEB_timing
        if (deltaV < 0.0) 
        deltaV = 0.0
        end

        # maturing energy
        deltaH = (((1.0 - Sardine.K_i) * pC - pJ) * model.DEB_timing)
        if deltaH < 0.0
            deltaH = 0.0
        end

        # update state variables
        Sardine.En = Endyn + deltaEn
        V = Vdyn + deltaV
        Sardine.Lw = (V ^ (1/3)) / model.del_M
        Sardine.H = Hdyn + deltaH
        Sardine.R = Rdyn + deltaR
        Sardine.Ww = (model.w *(model.d_V * V + model.w_E/ model.mu_E * (Sardine.En + Sardine.R)))
        Sardine.Scaled_En = Sardine.En / (Sardine.Em_i * (( Sardine.Lw * model.del_M)^3.0))
        Sardine.L = Sardine.Lw * model.del_M
        Sardine.pA = Sardine.f_i * Sardine.pAm_i * model.Tc_value * Sardine.s_M_i * ((Sardine.Lw * model.del_M)^2.0)
  
        # adjust acceleration factor
        if !Sardine.metamorph
           if Sardine.H <= model.Hb
                Sardine.s_M_i = 1.0
            elseif model.Hb < Sardine.H < model.Hj
                Sardine.s_M_i = (Sardine.Lw * model.del_M) / Sardine.Lb_i
            elseif Sardine.H >= model.Hj
                Sardine.Lj_i = Sardine.Lw * model.del_M
                Sardine.s_M_i = Sardine.Lj_i / Sardine.Lb_i
                Sardine.metamorph = true
            end
        end

        Sardine.pA = Sardine.f_i * Sardine.pAm_i * model.Tc_value * Sardine.s_M_i * ((Sardine.Lw * model.del_M)^2.0)
        Sardine.CI = 100 * Sardine.Ww / (Sardine.Lw^3)
    end
return
end

function juvemature!(Sardine, model)
    if !Sardine.Dead && (Sardine.H >= Sardine.Hp_i)
         Sardine.type = :adult
         Sardine.R = 0.0
         Sardine.pA = Sardine.f_i * Sardine.pAm_i * model.Tc_value * Sardine.s_M_i * ((Sardine.Lw * model.del_M)^2.0) #perchè non alla 2/3?
    end
    return
end

function juveaging!(Sardine, model)
    if !Sardine.Dead
    Sardine.Age += 1.0
    Sardine.Age_years = floor(Sardine.Age / 365.0)
    Sardine.t_puberty += 1.0
    end
return
end
                                  #####################
                                  #      ADULT 
                                  #####################


  function adult_step!(Sardine, model)
    adultdie!(Sardine, model)
    adultDEB!(Sardine, model)
    adultaging!(Sardine, model)
    adultspawn!(Sardine, model)
end


function adultdie!(Sardine, model)

    # Initialize deaths
    natural_deaths = 0.0
    total_deaths = 0.0
    fishing_deaths = 0.0

    if !Sardine.Dead

        # set the new AGE DEPENDENT MORTALITIES -- If Mf is not 0, it is added to M
        if floor(Sardine.Age / 365.0) == 0.0
            Mn = model.M0 
            Mf_agent = model.MF0_value / 365.0
        elseif floor(Sardine.Age / 365.0) == 1.0
            Mn = model.M1 
            Mf_agent = model.MF1_value / 365.0
        elseif floor(Sardine.Age / 365.0) == 2.0
            Mn = model.M2
            Mf_agent = model.MF2_value / 365.0
        elseif floor(Sardine.Age / 365.0) == 3.0
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
            total_deaths = Float64(rand(abmrng(model),Binomial(Int64(Sardine.Nind), 1 - exp(-M))))

            # Calculate the number of deaths due to natural causes
            natural_deaths = Float64(rand(abmrng(model),Binomial(total_deaths, 1 - exp(-(Mn / M)))))

            # The number of deaths due to fishing is the total deaths minus the natural deaths
            fishing_deaths = total_deaths - natural_deaths
        else
            total_deaths = Float64(rand(abmrng(model),Binomial(Int64(Sardine.Nind), 1 - exp(-Mn))))
            natural_deaths = total_deaths
        end

        # Update Sardine.Nind
        Sardine.Nind -= total_deaths

        # Update total mortality events
        model.fished += fishing_deaths
        model.fishedW += fishing_deaths * Sardine.Ww
        model.deadA_nat += natural_deaths
        model.natA_biom += natural_deaths * Sardine.Ww

        # differentiate mortality with age
        if floor(Sardine.Age / 365.0) == 0.0
            model.fished0 += fishing_deaths
            model.fished0_biom += fishing_deaths * Sardine.Ww
            model.deadA_nat0 += natural_deaths
            model.natA_biom0 += natural_deaths * Sardine.Ww
        elseif floor(Sardine.Age / 365.0) == 1.0
            model.fished1 += fishing_deaths
            model.fished1_biom += fishing_deaths * Sardine.Ww
            model.deadA_nat1 += natural_deaths
            model.natA_biom1 += natural_deaths * Sardine.Ww
        elseif floor(Sardine.Age / 365.0) == 2.0
            model.fished2 += fishing_deaths
            model.fished2_biom += fishing_deaths * Sardine.Ww
            model.deadA_nat2 += natural_deaths
            model.natA_biom2 += natural_deaths * Sardine.Ww
        elseif floor(Sardine.Age / 365.0) == 3.0
            model.fished3 += fishing_deaths
            model.fished3_biom += fishing_deaths * Sardine.Ww
            model.deadA_nat3 += natural_deaths
            model.natA_biom3 += natural_deaths * Sardine.Ww
        else
            model.fished4more += fishing_deaths
            model.fished4more_biom += fishing_deaths * Sardine.Ww
            model.deadA_nat4more += natural_deaths
            model.natA_biom4more += natural_deaths * Sardine.Ww
        end
    end

    if Sardine.Nind <= Sardine.Nind0 * model.death_threshold
        Sardine.Dead = true
        model.deadA_nat += Sardine.Nind
        Sardine.death_type = :decline
    end
    return
end


function adultDEB!(Sardine, model)

if !Sardine.Dead
    Sardine.f_i = model.f
    Vdyn = (Sardine.Lw * model.del_M) ^ 3.0
    Endyn = Sardine.En
    Hdyn = Sardine.Hp_i
    Rdyn = Sardine.R

    p_M_T = Sardine.pM_i * model.Tc_value 
    
    deltaV = 0.0
    deltaEn  = 0.0
    deltaH = 0.0
    deltaR = 0.0
    
    # Energy fluxes
    
    pA = (Sardine.f_i * Sardine.pAm_i * model.Tc_value * Sardine.s_M_i * (Vdyn ^ (2/3)))
    pS = p_M_T * Vdyn
    pC = ((Endyn/Vdyn) * (model.Eg * (model.v_rate * model.Tc_value) * Sardine.s_M_i * (Vdyn ^ (2/3)) + pS) / (Sardine.K_i * (Endyn/ Vdyn) + model.Eg))
    pJ = model.k_J * Hdyn  * model.Tc_value # should not take into account the temperature?
    deltaEn = (pA - pC) * model.DEB_timing
    
    deltaV = ((Sardine.K_i * pC - pS) / model.Eg) * model.DEB_timing #pG
    if (deltaV < 0.0) 
        deltaV = 0.0
    end
    
    #starvation
    if ((Sardine.K_i * pC) < pS)
        if (Rdyn < ((pS - (Sardine.K_i * pC)) * model.DEB_timing))
            model.deadA_starved += Sardine.Nind
            model.starvedA_biom += Sardine.Nind * Sardine.Ww

            #keep track of the ages
            if (floor(Sardine.Age / 365.0 ) == 0.0)
                model.deadA_starved0 += Sardine.Nind
                model.starvedA_biom0 += Sardine.Nind * Sardine.Ww
            elseif (floor(Sardine.Age / 365.0 ) == 1.0)
                model.deadA_starved1 += Sardine.Nind
                model.starvedA_biom1 += Sardine.Nind * Sardine.Ww
            elseif (floor(Sardine.Age / 365.0 ) == 2.0)
                model.deadA_starved2 += Sardine.Nind
                model.starvedA_biom2 += Sardine.Nind * Sardine.Ww
            elseif (floor(Sardine.Age / 365.0 ) == 3.0)
                model.deadA_starved3 += Sardine.Nind
                model.starvedA_biom3 += Sardine.Nind * Sardine.Ww
            else
                model.deadA_starved4more += Sardine.Nind
                model.starvedA_biom4more += Sardine.Nind * Sardine.Ww
            end

            Sardine.Dead = true
            Sardine.death_type = :starved
            return
        else
            #take energy from repro reserve in case of starvation
            Rdyn = (Rdyn - (pS - (Sardine.K_i * pC)) * model.DEB_timing)
        end
    end

    #maturing energy
    deltaR = (((1- Sardine.K_i)* pC - pJ)* model.DEB_timing)

    if (deltaR < 0.0)
        deltaR = 0.0
    end
    
    Sardine.En = Endyn + deltaEn
    V = Vdyn + deltaV
    Sardine.Lw = (V ^ (1/3)) / model.del_M
    Sardine.H = Hdyn + deltaH
    Sardine.R = Rdyn + deltaR
    Sardine.Ww = (model.w *(model.d_V * V + model.w_E/ model.mu_E * (Sardine.En + Sardine.R)))
    Sardine.Wg = (model.w * (model.w_E / model.mu_E) * Sardine.R)
    Sardine.L = Sardine.Lw * model.del_M
    Sardine.Scaled_En= Sardine.En / (Sardine.Em_i * (( Sardine.Lw * model.del_M)^3.0))
    Sardine.pA = Sardine.f_i * Sardine.pAm_i * model.Tc_value * Sardine.s_M_i * ((Sardine.Lw * model.del_M)^2.0)
    Sardine.CI = 100 * Sardine.Ww / (Sardine.Lw^3)
    Sardine.GSI = (model.w * (model.w_E / model.mu_E) * Sardine.R) / Sardine.Ww * 100
    Sardine.dryGSI = ((model.w_E / model.mu_E) * Sardine.R) / (Sardine.Ww)* 100
    Sardine.Scaled_En = Sardine.En / (Sardine.Em_i * (( Sardine.Lw * model.del_M)^3.0))

end
return
end

function adultaging!(Sardine, model) 
    if !Sardine.Dead 
        Sardine.Age += 1.0
        Sardine.Age_years = floor(Sardine.Age / 365.0)
    end
    return
end

function adultspawn!(Sardine, model)

    Sardine.reproduction = :nonspawner
    Sardine.superind_Neggs = 0.0
    reprostart = model.repro_start + rand(abmrng(model), -14:14) 
    reproend = model.repro_end

#1st condition to reproduce not being dead
if  ((reprostart <= model.day_of_the_year <= 365.0) || (1.0 <= model.day_of_the_year <= reproend))
          
            Neggs_value_single = Float64((model.fecundity + randn(abmrng(model)) * 50) * (Sardine.Ww - Sardine.Wg))
            superind_Neggs_value = Neggs_value_single * ceil((Sardine.Nind/2.0)) 
            
            # Then determine the energy content of the eggs from maternal effects
            Sardine.maternal_EggEn = (Float64(((model.E0_max - model.E0_min) / (1.0 - model.ep_min)) * (Sardine.Scaled_En - model.ep_min)) + model.E0_min) + randn() * 0.1 * model.E0_min #NOISE
            
            spawned_en = Neggs_value_single *  Sardine.maternal_EggEn 

            if (spawned_en <= Sardine.R * (model.KappaR + (randn() * 0.01 * model.KappaR))) # NOISE
                Sardine.superind_Neggs = superind_Neggs_value
                Sardine.reproduction = :spawner
                Sardine.R = Float64(Sardine.R - spawned_en) 
                Sardine.spawned += 1.0 
            else
                Sardine.superind_Neggs = 0.0
                Sardine.reproduction = :nonspawner
            end
    end
        return
end
