# Module Generate_Agents
# Age at length and parameters are taken from AmP and DEB portal with 20°C as reference.
# All rates and ages at length change with varying parameters.

function generate_EggMass(No_Egg, model, Nind = missing, maternal_EggEn = missing, En = missing, 
    Hp = missing, pM = missing, pAm = missing, zoom = missing)
    # Initialize default agent properties for EggMass
    agent_type = :eggmass
    agent_Age = 0.0
    agent_Age_years = 0.0
    agent_L = model.L0
    agent_H = 0.0
    agent_spawned = 0.0
    agent_Dead = false
    agent_reproduction = :nonspawner
    agent_CI = 0.0
    agent_GSI = 0.0
    agent_dryGSI = 0.0
    agent_Lj_i = 0.0
    agent_metamorph = false
    agent_Wg = 0.0
    agent_death_type = :alive
    
    # Set maternal egg energy
    agent_maternal_EggEn = ismissing(maternal_EggEn) ? Float64(model.E0) : Float64(maternal_EggEn)

    agent_f_i = model.f
    agent_t_puberty = 0.0
    agent_Lw = 0.0
    agent_Ww = 0.0
    agent_R = 0.0
    agent_Scaled_En = 0.0
    agent_s_M_i = 1.0
    agent_pA = 0.0
    agent_Lb_i = 0.0
    agent_superind_Neggs = 0.0

    # Generate EggMass agents
    for _ in 1:No_Egg
        # NO INHERITANCE - thinning model
        if !model.inheritance
            agent_zoom_i = 1.0
            agent_Hp_i = model.Hp * exp(rand(abmrng(model),Normal(0.0, 0.1)))
            agent_pM_i = model.p_M * exp(rand(abmrng(model),Normal(0.0, 0.1)))
            agent_pAm_i = model.p_Am 
            agent_K_i = model.Kappa_value
        else
            if model.zoom_factor
                agent_zoom_i = ismissing(zoom) ? rand(abmrng(model),Normal(1, 0.066)) : zoom #ismissing(zoom) ? rand(abmrng(model),Normal(1, 0.066)) : zoom
                agent_Hp_i = ismissing(Hp) ? model.Hp * agent_zoom_i  : Hp
                agent_pAm_i = ismissing(pAm) ?  model.p_Am * agent_zoom_i : pAm
                agent_pM_i = ismissing(pM) ? rand(abmrng(model),Normal(model.p_M, 2*0.4))  : pM
                agent_K_i = model.Kappa_value
            else
                agent_zoom_i = 1.0
                agent_Hp_i = ismissing(Hp) ? rand(abmrng(model),Normal(model.Hp, 30))  : Hp
                agent_pAm_i = ismissing(pAm) ?  rand(abmrng(model),Normal(model.p_Am, 0.5)) : pAm
                agent_pM_i = ismissing(pM) ? rand(abmrng(model),Normal(model.p_M, 2))  : pM
                agent_K_i = model.Kappa_value
            end
        end

        
        agent_Lm_i = agent_K_i .* agent_pAm_i .* model.s_M ./ agent_pM_i
        agent_Em_i = agent_pAm_i / model.v_rate
        agent_g_i = model.Eg / (agent_K_i * agent_Em_i)
        agent_k_M_i = agent_pM_i / model.Eg
        
        agent_Nind = ismissing(Nind) ? (1e7 / 2) * (40 * 400) : Float64(floor(Nind))
        agent_Nind0 = agent_Nind

        # Set reserve energy
        agent_En = ismissing(En) ? model.E0 : Float64(En)

        # Add agent to the model
        add_agent!(
            Anchovy, model, agent_type, agent_reproduction, agent_Nind, agent_Nind0, agent_Age, agent_Age_years, agent_L, agent_H,
            agent_maternal_EggEn, agent_superind_Neggs, agent_En,  agent_Dead, agent_death_type,
            agent_f_i, agent_t_puberty, agent_Lw, agent_Ww,  agent_Wg, agent_R, agent_Scaled_En,
            agent_s_M_i, agent_pA, agent_pAm_i, agent_Lb_i, agent_Lj_i, agent_metamorph, agent_Hp_i, agent_pM_i,agent_K_i,
            agent_Lm_i, agent_Em_i, agent_g_i, agent_k_M_i, agent_zoom_i,
            agent_CI, agent_GSI, agent_dryGSI, agent_spawned
        )

    end
end

function generate_Juvenile(No_J, model, Nind = missing,  
    En = missing, Lb_i = model.Lb, Lw = missing, Ww = missing, Scaled_En = missing)

    # Initialize default agent properties for Juvenile
    agent_type = :juvenile
    agent_f_i = model.f
    agent_Lb_i = Lb_i
    agent_R = 0.0
    agent_spawned = 0.0
    agent_Dead = false
    agent_reproduction = :nonspawner
    agent_Wg = 0.0
    agent_death_type = :alive
    
    # Silenced features
    agent_maternal_EggEn = model.E0
    agent_superind_Neggs = 0.0  # EggMass

    # Generate Juvenile agents
    for _ in 1:No_J
        # NO INHERITANCE - thinning model
        if !model.inheritance
            agent_zoom_i = 1.0
            agent_Hp_i = model.Hp * exp(rand(abmrng(model),Normal(0.0, 0.1)))
            agent_pM_i = model.p_M * exp(rand(abmrng(model),Normal(0.0, 0.1)))
            agent_pAm_i = model.p_Am 
            agent_K_i = model.Kappa_value
        else
            if model.zoom_factor
                agent_zoom_i = rand(abmrng(model),Normal(1, 0.066)) #20% dalla media del valore
                agent_Hp_i =  model.Hp * agent_zoom_i
                agent_pAm_i =  model.p_Am * agent_zoom_i
                agent_pM_i = rand(abmrng(model),Normal(model.p_M, 2))
                agent_K_i = model.Kappa_value
            else
                agent_zoom_i = 1.0
                agent_Hp_i = rand(abmrng(model),Normal(model.Hp, 30))
                agent_pAm_i = rand(abmrng(model),Normal(model.p_Am, 0.5))
                agent_pM_i = rand(abmrng(model),Normal(model.p_M, 2))
                agent_K_i = model.Kappa_value
            end
        end

        agent_Lm_i = agent_K_i .* agent_pAm_i .* model.s_M ./ agent_pM_i
        agent_Em_i = agent_pAm_i / model.v_rate
        agent_g_i = model.Eg / (agent_K_i * agent_Em_i)
        agent_k_M_i = agent_pM_i / model.Eg

        agent_Nind = ismissing(Nind) ? 1e7 : Float64(floor(Nind))
        agent_Nind0 = agent_Nind

        # Set length-weight relationship
        agent_Lw = ismissing(Lw) ? clamp(round(randn() * 0.5 + 3.5, digits=2), 3.0, 4.0) : Lw
        agent_L = agent_Lw * model.del_M

        # Calculate age and time to puberty based on length
        agent_Age = model.Ap * (agent_Lw * model.del_M) / model.Lp
        agent_Age_years = agent_Age / 365.0
        agent_t_puberty = model.Ap * (agent_Lw * model.del_M) / model.Lp

        # Set weight
        agent_Ww = ismissing(Ww) ? (model.w * (model.d_V * ((agent_Lw * model.del_M) ^ 3.0))) : Ww

        # Calculate maturation energy
        agent_H = model.Hp * (agent_Lw * model.del_M) / model.Lp

        # Set reserve energy
        agent_En = ismissing(En) ? agent_f_i * agent_Em_i * ((agent_Lw * model.del_M)^3.0) : En

        # Calculate scaled energy reserve
        agent_Scaled_En = ismissing(Scaled_En) ? agent_En / (agent_Em_i * ((agent_Lw * model.del_M)^3.0)) : Scaled_En

        # Determine shape parameter

        if model.Hb >= agent_H
            agent_s_M_i = 1.0
            agent_metamorph = false
        elseif model.Hb < agent_H < model.Hj
            agent_s_M_i = agent_Lw * model.del_M / agent_Lb_i
            agent_metamorph = false
        else
            agent_Lj_i = model.Lj
            agent_s_M_i = model.Lj / agent_Lb_i
            agent_metamorph = true
        end
        
        # Calculate assimilation rate
        Tc_value = isa(model.Tc, Vector{Float64}) ? model.Tc[model.sim_timing] : model.Tc
        agent_pA = agent_f_i * model.p_Am * Tc_value * agent_s_M_i * ((agent_Lw * model.del_M)^2.0)

        # Calculate condition index and gonadosomatic index
        agent_CI = 100 * agent_Ww / (agent_Lw^3)
        agent_GSI = 0.0
        agent_dryGSI = 0.0

        # Add agent to the model
        add_agent!(
            Anchovy, model, agent_type, agent_reproduction, agent_Nind, agent_Nind0, agent_Age,agent_Age_years, agent_L, agent_H,
            agent_maternal_EggEn, agent_superind_Neggs, agent_En,  agent_Dead,agent_death_type,
            agent_f_i, agent_t_puberty, agent_Lw, agent_Ww,  agent_Wg, agent_R, agent_Scaled_En,
            agent_s_M_i, agent_pA, agent_pAm_i, agent_Lb_i, agent_Lj_i, agent_metamorph, agent_Hp_i, agent_pM_i, agent_K_i,
            agent_Lm_i, agent_Em_i, agent_g_i, agent_k_M_i, agent_zoom_i,
            agent_CI, agent_GSI, agent_dryGSI, agent_spawned
        )

    end
end

function generate_Adult(No_A, model, Nind = missing, Age = missing, t_puberty = missing, Lw = missing, Ww = missing, H = missing, R = missing, En = missing, Scaled_En = missing,  pA = missing, Lj = missing)

    # Initialize default agent properties for Adult
    agent_type = :adult
    agent_f_i = model.f
    agent_reproduction = :nonspawner
    agent_maternal_EggEn = model.E0
    agent_superind_Neggs = 0.0
    agent_Lb_i = model.Lb
    agent_spawned = 0.0
    agent_Dead = false
    agent_metamorph = true
    agent_Wg = 0.0
    agent_death_type = :alive


    # Set maturation energy
    agent_H = ismissing(H) ? model.Hp : H

    # Determine shape parameter
    agent_s_M_i = model.Lj / model.Lb

    # Generate Adult agents
    for _ in 1:No_A

        # NO INHERITANCE - thinning model
        if !model.inheritance
            agent_zoom_i = 1.0
            agent_Hp_i = model.Hp * exp(rand(abmrng(model),Normal(0.0, 0.1)))
            agent_pM_i = model.p_M * exp(rand(abmrng(model),Normal(0.0, 0.1)))
            agent_pAm_i = model.p_Am
            agent_K_i = model.Kappa_value
        else
            if model.zoom_factor
                agent_zoom_i = rand(abmrng(model),Normal(1, 0.066))
                agent_Hp_i = model.Hp * agent_zoom_i
                agent_pAm_i = model.p_Am * agent_zoom_i
                agent_pM_i = rand(abmrng(model),Normal(model.p_M, 2))
                agent_K_i = model.Kappa_value
            else
                agent_zoom_i = 1.0
                agent_Hp_i = rand(abmrng(model),Normal(model.Hp, 30))
                agent_pAm_i = rand(abmrng(model),Normal(model.p_Am, 0.5))
                agent_pM_i = rand(abmrng(model),Normal(model.p_M, 2))
                agent_K_i = model.Kappa_value
            end
        end


        agent_Lm_i = agent_K_i .* agent_pAm_i .* model.s_M ./ agent_pM_i
        agent_Em_i = agent_pAm_i / model.v_rate #maximum reserve density
        agent_g_i = model.Eg / (agent_K_i * agent_Em_i)
        agent_k_M_i = agent_pM_i / model.Eg

        # Set length-weight relationship
        agent_Lw = ismissing(Lw) ? clamp(round(randn() * 5.0 + 13.0, digits=2), 8.0, 18.0) : agent_Lw
        agent_L = agent_Lw * model.del_M

                
        agent_Age = if ismissing(Age)
            if agent_Lm_i isa Float64
                model.Am * agent_Lw * model.del_M / agent_Lm_i
            else
                model.Am * agent_Lw * model.del_M / agent_Lm_i[model.sim_timing]
            end
        else
            Age
        end
        agent_Age_years = agent_Age / 365.0

        if ismissing(Nind)
            agent_Nind = 1e7
            agent_Nind0 = calculate_Nind0(Float64(agent_Nind), Int64(floor(agent_Age)), 365 .* [model.M0, model.M1, model.M2, model.M3, model.M4], 365)
        else
            agent_Nind = Nind
            agent_Nind0 = calculate_Nind0(Float64(agent_Nind), Int64(floor(agent_Age)), 365 .* [model.M0, model.M1, model.M2, model.M3, model.M4], 365)
        end

        agent_Lj_i = if ismissing(Lj)
            model.Lj
        else
            Lj
        end

        # Set time to puberty
        agent_t_puberty = ismissing(t_puberty) ? model.Ap * (agent_Lw * model.del_M) / model.Lp : t_puberty
        
        # Set reproduction energy
        agent_R = ismissing(R) ? 0.0 : R

        # Set reserve energy
        agent_En = ismissing(En) ? agent_f_i * agent_Em_i * ((agent_Lw * model.del_M)^3.0) : En

        # Calculate weight
        agent_Ww = ismissing(Ww) ? (model.w * (model.d_V * ((agent_Lw * model.del_M) ^ 3.0) + model.w_E / model.mu_E * (agent_En + agent_R))) : Ww

        # Calculate scaled energy reserve density
        agent_Scaled_En = ismissing(Scaled_En) ? agent_En / (agent_Em_i * ((agent_Lw * model.del_M)^3.0)) : Scaled_En
        
        # Calculate assimilation rate
        Tc_value = isa(model.Tc, Vector{Float64}) ? model.Tc[model.sim_timing] : model.Tc
        agent_pA = ismissing(pA) ? agent_f_i * model.p_Am * Tc_value * agent_s_M_i * ((agent_Lw * model.del_M)^2.0) : pA

        # Calculate condition index and gonadosomatic index
        agent_CI = 100 * agent_Ww / (agent_Lw^3)
        agent_GSI = (model.w * (model.w_E / model.mu_E) * agent_R) / agent_Ww * 100  # Weight of the gonads as a percentage of total weight
        agent_dryGSI = ((model.w_E / model.mu_E) * agent_R) / agent_Ww * 100  # Weight of the gonads as a percentage of total weight

        # Add agent to the model
        add_agent!(
            Anchovy, model, agent_type, agent_reproduction, agent_Nind, agent_Nind0, agent_Age, agent_Age_years, agent_L, agent_H,
            agent_maternal_EggEn, agent_superind_Neggs, agent_En,  agent_Dead,agent_death_type,
            agent_f_i, agent_t_puberty, agent_Lw, agent_Ww, agent_Wg, agent_R, agent_Scaled_En,
            agent_s_M_i, agent_pA, agent_pAm_i, agent_Lb_i, agent_Lj_i, agent_metamorph, agent_Hp_i, agent_pM_i, agent_K_i,
            agent_Lm_i, agent_Em_i, agent_g_i, agent_k_M_i, agent_zoom_i,
            agent_CI, agent_GSI, agent_dryGSI, agent_spawned
        )
    end

    return
end

function generate_Adult_byAge(No_A, model, Age, Nind = missing)

    # Initialize default agent properties for Adult
    agent_type = :adult
    agent_f_i = model.f
    agent_reproduction = :nonspawner
    agent_maternal_EggEn = model.E0
    agent_superind_Neggs = 0.0
    agent_Lb_i = model.Lb
    agent_spawned = 0.0
    agent_Dead = false
    agent_metamorph = true
    agent_Wg = 0.0
    agent_Hp_i = model.Hp
    agent_death_type = :alive


    # Determine shape parameter
    agent_s_M_i = model.Lj / model.Lb

    # Set default values based on Age
    default_R = Age == 1.0 ? 470.0 : Age == 2.0 ? 900.0 : 500.0
    default_Lw_mean = Age == 1.0 ? 13.0 : Age == 2.0 ? 16.0 : 18.0
    default_Lw_range = Age == 1.0 ? (8.0, 15.0) : Age == 2.0 ? (12.0, 18.0) : (16.0, 22.0)
    default_Age_range = Age == 1.0 ? (365:365*2) : Age == 2.0 ? (365*2:365*3) : (365*3:365*4)

    # Generate Adult agents
    for _ in 1:No_A

        if !model.inheritance
            # NO INHERITANCE - thinning model
            agent_zoom_i = 1.0
            agent_Hp_i = model.Hp * exp(rand(abmrng(model),Normal(0.0, 0.1)))
            agent_H = agent_Hp_i
            agent_pM_i = model.p_M * exp(rand(abmrng(model),Normal(0.0, 0.1)))
            agent_pAm_i = model.p_Am
            agent_K_i = model.Kappa_value
        else
            if model.zoom_factor
                agent_zoom_i = rand(abmrng(model),Normal(1, 0.066))
                agent_Hp_i = model.Hp * agent_zoom_i
                agent_pAm_i = model.p_Am * agent_zoom_i
                agent_pM_i = rand(abmrng(model),Normal(model.p_M, 2))
                agent_K_i = model.Kappa_value
            else
                agent_zoom_i = 1.0
                agent_Hp_i = rand(abmrng(model),Normal(model.Hp, 30))
                agent_pAm_i = rand(abmrng(model),Normal(model.p_Am, 0.5))
                agent_pM_i = rand(abmrng(model),Normal(model.p_M, 2))
                agent_K_i = model.Kappa_value
            end
        end

        agent_H = agent_Hp_i

        
        agent_Lm_i = agent_K_i .* agent_pAm_i .* model.s_M ./ agent_pM_i
        agent_Em_i = agent_pAm_i / model.v_rate
        agent_g_i = model.Eg / (agent_K_i * agent_Em_i)
        agent_k_M_i = agent_pM_i / model.Eg

        # Set length-weight relationship
        agent_Lw = clamp(round(randn() * 2.0 + default_Lw_mean, digits=2), default_Lw_range...)
        agent_L = agent_Lw * model.del_M

        # Set age
        agent_Age = rand(abmrng(model),default_Age_range)
        agent_Age_years = floor(agent_Age / 365.0)

        if ismissing(Nind)
            agent_Nind = 1e7
            agent_Nind0 = calculate_Nind0(Float64(agent_Nind), Int64(floor(agent_Age)), 365 .* [model.M0, model.M1, model.M2, model.M3, model.M4], 365)
        else
            agent_Nind = Nind
            agent_Nind0 = calculate_Nind0(Float64(agent_Nind), Int64(floor(agent_Age)), 365 .* [model.M0, model.M1, model.M2, model.M3, model.M4], 365)
        end

        agent_Lj_i = model.Lj

        # Set time to puberty
        agent_t_puberty = model.Ap * (agent_Lw * model.del_M) / model.Lp

        # Set reproduction energy
        agent_R = max(rand(abmrng(model),Normal(default_R, 70.0)), 0.0)

        # Set reserve energy
        agent_En = agent_f_i * agent_Em_i * ((agent_Lw * model.del_M)^3.0)

        # Calculate weight
        agent_Ww = (model.w * (model.d_V * ((agent_Lw * model.del_M)^3.0) + model.w_E / model.mu_E * (agent_En + agent_R)))

        # Calculate scaled energy reserve
        agent_Scaled_En = agent_En / (agent_Em_i * ((agent_Lw * model.del_M)^3.0))

        # Calculate assimilation rate
        Tc_value = isa(model.Tc, Vector{Float64}) ? model.Tc[model.sim_timing] : model.Tc
        agent_pA = agent_f_i * agent_pAm_i * Tc_value * agent_s_M_i * ((agent_Lw * model.del_M)^2.0)

        # Calculate condition index and gonadosomatic index
        agent_CI = 100 * agent_Ww / (agent_Lw^3)
        agent_GSI = (model.w * (model.w_E / model.mu_E) * agent_R) / agent_Ww * 100
        agent_dryGSI = ((model.w_E / model.mu_E) * agent_R) / agent_Ww * 100  # Weight of the gonads as a percentage of total weight


        # Add agent to the model
        add_agent!(
            Anchovy, model, agent_type, agent_reproduction, agent_Nind, agent_Nind0, agent_Age,agent_Age_years, agent_L, agent_H,
            agent_maternal_EggEn, agent_superind_Neggs, agent_En,  agent_Dead, agent_death_type,
            agent_f_i, agent_t_puberty, agent_Lw, agent_Ww,  agent_Wg, agent_R, agent_Scaled_En,
            agent_s_M_i, agent_pA, agent_pAm_i, agent_Lb_i, agent_Lj_i, agent_metamorph, agent_Hp_i, agent_pM_i, agent_K_i,
            agent_Lm_i, agent_Em_i, agent_g_i, agent_k_M_i, agent_zoom_i,
            agent_CI, agent_GSI, agent_dryGSI, agent_spawned
        )
    end

    return
end



# DEPRECATED #

function generate_adult_pop(model, Lwclass = missing, Lw_biom = missing)
    # to be used with medias info on length class and biomass
    #USE ONLY TO Initialize MODEL FROM MEDIAS DATA

        # Initialize default agent properties for Adult - true for everyone
        agent_type = :adult
        agent_f_i = model.f
        agent_reproduction = :nonspawner
        agent_maternal_EggEn = model.E0
        agent_superind_Neggs = 0.0
        agent_Lb_i = model.Lb
        agent_spawned = 0.0
        agent_Dead = false
        agent_metamorph = true
        agent_H = model.Hp
        agent_s_M_i = model.Lj / model.Lb
        agent_Wg = 0.0
        agent_Hp_i = model.Hp
        agent_pM_i = model.p_M + randn() * 0.01 * model.p_M
    

    #Lw class - biom relationship - based on MEDIAS
    agent_R =  0.0
    # Set reserve energy
    agent_En_class = agent_f_i * agent_Em_i * ((Lwclass * model.del_M)^3.0)
    agent_Ww_class = (model.w * (model.d_V * ((Lwclass * model.del_M) ^ 3.0) + model.w_E / model.mu_E * (agent_En_class + agent_R))) 

    Nind = Int64(floor(Lw_biom/agent_Ww_class))

    # Define No_A based on Nind if Nind is greater than 1e7
    if Nind > 1e7
        No_A = ceil(Int, Nind / 1e7)
        agent_Nind = ceil(Int, Nind / No_A)
        agent_Nind0 = agent_Nind
    else
        No_A = 1
        agent_Nind = Nind
        agent_Nind0 = Nind # Adjust the divisor as needed to define the number of superindividuals
    end

    # Generate Adult agents
    for _ in 1:No_A
        # Set number of individuals in the superindividual
        agent_Lw = clamp(round(randn() * 0.5 + Lwclass, digits=2), Lwclass - 0.5, Lwclass+0.5)
        agent_L = agent_Lw * model.del_M
        # Set reproduction energy
        agent_En = agent_f_i * agent_Em_i * ((agent_Lw * model.del_M)^3.0)

        agent_Ww = (model.w * (model.d_V * ((agent_Lw * model.del_M) ^ 3.0) + model.w_E / model.mu_E * (agent_En + agent_R))) 

        # Calculate age based on length
        agent_Age = 
            if model.Lm isa Float64
            model.Am * agent_Lw * model.del_M / model.Lm
            else
            model.Am * agent_Lw * model.del_M / model.Lm[model.sim_timing]
            end


        
        # Calculate scaled energy reserve
        agent_Scaled_En = agent_En / (agent_Em_i * ((agent_Lw * model.del_M)^3.0))

         # depending on age, above, i put the Nind to ensure correct lifespan but I assume that they were 1e7 at age 0+
        
        agent_Lj_i = 
            model.Lj


        # Set time to puberty
        agent_t_puberty = model.Ap * (agent_Lw * model.del_M) / model.Lp
        
        # Calculate assimilation rate
        Tc_value = isa(model.Tc, Vector{Float64}) ? model.Tc[model.sim_timing] : model.Tc
        agent_pA = agent_f_i * model.p_Am * Tc_value * agent_s_M_i * ((agent_Lw * model.del_M)^2.0)

        # Calculate condition index and gonadosomatic index
        agent_CI = 100 * agent_Ww / (agent_Lw^3)
        agent_GSI = (model.w * (model.w_E / model.mu_E) * agent_R) / agent_Ww * 100  # Weight of the gonads as a percentage of total weight
        agent_dryGSI = ((model.w_E / model.mu_E) * agent_R) / agent_Ww * 100  # Weight of the gonads as a percentage of total weight

        # Add agent to the model
        add_agent!(
            Anchovy, model, agent_type, agent_reproduction, agent_Nind, agent_Nind0, agent_Age, agent_L, agent_H,
            agent_maternal_EggEn, agent_superind_Neggs, agent_En,  agent_Dead,
            agent_f_i, agent_t_puberty, agent_Lw, agent_Ww,  agent_Wg, agent_R, agent_Scaled_En,
            agent_s_M_i, agent_pA, agent_Lb_i, agent_Lj_i, agent_metamorph, agent_Hp_i, agent_pM_i,
            agent_CI, agent_GSI, agent_dryGSI, agent_spawned
        )
    end
    return
end


function generate_juvenile_pop(model, Lwclass = missing, Lw_biom = missing)
    # to be used with medias info on length class and biomass

        # Initialize default agent properties for Adult - true for everyone
        agent_type = :juvenile
        agent_f_i = model.f
        agent_reproduction = :nonspawner
        agent_maternal_EggEn = model.E0
        agent_superind_Neggs = 0.0
        agent_Lb_i = model.Lb
        agent_spawned = 0.0
        agent_Dead = false
        agent_Wg = 0.0
        agent_Hp_i = model.Hp + randn() * 0.01 * model.Hp
        agent_pM_i = model.p_M + randn() * 0.01 * model.p_M
    

    #Lw class - biom relationship - based on MEDIAS
    agent_R =  0.0
    # Set reserve energy
    agent_En_class = agent_f_i * agent_Em_i * ((Lwclass * model.del_M)^3.0)
    agent_Ww_class = (model.w * (model.d_V * ((Lwclass * model.del_M) ^ 3.0) + model.w_E / model.mu_E * (agent_En_class + agent_R))) 
    Nind = Int64(floor(Lw_biom/agent_Ww_class))

    # Define No_A based on Nind if Nind is greater than 1e7
    if Nind > 1e7
        No_J = ceil(Int, Nind / 1e7)
        agent_Nind = ceil(Int, Nind / No_J)
        agent_Nind0 = agent_Nind
    else
        No_J = 1
        agent_Nind = Nind
        agent_Nind0 = Nind # Adjust the divisor as needed to define the number of superindividuals
    end

    # Generate Adult agents
    for _ in 1:No_J
        # Set number of individuals in the superindividual
        agent_Lw = clamp(round(randn() * 0.5 + Lwclass, digits=2), Lwclass - 0.5, Lwclass+0.5)
        agent_L = agent_Lw * model.del_M
        # Set reproduction energy
        agent_En = agent_f_i * agent_Em_i * ((agent_Lw * model.del_M)^3.0)
        # Calculate scaled energy reserve
        agent_Scaled_En = agent_En / (agent_Em_i * ((agent_Lw * model.del_M)^3.0))
        agent_Ww = (model.w * (model.d_V * ((agent_Lw * model.del_M) ^ 3.0) + model.w_E / model.mu_E * (agent_En + agent_R))) 
        agent_H = model.Hp * (agent_Lw * model.del_M) / model.Lp

        if model.Hb >= agent_H
            agent_s_M_i = 1.0
            agent_metamorph = false
        elseif model.Hb < agent_H < model.Hj
            agent_s_M_i = agent_Lw * model.del_M / agent_Lb_i
            agent_metamorph = false
        else
            agent_Lj_i = model.Lj
            agent_s_M_i = model.Lj / agent_Lb_i
            agent_metamorph = true
        end

        # Calculate age based on length
        agent_Age = model.Ap * (agent_Lw * model.del_M) / model.Lp
        agent_t_puberty = model.Ap * (agent_Lw * model.del_M) / model.Lp
        
        # Calculate assimilation rate
        Tc_value = isa(model.Tc, Vector{Float64}) ? model.Tc[model.sim_timing] : model.Tc
        agent_pA = agent_f_i * model.p_Am * Tc_value * agent_s_M_i * ((agent_Lw * model.del_M)^2.0)

        # Calculate condition index and gonadosomatic index
        agent_CI = 100 * agent_Ww / (agent_Lw^3)
        agent_GSI = (model.w * (model.w_E / model.mu_E) * agent_R) / agent_Ww * 100  # Weight of the gonads as a percentage of total weight
        agent_dryGSI = ((model.w_E / model.mu_E) * agent_R) / agent_Ww * 100  # Weight of the gonads as a percentage of total weight
        # Add agent to the model
        add_agent!(
            Anchovy, model, agent_type, agent_reproduction, agent_Nind, agent_Nind0, agent_Age, agent_L, agent_H,
            agent_maternal_EggEn, agent_superind_Neggs, agent_En,  agent_Dead,
            agent_f_i, agent_t_puberty, agent_Lw, agent_Ww, agent_Wg, agent_R, agent_Scaled_En,
            agent_s_M_i, agent_pA, agent_Lb_i, agent_Lj_i, agent_metamorph, agent_Hp_i, agent_pM_i,
            agent_CI, agent_GSI,agent_dryGSI,  agent_spawned
        )
    end
    return
end
