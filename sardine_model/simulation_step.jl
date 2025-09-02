###############
#   Wraps     #
###############
function parallel_eggmass_step!(Sardine, model)
    eggDEB!(Sardine, model)
    eggaging!(Sardine, model)
    egghatch!(Sardine, model) # egghatch non comporta più un generate_fx() con i superindividui quindi può andare in paralelo
end

function parallel_juvenile_step!(Sardine, model)
    juvedie!(Sardine, model)
    juveDEB!(Sardine, model)
    juvemature!(Sardine,model)
    juveaging!(Sardine, model)
end

function parallel_adult_step!(Sardine, model)
    adultdie!(Sardine, model)
    adultDEB!(Sardine, model)
    adultaging!(Sardine, model)
end


function parallel_sardine_step!(Sardine, model)
    if Sardine.type == :eggmass
        parallel_eggmass_step!(Sardine, model)  # DEB + aging + hatch
    elseif Sardine.type == :juvenile
        parallel_juvenile_step!(Sardine, model)  # Die + DEB + mature + aging
    elseif Sardine.type == :adult
        parallel_adult_step!(Sardine, model)    # Die + DEB + aging 
    end
end

###################
## ENVIRONMENT   ##
###################

function evolve_environment!(model)
    # Day counter
    if model.day_of_the_year == 365.0
        model.day_of_the_year = 1.0
        model.year += 1.0
        model.fishedW = 0.0
        model.fished0 = 0.0
        model.fished1 = 0.0
        model.fished2 = 0.0
        model.fished3 = 0.0
        model.fished4more = 0.0
    else
        model.day_of_the_year += 1.0
    end


    # Increase simulation timing
    model.sim_timing += 1

    # Update time-dependent parameters
    #update_alpha!(model, model.alpha)
    update_Tc!(model, model.Tc)
    update_Kappa!(model, model.Kappa)
    update_Xmax!(model, model.Xmax)
    update_MF0!(model, model.M_f0)
    update_MF1!(model, model.M_f1)
    update_MF2!(model, model.M_f2)
    update_MF3!(model, model.M_f3)
    update_MF4!(model, model.M_f4)

    max_assimilation = calculate_max_assimilation(model)
    
    if ismissing(max_assimilation) || max_assimilation == 0.0 || isnan(max_assimilation)
        f = 0.0
    else
        # Ratio between available food and what is consumed based on size and Tc
        f = (model.Xmax_value * model.Wv * model.KappaX) / max_assimilation
    end

    ## Ensure that f is bounded between 0 and 1
    model.f = max(0, min(f, 1.0))

    adults_juve = filter(a -> a.type == :adult || a.type == :juvenile, collect(values(allagents(model))))

    ## If there are no adults or juveniles, set f to 0.8.
    ## This prevents numerical instability when there are no agents in the model that feed exogenously.
    if isempty(adults_juve)
        model.f = 0.8 
    end


    return
end

function update_outputs!(model)
    # derive population structure
    agents = collect(values(allagents(model)))
    eggsmass = filter(a -> a.type == :eggmass, agents)
    juveniles = filter(a -> a.type == :juvenile, agents)
    adults = filter(a -> a.type == :adult, agents)
    model.Nsuperind = length(agents)
    model.Neggs = isempty(eggsmass) ? 0.0 : length(eggsmass)
    model.Njuveniles = isempty(juveniles) ? 0.0 : length(juveniles)
    Ad0 = filter(a -> a.type == :adult && a.Age_years < 1.0, adults)
    Ad1 = filter(a -> a.type == :adult && a.Age_years >= 1.0 && a.Age_years < 2.0, adults)
    Ad2 = filter(a -> a.type == :adult && a.Age_years >= 2.0 && a.Age_years < 3.0, adults)
    Ad3 = filter(a -> a.type == :adult && a.Age_years >= 3.0 && a.Age_years < 4.0, adults)
    Ad4more = filter(a -> a.type == :adult && a.Age_years >= 4.0, adults)
    model.Nadult0 = isempty(Ad0) ? 0.0 : length(Ad0)
    model.Nadult1 = isempty(Ad1) ? 0.0 : length(Ad1)
    model.Nadult2 = isempty(Ad2) ? 0.0 : length(Ad2)
    model.Nadult3 = isempty(Ad3) ? 0.0 : length(Ad3)
    model.Nadult4more = isempty(Ad4more) ? 0.0 : length(Ad4more)
    model.mLw0 = isempty(Ad0) ? 0.0 : mean([a.Lw for a in Ad0])
    model.mLw1 = isempty(Ad1) ? 0.0 : mean([a.Lw for a in Ad1])
    model.mLw2 = isempty(Ad2) ? 0.0 : mean([a.Lw for a in Ad2])
    model.mLw3 = isempty(Ad3) ? 0.0 : mean([a.Lw for a in Ad3])
    model.mLw4more = isempty(Ad4more) ? 0.0 : mean([a.Lw for a in Ad4more])
    model.mNind_in_Ad0 = isempty(Ad1) ? 0.0 : mean([a.Nind for a in Ad0])
    model.mNind_in_Ad1 = isempty(Ad1) ? 0.0 : mean([a.Nind for a in Ad1])
    model.mNind_in_Ad2 = isempty(Ad2) ? 0.0 : mean([a.Nind for a in Ad2])
    model.mNind_in_Ad3 = isempty(Ad3) ? 0.0 : mean([a.Nind for a in Ad3])
    model.mNind_in_Ad4more = isempty(Ad4more) ? 0.0 : mean([a.Nind for a in Ad4more])
    model.mNind_in_Eggs = isempty(eggsmass) ? 0.0 : mean([a.Nind for a in eggsmass])
    model.mNind_in_Juv = isempty(juveniles) ? 0.0 : mean([a.Nind for a in juveniles])
    
    

    adults_juv = filter(a -> a.type == :adult || a.type == :juvenile, agents)
    if !isempty(adults_juv)
        # B plot: take into account Nind
        model.TotB = calculate_sum_prop(model,:Ww, Nind = true)
        model.JuvB = calculate_sum_prop(model,:Ww, type = :juvenile, Nind = true)
        model.AdB = calculate_sum_prop(model,:Ww, type = :adult, Nind = true)


        # Mean weight (Ww) plot
        model.meanAdWw = calculate_mean_prop(model,:Ww, type = [:adult])
        model.sdAdWw =  calculate_sd_prop(model,:Ww, type = [:adult])
        model.meanJuvWw = calculate_mean_prop(model,:Ww, type = [:juvenile])
        model.sdJuvWw = calculate_sd_prop(model,:Ww, type = [:juvenile])
        

        # Mean length (Lw) plot
        model.meanAdL = calculate_mean_prop(model, :Lw, type = [:adult])
        model.sdAdL = calculate_sd_prop(model, :Lw, type = [:adult])
        model.meanJuvL = calculate_mean_prop(model, :Lw, type = [:juvenile])
        model.sdJuvL = calculate_sd_prop(model, :Lw, type = [:juvenile])

        # Mean time to puberty plot
        model.mean_tpuberty = calculate_mean_prop(model, :t_puberty, type = [:adult])
        model.sd_tpuberty = calculate_sd_prop(model, :t_puberty, type = [:adult])

        # Mean juvenile maturation energy plot
        model.mean_Hjuve = calculate_mean_prop(model, :H, type = [:juvenile])
        model.sd_Hjuve = calculate_sd_prop(model, :H, type = [:juvenile])

        # mean scattered parameters
        model.meanHp = calculate_mean_prop(model, :Hp_i,  type = [:adult, :juvenile])
        model.sdHp = calculate_sd_prop(model, :Hp_i)#default type adult and juve
    end

    return
end

function reset_variables(model)
        # natural mortality
        #adults
        model.deadA_nat = 0.0
        model.deadA_nat1 = 0.0
        model.deadA_nat2 = 0.0
        model.deadA_nat3 = 0.0
        model.deadA_nat4more = 0.0
        model.natA_biom = 0.0
        model.natA_biom0 = 0.0
        model.natA_biom1 = 0.0
        model.natA_biom2 = 0.0
        model.natA_biom3 = 0.0
        model.natA_biom4more = 0.0
    
            #juvenile
        model.deadJ_nat = 0.0
        model.deadJ_nat0 = 0.0
        model.deadJ_nat1 = 0.0
        model.natJ_biom = 0.0
        model.natJ_biom0 = 0.0
        model.natJ_biom1 = 0.0
    
        # starving mortality
            # adult
        model.deadA_starved = 0.0
        model.deadA_starved0 = 0.0
        model.deadA_starved1 = 0.0
        model.deadA_starved2 = 0.0
        model.deadA_starved3 = 0.0
        model.deadA_starved4more = 0.0
        model.starvedA_biom = 0.0
        model.starvedA_biom0 = 0.0
        model.starvedA_biom1 = 0.0
        model.starvedA_biom2 = 0.0
        model.starvedA_biom3 = 0.0
        model.starvedA_biom4more = 0.0
            # juvenile
        model.starvedJ_biom = 0.0
        model.starvedJ_biom0 = 0.0
        model.starvedJ_biom1 = 0.0
        model.deadJ_starved = 0.0
        model.deadJ_starved0 = 0.0
        model.deadJ_starved1 = 0.0
    return
end

#################
#     Scheduler #
#################

mutable struct scheduler_Adults end

function (sEA::scheduler_Adults)(model::ABM)
    ids = [agent.id for agent in values(allagents(model))] 
    ids = filter!(id -> hasid(model, id) && (model[id].type == :adult), ids)
    return ids
end

sEA = scheduler_Adults()

function complex_step!(model)

    reset_variables(model)

    remove_all!(model, is_dead)

    # Parallel processing for Sardine agents
    Threads.@threads for sardine in collect(values(allagents(model)))
        parallel_sardine_step!(sardine, model)
    end 

    # Remove all dead agents
    remove_all!(model, is_dead)

    # Get IDs of adult agents
    sEA_ids = sEA(model)
    adult_ids = filter!(id -> hasid(model, id) && model[id].type == :adult, copy(sEA_ids))

    # Handle spawning for adult agents
    for sardine in adult_ids
        adultspawn!(model[sardine], model) 
    end

    # Filter spawners for creating new EggMass agents
    spawners = filter!(id -> hasid(model, id) && model[id].reproduction == :spawner, copy(sEA_ids))

    remove_all!(model, is_dead)

    if !isempty(spawners)
        # Create new born daily superindividuals
        prop_values = [getfield(model[agent], :superind_Neggs) for agent in spawners]
        mean_Egg_energy = mean([getfield(model[agent], :maternal_EggEn) for agent in spawners])
        
        if !model.inheritance 
            # NO INHERITANCE
            Hp_value = model.Hp * exp(rand(abmrng(model),Normal(0.0, 0.1)))
            pM_value = model.p_M * exp(rand(abmrng(model),Normal(0.0, 0.1)))
            pAm_value = model.p_Am
            K_value = model.Kappa_value
            zoom_value = 1.0
        else
            # INHERITANCE

            if model.zoom_factor
                zoom_i_values = [getfield(model[agent], :zoom_i) for agent in spawners]
                pM_values = [getfield(model[agent], :pM_i) for agent in spawners]
                    if model.weight
                        weights = [getfield(model[agent], :superind_Neggs) for agent in spawners]
                        zoom_mean = sum(zoom_i_values .* weights) / sum(weights)
                    else
                        zoom_mean = mean(zoom_i_values)
                    end
                zoom_mean = mean(zoom_i_values)
                zoom_value = rand(abmrng(model),Normal(zoom_mean, sqrt((0.066^2)/2)))
                Hp_value = model.Hp * zoom_value
                pAm_value = model.p_Am * zoom_value
                K_value = model.Kappa_value 
                pM_value = rand(abmrng(model), Normal(mean(pM_values), sqrt((15^2)/2)))

            else
                Hp_values = [getfield(model[agent], :Hp_i) for agent in spawners]
                pM_values = [getfield(model[agent], :pM_i) for agent in spawners]
                pAm_values = [getfield(model[agent], :pAm_i) for agent in spawners]

                    if !model.weight
                        # calculate arithmetic mean
                        mean_Hp = mean(Hp_values)
                        mean_pM = mean(pM_values)
                        mean_pAm_i = mean(pAm_values)
                    else
                        # calculate weighted means if considering weights as the number of eggs spawned by each spawner
                        weights = [getfield(model[agent], :superind_Neggs) for agent in spawners]
                        mean_Hp = sum(Hp_values .* weights) / sum(weights)
                        mean_pM = sum(pM_values .* weights) / sum(weights)
                        mean_pAm_i = sum(pAm_values .* weights) / sum(weights)
                    end

                ## quantitative genetics: infinite loci (Reed 2011; Dunlop 2007) 
                Hp_value = rand(abmrng(model),Normal(mean_Hp, sqrt((500^2)/2))) # the initial variance 4, std 2: we assume half of initial variance = 2, std = sqrt(2.0)
                pAm_value = rand(abmrng(model),Normal(mean_pAm_i, sqrt((15^2)/2)))
                pM_value = rand(abmrng(model),Normal(mean_pM, sqrt((15^2)/2)))
                K_value = model.Kappa_value
                zoom_value = 1.0
            end
        end

        tot_Neggs = sum(prop_values)
        Neggmass = 1
        percapita = tot_Neggs
        
        generate_EggMass(Neggmass, model, percapita, mean_Egg_energy, mean_Egg_energy, Hp_value, pM_value, pAm_value, zoom_value)
    end

    remove_all!(model, is_dead)

    # Update model outputs
    update_outputs!(model)

    # Evolve the environment
    evolve_environment!(model)

end