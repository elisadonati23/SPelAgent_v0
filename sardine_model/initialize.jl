function model_initialize_parallel(;
    inheritance,
    weight,
    thinning,
    zoom_factor,
    No_A,
    No_J, 
    No_Egg, 
    M_f0,
    M_f1,
    M_f2,
    M_f3,
    M_f4, 
    Wv,
    day_of_the_year,
    Xmax,
    Kappa,
    Temp,
    M_egg,
    M0y,
    M1y,
    M2y,
    M3y,
    M4y
) 

    # Generate model parameters using the provided inputs

    properties = create_params(
        inheritance,
        weight,
        thinning,
        zoom_factor,
        No_A,
        No_J,
        No_Egg,
        M_f0,
        M_f1,
        M_f2,
        M_f3,
        M_f4, 
        Wv,
        day_of_the_year,
        Xmax,
        Kappa,
        Temp,
        M_egg,
        M0y,
        M1y,
        M2y,
        M3y,
        M4y
    )


    # Create the Agent-Based Model (ABM) for Sardines
    model = ABM(
        Sardine;
        properties = properties,
        model_step! = complex_step!
    )

    # Add agents to the model: Adults, Juveniles, and EggMass
    generate_Adult(No_A, model)
    generate_Juvenile(No_J, model)
    generate_EggMass(No_Egg, model)


    # Calculate the mean length-weight (Lw) for initialization
    mean_Lw = calculate_mean_prop(model, :Lw)

    agents = collect(values(allagents(model)))
    adults_juve = filter(a -> a.type == :adult || a.type == :juvenile, agents)
                        
    # If there are no adults or juveniles, set f to 0.8.
    # This prevents numerical instability when there are no agents in the model that feed exogenously.
    # Determine the initial value of functional response (f)
    if isempty(adults_juve)
        model.f = 0.8
    else
        # sim_timing = 1 at initialization, so directly index Xmax and Tc
        f = (model.Xmax[model.sim_timing] * model.Wv * model.KappaX) / 
            (model.p_Am * model.Tc[model.sim_timing] * model.s_M * ((mean_Lw * model.del_M) ^ 2))

        # Ensure f is within the range [0, 1]
        model.f = max(0, min(f, 1.0))
        #model.f_adult = model.f
    end

        # Assign the calculated f to all agents in the model
        for agent in agents
            #no distinction at beginning between adults and juveniles for f
            agent.f_i = model.f
        end

    return model
end
