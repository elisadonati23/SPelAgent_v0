
#schedulers
include("dependencies.jl")

#set seed for reproducibility
Random.seed!(1234)

include("agents.jl")
include("params.jl")
include("fx.jl")
include("generate.jl")
include("initialize.jl")
include("agent_step!.jl")
include("simulation_step.jl")
include("timeseries.jl")


# -------------------------------------------
# set up
# -------------------------------------------
# NB: Wv is the water basin volume: set to 1 since timeseries were calculated as total biomass (mgC) in the basin
# Wv allows to use concentrations (mgC/m3) instead of total biomass (mgC) by setting the volume of water basin
# to the actual volume of the basin

 model = 
     model_initialize_parallel(inheritance = false, weight = false, thinning = true, zoom_factor =  false,
     No_A = 0.0, No_J = 0.0, No_Egg = 100.0,
     M_f0 = 0.0, M_f1 = 0.0, M_f2 = 0.0, M_f3 = 0.0, M_f4 = 0.0,
     Wv = 1.0, day_of_the_year = 1.0, Xmax = clima_X_tot * 0.005, Kappa = 0.883, Temp = clima_temp,
     M_egg = 0.9998, M0y = 1.08, M1y = 0.86, M2y = 0.69, M3y = 0.62, M4y = 0.48)
 
 name_sim = "test"
 length_sim = 365*3


adata = [:type, :reproduction, :Nind, :Nind0, :Age, :Age_years, :L, :H, :maternal_EggEn,
    :superind_Neggs, :En, :Dead, :death_type, :f_i, :t_puberty, :Lw, :Ww, :Wg, :R,
    :Scaled_En, :s_M_i, :pA, :pAm_i, :Lb_i, :Lj_i, :metamorph, :Hp_i, :pM_i, :K_i,
    :Lm_i, :Em_i, :g_i, :k_M_i, :zoom_i, :CI, :GSI, :dryGSI, :spawned]


mdata = [:day_of_the_year, :year, :Nsuperind, 
:Tc_value, :Xmax_value, :TotB,:JuvB,:AdB, :f, 
:deadJ_nat, :starvedJ_biom,:starvedA_biom,:natJ_biom, :natA_biom,
:deadJ_starved, :deadA_nat, :deadA_starved,
:Nadult0, :Nadult1, :Nadult2, :Nadult3, :Nadult4more,:Neggs, :Njuveniles,
:fished, :fishedW, :fished0, :fished1, :fished2, :fished3, :fished4more,
:meanJuvL, :sdJuvL, :meanAdL, :sdAdL, :mean_tpuberty, :sd_tpuberty, 
:meanJuvWw, :sdJuvWw, :meanAdWw, :sdAdWw, :mean_Hjuve, :sd_Hjuve, 
:deadA_nat0, :deadA_nat1, :deadA_nat2, :deadA_nat3, :deadA_nat4more, 
:natJ_biom0, :natJ_biom1, :deadJ_nat0, :deadJ_nat1, 
:natA_biom0, :natA_biom1, :natA_biom2, :natA_biom3, :natA_biom4more, 
:starvedA_biom0, :starvedA_biom1, :starvedA_biom2, :starvedA_biom3, :starvedA_biom4more, 
:fished0_biom, :fished1_biom, :fished2_biom, :fished3_biom, :fished4more_biom]

# ------------------------------------------------
# RUN SIMULATION
# ------------------------------------------------

    df_agent = init_agent_dataframe(model, adata)
    df_model = init_model_dataframe(model, mdata)
    df_agent = run!(model, length_sim; adata, mdata) #10 y left of spin up + forcings 

    results = []
    push!(results, df_agent)

    p1 = diagnostic_plots_pt1(results[1][1], results[1][2], model)
    p2 = diagnostic_plots_pt2(results[1][2], model)
