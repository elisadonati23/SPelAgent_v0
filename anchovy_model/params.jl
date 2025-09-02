# Parameters for the anchovy DEB model
# These parameters are based on the DEB theory and literature values for anchovies
# They include physiological rates, shape coefficients, and other relevant parameters
# Full description of each parameter can be found in the associated documentation ODD protocol

function create_params(
    inheritance::Bool, # if true inheritance of parameters is used
    weight::Bool, # weighted inheritance of parameters
    thinning::Bool, #if true applied to juveniles
    zoom_factor::Bool, #true when inheritance of parameters is used optionally
    No_A, # Number of adult superindividuals
    No_J, # Number of juvenile superindividuals
    No_Egg, # Number of eggs
    M_f0::Union{Float64, Vector{Float64}},# Fishing mortality (from 0 to 4 /year)
    M_f1::Union{Float64, Vector{Float64}},
    M_f2::Union{Float64, Vector{Float64}},
    M_f3::Union{Float64, Vector{Float64}},
    M_f4::Union{Float64, Vector{Float64}},
    Wv, # Water basin volume (L)
    day_of_the_year,
    Xmax::Union{Float64, Vector{Float64}}, # amount of food in joules; if total, set Wv to 1.0
    Kappa::Union{Float64, Vector{Float64}}, # allocation fraction to soma
    Temp::Union{Float64, Vector{Float64}}, # temperature in degree Celsius
    M_egg::Float64, # fraction of eggs dying before hatching
    M0y::Float64, # Natural mortality (from 0 to 4 /year)
    M1y::Float64,
    M2y::Float64,
    M3y::Float64,
    M4y::Float64
)
Tc = exp.(9800.0 ./ 293.0 .- 9800.0 ./ (Temp .+ 273.0)) # Temperature correction factor
Hp = 244.0
p_Am = 11.1372
p_M = 54.67 

M0 = M0y / 365.0
M1 = M1y / 365.0
M2 = M2y / 365.0
M3 = M3y / 365.0
M4 = M4y / 365.0

model_parameters = Dict(
    :inheritance => inheritance,
    :weight => weight,
    :thinning => thinning,
    :zoom_factor => zoom_factor,
    :No_A => Float64(No_A),
    :No_J => Float64(No_J),
    :No_Egg => Float64(No_Egg),
    :Temp => Temp,
    :Tc => Tc,
    :Tc_value => Tc[1],
    :Kappa_value => Kappa isa Vector ? Kappa[1] : Kappa,
    :M_f0 => M_f0,
    :M_f1 => M_f1,
    :M_f2 => M_f2,
    :M_f3 => M_f3,
    :M_f4 => M_f4,
    :MF0_value => M_f0 isa Vector ? M_f0[1] : M_f0,
    :MF1_value => M_f1 isa Vector ? M_f1[1] : M_f1,
    :MF2_value => M_f2 isa Vector ? M_f2[1] : M_f2,
    :MF3_value => M_f3 isa Vector ? M_f3[1] : M_f3,
    :MF4_value => M_f4 isa Vector ? M_f4[1] : M_f4,
    :L50 => fit_selectivity([3.9, 12.0, 16.0, 18.0, 20.5], [0.05, 0.25, 0.75, 1.43, 1.48])[:L50],
    :slope_sel => fit_selectivity([3.9, 12.0, 16.0, 18.0, 20.5], [0.05, 0.25, 0.75, 1.43, 1.48])[:slope],
    :Wv => Wv,
    :day_of_the_year => day_of_the_year,
    :year => 1.0,
    :f => 1.0,
    :Xmax => Xmax,
    :Xmax_value => Xmax isa Vector ? Xmax[1] : Xmax,
    :M_egg => M_egg,
    :M0 => M0, #daily natural mortalities
    :M1 => M1,
    :M2 => M2,
    :M3 => M3,
    :M4 => M4,
    :death_threshold => Float64(calculate_dead_threshold_binomial(7, 365 .* [M0, M1, M2, M3, M4], 365)),
    :r_food => 0.5,
    :DEB_timing => 1.0,
    :sim_timing => 1,
    :mean_batch_eggs => 0.0,
    :mean_spawning_events => 0.0,
    :fished => 0.0,
    :fishedW => 0.0,
    :repro_start => 90.0,
    :repro_end => 270.0,
    :peak1_anchovy => 180,
    :peak2_anchovy => missing,
    :fecundity => 450.0,
    :total_repro_anchovy => 15,
    :std_dev => 60,
    :repro_period => vcat(90.0:270.0),
    :spawn_period => days_between_dates(90.0, 270.0),
    :Sex_ratio => 0.5,
    :p_Am => p_Am,
    :v_rate => 0.02,
    :Kappa => Kappa,
    :KappaX => 0.8,
    :KappaR => 0.95,
    :kap_G => 0.824129,
    :Fm => 6.5,
    :del_M => 0.1656,
    :k_J => 0.002,
    :s_M => 17.3829,
    :p_M => 54.67,
    :Eg => 5077.00,
    :d_V => 0.2,
    :mu_V => 500000.0,
    :mu_E => 550000.0,
    :w_V => 23.9,
    :w_E => 23.9,
    :w => 5.0,
    :Hb => 0.0001223,
    :Hj => 0.6741,
    :Hp => Hp,
    :Lb =>  0.0133472,
    :Lj =>  0.232013,
    :Lp => 1.50,
    :Ab => 6.0,
    :Ap => 292.0,
    :Am => 1825.0,
    :E0 =>  0.0137527,
    :ep_min => 0.30,
    :E0_min => 0.004,
    :E0_max => 0.0137527 ,
    :W0 => 2.98e-6,
    :L0 => 0.001,
    :Ta => 9800.0,
    :Tr => 293.0,
    :Xall => Xmax isa Vector ? Xmax[1] : Xmax,

    # Initialize output variables
    :year => 1.0,
    :dead_eggmass => 0,
    :mean_batch_eggs => 0.0,
    :mean_spawning_events => 0.0,
    
    #output
    :TotB => 0.0,
    :JuvB => 0.0,
    :AdB => 0.0,
    :meanAdWw => 0.0,
    :sdAdWw => 0.0,
    :meanFAdWw => 0.0,
    :sdFAdWw => 0.0,
    :meanJuvWw => 0.0,
    :sdJuvWw => 0.0,
    :meanAdL => 0.0,
    :sdAdL => 0.0,
    :meanJuvL => 0.0,
    :sdJuvL => 0.0,
    :mean_tpuberty => 0.0,
    :sd_tpuberty => 0.0,
    :mean_Lw_puberty => 0.0,
    :sd_Lw_puberty => 0.0,
    :mean_Ww_puberty => 0.0,
    :sd_Ww_puberty => 0.0,
    :mean_Hjuve => 0.0,
    :sd_Hjuve => 0.0,

    # Population structure
    :Nsuperind => No_A + No_J + No_Egg,
    :Neggs => No_Egg,
    :Njuveniles => No_J,
    :Nadults => No_A,
    :Nadult0 => 0.0,
    :Nadult1 => 0.0,
    :Nadult2 => 0.0,
    :Nadult3 => 0.0,
    :Nadult4more => 0.0,
    :mLw0 => 0.0,
    :mLw1 => 0.0,
    :mLw2 => 0.0,
    :mLw3 => 0.0,
    :mLw4more => 0.0,
    :mNind_in_Ad0 => 0.0,
    :mNind_in_Ad1 => 0.0,
    :mNind_in_Ad2 => 0.0,
    :mNind_in_Ad3 => 0.0,
    :mNind_in_Ad4more => 0.0,
    :mNind_in_Eggs => 0.0,
    :mNind_in_Juv => 0.0,

   # natural mortality
        #adults
    :deadA_nat => 0.0,
    :deadA_nat0 => 0.0,
    :deadA_nat1 => 0.0,
    :deadA_nat2 => 0.0,
    :deadA_nat3 => 0.0,
    :deadA_nat4more => 0.0,
    :natA_biom => 0.0,
    :natA_biom0 => 0.0,
    :natA_biom1 => 0.0,
    :natA_biom2 => 0.0,
    :natA_biom3 => 0.0,
    :natA_biom4more => 0.0,
    
            #juvenile
    :deadJ_nat => 0.0,
    :deadJ_nat0 => 0.0,
    :deadJ_nat1 => 0.0,
    :natJ_biom => 0.0,
    :natJ_biom0 => 0.0,
    :natJ_biom1 => 0.0,
    
        # starving mortality
            # adult
    :deadA_starved => 0.0,
    :deadA_starved0 => 0.0,
    :deadA_starved1 => 0.0,
    :deadA_starved2 => 0.0,
    :deadA_starved3 => 0.0,
    :deadA_starved4more => 0.0,
    :starvedA_biom => 0.0,
    :starvedA_biom0 => 0.0,
    :starvedA_biom1 => 0.0,
    :starvedA_biom2 => 0.0,
    :starvedA_biom3 => 0.0,
    :starvedA_biom4more => 0.0,
            # juvenile
    :starvedJ_biom => 0.0,
    :starvedJ_biom0 => 0.0,
    :starvedJ_biom1 => 0.0,
    :deadJ_starved => 0.0,
    :deadJ_starved0 => 0.0,
    :deadJ_starved1 => 0.0,
    
    
        # fishing mortality
    :fished => 0.0,
    :fishedW => 0.0,
    :fished0 => 0.0,
    :fished1 => 0.0,
    :fished2 => 0.0,
    :fished3 => 0.0,
    :fished4more => 0.0,
    :fished0_biom => 0.0,
    :fished1_biom => 0.0,
    :fished2_biom => 0.0,
    :fished3_biom => 0.0,
    :fished4more_biom => 0.0
)

    return model_parameters
         
end