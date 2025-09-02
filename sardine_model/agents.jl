@agent struct Sardine(NoSpaceAgent)
    # Basic characteristics
    type::Symbol              # :eggmass, :juvenile, :adult
    reproduction::Symbol      # :spawner, :nonspawner
    Nind::Float64             # Number of individuals in the superindividual - current
    Nind0::Float64            # Number of individuals in the superindividual - initial
    Age::Float64              # Age in days
    Age_years::Float64        # Age in years
    L::Float64                # Structural length (assumed to be close to 0 for eggs from DEB Theory)
    H::Float64                # Maturation energy
    maternal_EggEn::Float64   # Energy of the egg due to maternal effect (E0)
    superind_Neggs::Float64   # Number of eggs produced by a superindividual
    En::Float64               # Reserve energy
    Dead::Bool                # Indicates if the sardine is dead
    death_type::Symbol         # :decline, :starvation

    # Features for Juveniles and Adults
    f_i::Float64              # Individual functional response
    t_puberty::Float64        # Time to puberty
    Lw::Float64               # Length-weight relationship
    Ww::Float64               # Weight
    Wg::Float64               # Gonad weight
    R::Float64                # Reproduction energy
    Scaled_En::Float64        # Scaled energy reserve

    s_M_i::Float64            # Shape parameter
    pA::Float64               # Assimilation rate
    pAm_i::Float64            # Maximum assimilation rate (individual)
    Lb_i::Float64             # Length at birth (individual)
    Lj_i::Float64             # Length at metamorphosys (individual)
    metamorph::Bool           # Indicates if the sardine has metamorphosed -- In DEB meaning
    Hp_i::Float64             # Maturation energy at puberty (individual) -- to remove juvenile cycling
    pM_i::Float64             # Maintenance rate (individual) -- to remove juvenile cycling
    K_i::Float64              # Individual K-rule value
    #Compound parameters
    Lm_i::Float64 # Maximum length
    Em_i::Float64 # Maximum reserve density
    g_i::Float64
    k_M_i::Float64

    zoom_i::Float64 #zoom factor

    CI::Float64               # Condition Index
    GSI::Float64              # Gonadosomatic Index
    dryGSI::Float64            # Gonadosomatic Index dry weight

    # Features specific to Adults
    spawned::Float64          # Number of times the sardine has spawned
end