# structure containing all parameters defining the model

mutable struct Options
    # time stepping
    tStep::Int  # time step size
    N::Int  # number of time steps
    N100::Int # timestep index of year 100
    # labour
    Lₐ::Int # asymptotic popultation in million people
    l_g::Float64  # population growth rate
    # factor productivity
    δₐ::Float64 # Decline rate of TFP per 5 years
    gₐ::Float64  # Initial growth rate for TFP per 5 years
    # carbon intensity
    δσ::Float64 #Decline rate of decarbonization per period
    gσ::Float64 #Initial growth of sigma (continuous per year)
    # carbon effects
    f₀::Float64 #  2010 forcings of non-CO2 GHG (Wm-2)
    f₁::Float64 # 2100 forcings of non-CO2 GHG (Wm-2)
    E_L0::Float64  # Carbon emissions from land 2010 (GtCO2 per year)
    δEL::Float64 # Decline rate of land emissions (per period)
    # climate model
    Φ_T::Array{Float64,2}  # temperature transfer coefficients
    b_T::Array{Float64,1}  # temperature forcings coefficient
    M_AT1750::Float64  # reference carbon resevoir
    η::Float64 #Forcings of equilibrium CO2 doubling (Wm-2)
    # CO2 model
    Φ_M::Array{Float64,2}  # carbon stock transfer coefficients
    b_M::Array{Float64,1}  # carbon forcings coefficient
    # economic model
    a₂::Float64  # linear damages coefficient
    a₃::Float64  # damages exponential coefficient
    θ₂::Float64  # mitigation effort exponent
    δₖ::Float64 # Depreciation rate on capital (per year)
    γ::Float64  # Capital elasticity in production function (Capital share)
    # cost of mititgation
    pb::Float64  # Cost of backstop 2005$ per tCO2 2010
    δpb::Float64 #Initial cost decline backstop cost per period
    ρ::Float64  # social discount rate
    μ_ubound::Float64  # Upper limit on control rate after 2150
    fosslim::Float64 #Maximum cumulative extraction fossil fuels (GtC)
    α::Float64  # Elasticity of marginal utility of consumption
end

# structure containing all initial values
struct InitialValues
    μ₀::Float64 #Initial emissions control rate for base case 2010
    S₀::Float64  # initial saving rate
    L₀::Int  # Initial population
    A₀::Float64  # Initial factor productivity
    σ₀::Float64  # Initial carbon intensity
    T₀::Array{Float64,1}  # initial temperature
    M₀::Array{Float64,1}  # initial carbon stock
    K₀::Float64  # initial capital available
end

# struture containing all exogenous variables, which will be pre-computed

mutable struct ExogenousVariables
    L::Array{Float64, 1}  # labour
    A::Array{Float64, 1}  # factor productivity
    σ::Array{Float64, 1}  # carbon intensity
    F_EX::Array{Float64, 1}  # forcings besides carbon
    E_LAND::Array{Float64, 1}  # emissions from deforestation
    θ₁::Array{Float64, 1}  # cost of mitigation
end

# structure containing the variables for optimization

mutable struct Variables
    μ::Array{VariableRef, 1}  # mitigation efforts
    S::Array{VariableRef, 1}  # savings rate
    Ω::Array{VariableRef, 1}  # damages
    Λ::Array{VariableRef, 1}  # mitigation effort costs
    Y::Array{VariableRef, 1}  # gross WP
    Q::Array{VariableRef, 1}  # net WP
    C::Array{VariableRef, 1}  # consumption
    E::Array{VariableRef, 1}  # ausstoss
    K::Array{VariableRef, 1}  # capital
    Tₐₜ::Array{VariableRef, 1}  # temperature
    Tₗₒ::Array{VariableRef, 1}  # temperature
    Mₐₜ::Array{VariableRef, 1}  # carbon stock
    Mᵤₚ::Array{VariableRef, 1}  # carbon stock
    Mₗₒ::Array{VariableRef, 1}  # carbon stock
    U::Array{VariableRef, 1}  # carbon stock
    W::VariableRef  # social welfare
end

# structure containing the results

mutable struct Results
    years::Array{Int64, 1}
    μ::Array{Float64, 1}  # mitigation efforts
    S::Array{Float64, 1}  # savings rate
    Ω::Array{Float64, 1}  # damages
    Λ::Array{Float64, 1}  # mitigation effort costs
    Y::Array{Float64, 1}  # gross WP
    Q::Array{Float64, 1}  # net WP
    C::Array{Float64, 1}  # consumption
    E::Array{Float64, 1}  # ausstoss
    K::Array{Float64, 1}  # capital
    Tₐₜ::Array{Float64, 1}  # temperature
    Tₗₒ::Array{Float64, 1}  # temperature
    Mₐₜ::Array{Float64, 1}  # carbon stock
    Mᵤₚ::Array{Float64, 1}  # carbon stock
    Mₗₒ::Array{Float64, 1}  # carbon stock
    U::Array{Float64, 1}  # carbon stock
    W::Float64  # social welfare
end

# handy simulation object, that has all interesting stuff in one place

mutable struct DiceSimulation
    scenario::String
    options::Options
    initialValues::InitialValues
    exogenousVariables::ExogenousVariables
    #variables::Variables
    results::Results
    DiceSimulation(sc, opts, iv) = new(sc, opts, iv, computeExogenousVariables(opts, iv))
    DiceSimulation() = new()
end










mutable struct VariablesSim
    Ω::Array{Float64, 1}  # damages
    Λ::Array{Float64, 1}  # mitigation effort costs
    Y::Array{Float64, 1}  # gross WP
    C::Array{Float64, 1}  # consumption
    E::Array{Float64, 1}  # ausstoss
    K::Array{Float64, 1}  # capital
    T::Array{Float64, 2}  # temperature
    M::Array{Float64, 2}  # carbon stock
end
