
# simulation function (no optimization)

function simulate!(ds)
    # unpack
    opts = ds.options
    exVar = ds.exogenousVariables
    iv = ds.initialValues
    #init
    Ω, Λ, Y, Q, C, E, K, U = (zeros(opts.N) for _=1:12)
    T = zeros(2, opts.N)
    M = zeros(3, opts.N)
    μ = iv.μ₀*ones(opts.N)
    S = iv.S₀*ones(opts.N)
    T[:,1] = iv.T₀
    M[:,1] = iv.M₀
    K[1] = iv.K₀
    Φ_K = (1-opts.δₖ)^opts.tStep
    for t  = 1:opts.N
        Ω[t] = (1-1/(1 + opts.a₂*T[1, t]^opts.a₃))  # damages
        Λ[t] = exVar.θ₁[t]*μ[t]^opts.θ₂  # mitigation
        Y[t] = exVar.A[t]*K[t]^opts.γ*(exVar.L[t]/1000)^(1-opts.γ)
        Q[t] = (1-Ω[t])*(1-Λ[t])*Y[t]
        C[t] = (1-S[t])*Y[t]
        E[t] = exVar.σ[t]*(1-μ[t])*Y[t] + exVar.E_LAND[t]
        U[t] = exVar.L[t]*((C[t]*1000.0/exVar.L[t])^(1-opts.α)-1)/(1-opts.α)
        if t<opts.N
            K[t+1] = Φ_K*K[t] + opts.tStep*Q[t]*S[t]  # [CAP]
            M[:,t+1] = opts.Φ_M*M[:,t] +  opts.b_M*E[t] # [CAR]
            T[:,t+1] = opts.Φ_T*T[:,t] + opts.b_T*(opts.η * log2(M[1,t]/opts.M_AT1750)+exVar.F_EX[t])  # [CLI]
        end
    end
    # compute social welfare
    W = sum(U[i]/((1+opts.ρ)^(opts.tStep*i)) for i=1:opts.N)
    # save to results
    years = 2005 .+ (opts.tStep*(1:opts.N))
    ds.results = Results(years, μ, S, Ω, Λ, Y, Q, C, E, K, T[1,:], T[2,:], M[1,:], M[2,:], M[3,:], U, W)
end

# optimization function

function optimization!(ds::DiceSimulation)

    # for the baseline scenario, no optimization is needed, hence we simply simulate
    if ds.scenario =="baseline"
        simulate!(ds)
    else
        # unpack input
        opts = ds.options
        exVar = ds.exogenousVariables
        iv = ds.initialValues
        N = opts.N
        μ_ubound = opts.μ_ubound*ones(N)
        cca_ubound = opts.fosslim
        # init
        model =  Model()
        # define Variables
        @variable(model, 0.0 <= μ[i=1:N] <= μ_ubound[i]); # Emission control rate GHGs
        @variable(model, FORC[1:N]); # Increase in radiative forcing (watts per m2 from 1900)
        @variable(model, 0.0 <= Tₐₜ[1:N] <= 40.0); # Increase temperature of atmosphere (degrees C from 1900)
        @variable(model, -1.0 <= Tₗₒ[1:N] <= 20.0); # Increase temperatureof lower oceans (degrees C from 1900)
        @variable(model, Mₐₜ[1:N] >= 10.0); # Carbon concentration increase in atmosphere (GtC from 1750)
        @variable(model, Mᵤₚ[1:N] >= 100.0); # Carbon concentration increase in shallow oceans (GtC from 1750)
        @variable(model, Mₗₒ[1:N] >= 1000.0); # Carbon concentration increase in lower oceans (GtC from 1750)
        @variable(model, E[1:N]); # Total CO2 emissions (GtCO2 per year)
        @variable(model, Eind[1:N]); # Industrial emissions (GtCO2 per year)
        @variable(model, C[1:N] >= 2.0); # Consumption (trillions 2005 US dollars per year)
        @variable(model, K[1:N] >= 1.0); # Capital stock (trillions 2005 US dollars)
        @variable(model, CPC[1:N] >= 0.01); #  Per capita consumption (thousands 2005 USD per year)
        @variable(model, I[1:N] >= 0.0); # Investment (trillions 2005 USD per year)
        @variable(model, S[1:N]>= 0.0); # Gross savings rate as fraction of gross world product
        @variable(model, RI[1:N]); # Real interest rate (per annum)
        @variable(model, Y[1:N] >= 0.0); # Gross world product (trillions 2005 USD per year)
        @variable(model, Q[1:N] >= 0.0); # Gross world product net of abatement and damages (trillions 2005 USD per year)
        @variable(model, YGROSS[1:N] >= 0.0); # Gross world product GROSS of abatement and damages (trillions 2005 USD per year)
        @variable(model, YNET[1:N]); # Output net of damages equation (trillions 2005 USD per year)
        @variable(model, DAMAGES[1:N]); # Damages (trillions 2005 USD per year)
        @variable(model, Ω[1:N]>=0.0); # Damages as fraction of gross output
        @variable(model, Λ[1:N]>=0.0); # Cost of emissions reductions  (trillions 2005 USD per year)
        @variable(model, MCABATE[1:N]); # Marginal cost of abatement (2005$ per ton CO2)
        @variable(model, CCA[1:N] <= cca_ubound); # Cumulative industrial carbon emissions (GTC)
        @variable(model, U[1:N]); # utility function
        @variable(model, W); # Welfare function

        # pack them into an object
        vars = Variables(μ, S, Ω, Λ, Y, Q, C, E, K, Tₐₜ, Tₗₒ, Mₐₜ, Mᵤₚ, Mₗₒ, U, W)
        # Temperature-climate equations
        @NLconstraint(model, [i=1:N-1], vars.Tₐₜ[i+1] == opts.Φ_T[1,1]*vars.Tₐₜ[i] +  opts.Φ_T[1,2]*vars.Tₗₒ[i] + opts.b_T[1]*(opts.η * log2(vars.Mₐₜ[i]/opts.M_AT1750)+exVar.F_EX[i]))
        @NLconstraint(model, [i=1:N-1], vars.Tₗₒ[i+1] == opts.Φ_T[2,1]*vars.Tₐₜ[i] +  opts.Φ_T[2,2]*vars.Tₗₒ[i])
        # carbon stock equations
        @constraint(model, [i=1:N-1], vars.Mₐₜ[i+1] == opts.Φ_M[1,1]*vars.Mₐₜ[i] + opts.Φ_M[1,2]*vars.Mᵤₚ[i] + opts.Φ_M[1,3]*vars.Mₗₒ[i]  + opts.b_M[1]*vars.E[i] )
        @constraint(model, [i=1:N-1], vars.Mᵤₚ[i+1] == opts.Φ_M[2,1]*vars.Mₐₜ[i] + opts.Φ_M[2,2]*vars.Mᵤₚ[i] + opts.Φ_M[2,3]*vars.Mₗₒ[i] )
        @constraint(model, [i=1:N-1], vars.Mₗₒ[i+1] == opts.Φ_M[3,1]*vars.Mₐₜ[i] + opts.Φ_M[3,2]*vars.Mᵤₚ[i] + opts.Φ_M[3,3]*vars.Mₗₒ[i] )
        # Capital balance equation
        @NLconstraint(model, [i=1:N-1], vars.K[i+1] == (1-opts.δₖ)^opts.tStep * K[i] + opts.tStep*vars.Q[i]*vars.S[i] )
        # Damages  euqation
        @NLconstraint(model, [i=1:N], vars.Ω[i] == 1 - 1/(1 + opts.a₂*Tₐₜ[i]^opts.a₃) )
        # Mitigation effort cost equation
        @NLconstraint(model, [i=1:N], vars.Λ[i] == exVar.θ₁[i]*vars.μ[i]^opts.θ₂)
        # Output equation
        @NLconstraint(model, [i=1:N], vars.Y[i] == exVar.A[i]*vars.K[i]^opts.γ*(exVar.L[i]/1000)^(1-opts.γ))
        # Net economic output equation
        @NLconstraint(model, [i=1:N], vars.Q[i] == (1 - vars.Ω[i])*(1-vars.Λ[i])*vars.Y[i])
        # consumption
        @NLconstraint(model, [i=1:N], vars.C[i] == (1-vars.S[i])*vars.Q[i])
        # emissions
        @NLconstraint(model, [i=1:N], vars.E[i] == exVar.σ[i]*(1-vars.μ[i])*vars.Y[i] + exVar.E_LAND[i])
        # Instantaneous utility function equation
        @NLconstraint(model, [i=1:N], vars.U[i] == exVar.L[i]*((vars.C[i]*1000.0/exVar.L[i])^(1-opts.α)-1)/(1-opts.α))
        # Social welfare function
        #scale1 = 0.016408662 #Multiplicative scaling coefficient
        #scale2 = -3855.106895 #Additive scaling coefficient
        #@constraint(model, vars.W == scale1 *sum(vars.U[i]/((1+opts.ρ)^(opts.tStep*i)) for i=1:N)+ scale2)
        @constraint(model, vars.W == sum(vars.U[i]/((1+opts.ρ)^(opts.tStep*i)) for i=1:N))


        # set initial conditions
        JuMP.fix(vars.Mₐₜ[1], iv.M₀[1]; force=true)
        JuMP.fix(vars.Mᵤₚ[1], iv.M₀[2]; force=true)
        JuMP.fix(vars.Mₗₒ[1], iv.M₀[3]; force=true)
        JuMP.fix(vars.Tₗₒ[1], iv.T₀[2]; force=true)
        #JuMP.fix(vars.Tₐₜ[1], iv.T₀[1]; force=true)
        # JuMP.fix(vars.K[1], iv.K₀; force=true);
        @NLconstraint(model, vars.K[1] == iv.K₀)
        @constraint(model, vars.Tₐₜ[1] == iv.T₀[1])

        if ds.scenario=="twoDegree"# && N<=19
            for i in 2:N
                JuMP.set_upper_bound(vars.Tₐₜ[i], 2.0);
            end
#        elseif ds.scenario=="twoDegree" && N>19
#            throw("Only works until 2100")
        end

        # Objective function
        @objective(model, Max, vars.W);
        # set optimizer
        set_optimizer(model, Ipopt.Optimizer)
        set_optimizer_attribute(model, "print_level", 5)
        set_optimizer_attribute(model,  "sb",  "yes")
        set_optimizer_attribute(model,  "start_with_resto",  "yes")
        set_optimizer_attribute(model,  "expect_infeasible_problem",  "yes")
        set_optimizer_attribute(model,  "max_iter",  5000)
        set_optimizer_attribute(model, "tol", 0.1)
        # actual optimization
        optimize!(model)
        # save results to ds
        ds.results= model_results(vars, opts)
        return ds
    end
end
# Equation for MC abatement
#@NLconstraint(model, [i=1:N], vars.MCABATE[i] == params.pbacktime[i] * vars.μ[i]^(opts.θ₂-1));
# Carbon price equation from abatement
#@NLconstraint(model, [i=1:N], vars.CPRICE[i] == params.pbacktime[i] * (vars.μ[i]/params.partfract[i])^(opts.θ₂-1));
