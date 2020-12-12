# function to set the initial values. If called without argument, standard 2013
# values are used
function setInitialValues(;
    μ₀::Float64 = 0.039,
    S₀::Float64 =0.25,
    L₀::Int = 6838,
    A₀::Float64 = 3.8,
    σ₀::Float64 = 0.5491,
    T₀::Array{Float64,1} = [0.8; 0.0068],
    M₀::Array{Float64,1} = [830.4; 1527;10010],
    K₀::Float64 = 135.0
    )
    InitialValues(μ₀,S₀, L₀, A₀, σ₀, T₀, M₀, K₀)
end

# function to set parameters of the model. If called without argument, standard
# 2013 values are used
function setOptions(;
    tStep = 5,
    N = 60,
    N100 = 19,
    Lₐ = 10500, # in million people
    l_g = 0.134,
    δₐ = 0.006,
    gₐ = 0.079,
    δσ = 0.001,
    gσ = 0.01,
    f₀ = 0.25,
    f₁ = 0.7,
    E_L0 = 3.3,
    δEL = 0.2,
    Φ11 = 0.8630,
    Φ12 = 0.0086,
    Φ21 = 0.025,
    Φ22 = 0.975,
    Φ_T = [Φ11 Φ12; Φ21 Φ22],
    ξ1 = 0.098,
    b_T = [ξ1; 0],
    M_AT1750 = 588,
    η = 3.8,
    ζ11 = 0.912,
    ζ12 = 0.03833,
    ζ21 = 0.088,
    ζ22 = 0.9592,
    ζ23 = 0.0003375,
    ζ32 = 0.00250,
    ζ33 = 0.9996625,
    ξ2 = 5/3.666,
    Φ_M = [ζ11 ζ12 0; ζ21 ζ22 ζ23; 0 ζ32 ζ33],
    b_M = [ξ2;0;0],
    a₂ = 0.00267,
    a₃ = 2,
    θ₂ = 2.8,
    δₖ = 0.1,
    γ = 0.3,
    pb = 344,
    δpb = 0.025,
    ρ = 0.015,
    limμ = 1.2,
    fosslim = 6000.0,  #Maximum cumulative extraction fossil fuels (GtC)
    α = 1.45)
    Options(tStep, N, N100, Lₐ, l_g, δₐ, gₐ, δσ, gσ, f₀, f₁, E_L0, δEL, Φ_T, b_T, M_AT1750, η, Φ_M, b_M, a₂, a₃,  θ₂, δₖ, γ, pb, δpb, ρ, limμ, fosslim, α)
end

# function that precomputes the exogenous variables, as they do not depend on
# the control variables
function computeExogenousVariables(c, iv)
    # init
    L, A, σ, F_EX, E_LAND, θ₁ = (zeros(c.N) for _ = 1:6)
    L[1] = iv.L₀
    A[1] = iv.A₀
    σ[1] = iv.σ₀
    # variables depending on previous step
    for t = 1:c.N-1
        L[t+1] = L[t]* ((1 + c.Lₐ)/(1+L[t]))^c.l_g  # [POP]
        A[t+1] = A[t]/ (1-c.gₐ*exp(-c.δₐ*c.tStep*(t-1)))  # [TFP]
        σ[t+1] = σ[t]*exp(-c.gσ*c.tStep*(1-c.δσ)^(c.tStep*(t-1)) )  # [EI]
    end
    # variables not depending on previous step
    F_EX = [c.f₀ + min(c.f₁-c.f₀, (c.f₁-c.f₀)*(t-1)/c.N100) for t=1:c.N]
    E_LAND = [c.E_L0*(1-c.δEL)^(t-1) for t=1:c.N]
    θ₁ = [c.pb/(1000*c.θ₂) *(1-c.δpb).^(t-1)* σ[t] for t=1:c.N]
    ExogenousVariables(L, A, σ, F_EX, E_LAND, θ₁)
end

# simple plot function. The option "Scenario", plots all interesting variables
function plotDice(ds, plotvar)

    # set times
    actualTimes = range(2010,step=ds.options.tStep, length=ds.options.N)
    if plotvar=="L"
        fn = plot(actualTimes, ds.exogenousVariables.L, ylabel="population (in Million)",xlabel="years",title="Labour",legend=false)
        #savefig(fn, "C:/Users/andre/Google Drive/03_Dokumente/DICE Modell/pictures/labour.tex")
    elseif plotvar == "A"
        fn = plot(actualTimes, ds.exogenousVariables.LA, ylabel="Efficiency",xlabel="years",title="Total Factor Productivity",legend=false)
        #savefig(fn, "C:/Users/andre/Google Drive/03_Dokumente/DICE Modell/pictures/facProd.tex")
    elseif plotvar == "σ"
        fn = plot(actualTimes, ds.exogenousVariables.Lσ, ylabel="intensity",xlabel="years",title="Carbon Intensity",label = "sigma")
        #savefig(fn, "C:/Users/andre/Google Drive/03_Dokumente/DICE Modell/pictures/carbonIntensity.tex")
    elseif plotvar == "E_LAND"
        fn = plot(actualTimes, ds.exogenousVariables.LE_LAND, ylabel="climate forcing",xlabel="years",title="Land Conversion",label = "E")
        #savefig(fn, "C:/Users/andre/Google Drive/03_Dokumente/DICE Modell/pictures/landConversion.tex")
    elseif plotvar == "F_EX"
        fn = plot(actualTimes, ds.exogenousVariables.LF_EX, ylabel="climate forcing",xlabel="years",title="GHG Except Carbon",label = "FEX")
        # savefig(fn, "C:/Users/andre/Google Drive/03_Dokumente/DICE Modell/pictures/externalForcing.tex")
    elseif plotvar == "θ₁"
        fn = plot(actualTimes, ds.exogenousVariables.Lθ₁, ylabel="percentage cost",xlabel="years",title="Cost Of Mitigation",label = "theta1")
        # savefig(fn, "C:/Users/andre/Google Drive/03_Dokumente/DICE Modell/pictures/mitigationCost.tex")
    elseif plotvar == "damages"
        t = -3:0.1:10
        fn = plot(t, 1 .- 1 ./(1 .+ ds.options.a₂ .*t.^2), xlabel="temperature",ylabel="damages (in p of GWP)", title="Damages")
    elseif plotvar == "Scenario"
        # scenario plot
        fR1 = plot(actualTimes, ds.results.μ, ylabel="effort taken",xlabel="years",title="Mitigation Efforts",label = "mu")
        fR2 = plot(actualTimes,  ds.results.S, ylabel="rate in percent",xlabel="years",title="Saving Rate",label = "s", ylim = [0,1])
        f1 = plot(actualTimes, ds.results.Tₐₜ, ylabel="degree C (compared to 1900)",xlabel="years",title="Temperature",label = "TAT")
        plot!(actualTimes,  ds.results.Tₗₒ, label = "TLO")
        f2 = plot(actualTimes,  ds.results.Mₐₜ,  ylabel="carbon (GTC)",xlabel="years",title="Carbon Stock",label = "MAT")
        plot!(actualTimes,  ds.results.Mᵤₚ, label = "MUP")
        f3 = plot(actualTimes,  ds.results.Y,  ylabel="USD (trillion)",xlabel="years",title="Gross World Product",label = "Y")
        plot!(actualTimes,  ds.results.C, label ="Q")
        plot!(actualTimes,  ds.results.Λ, label = "Λ")
        plot!(actualTimes,  ds.results.Ω, label = "Ω")
        f4 = plot(actualTimes,  ds.results.E,  ylabel="carbon (GTC)",xlabel="years",title="Emissions",label = "E")
        plot!(actualTimes,  ds.exogenousVariables.E_LAND, label = "ELand")
        l = @layout [a b; c d; e f]
        fn = plot(fR1, fR2, f1, f2, f3, f4, layout = l)
        gui(fn)
        #savefig(fG, "C:/Users/andre/Google Drive/03_Dokumente/DICE Modell/pictures/baseLineScenario2.tikz")
    else
        throw("plot argument not defined")
    end
    display(fn)
end

# function that stores the JuMP variables, that contain the optimization
# results in variables for output
function model_results(vars::Variables, opts::Options)
    years = 2005 .+ (opts.tStep*(1:opts.N))
    μ = value.(vars.μ);
    S = value.(vars.S);
    Ω = value.(vars.Ω);
    Λ = value.(vars.Λ);
    Y = value.(vars.Y);
    Q = value.(vars.Q);
    C = value.(vars.C);
    E = value.(vars.E);
    K = value.(vars.K);
    Tₐₜ = value.(vars.Tₐₜ);
    Tₗₒ = value.(vars.Tₗₒ);
    Mₐₜ = value.(vars.Mₐₜ);
    Mᵤₚ = value.(vars.Mᵤₚ);
    Mₗₒ = value.(vars.Mₗₒ);
    #Eind = value.(vars.Eind);
    U = value.(vars.U);
    W = value.(vars.W);
    #1000 .* shadow_price.(eqs.eeq)./shadow_price.(eqs.cc)
    Results(years, μ, S, Ω, Λ, Y, Q, C, E, K, Tₐₜ, Tₗₒ, Mₐₜ, Mᵤₚ, Mₗₒ, U, W)
end

# standard simulation function
function runScenario(sc)
    opts = setOptions()  # set parameters
    iv = setInitialValues()  # set initial values
    if ~ (sc in ["baseline", "optimal", "2degree"])
        throw("scenario not known. Possible scenarios are baseline, optimal, and 2degree")
    end
    ds = DiceSimulation(sc, opts, iv)   # construct simulation object
    optimization!(ds)  # optimize using Ipopt and MUMPS
    plotDice(ds, "Scenario")  # plot the scenario
end
