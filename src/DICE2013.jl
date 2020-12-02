module DICE2013

using JuMP
using Ipopt
using Plots


include("types.jl")
include("functions.jl")
include("optimization.jl")

export setOptions, setInitialValues, optimization!, plotDice, DiceSimulation

end # module
