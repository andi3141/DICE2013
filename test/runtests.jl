using DICE2013
using JuMP
using Ipopt


opts = setOptions()
iv = setInitialValues()
dsBase = DiceSimulation("baseline", opts, iv)
optimization!(dsBase)
plotDice(dsBase, "Scenario")
