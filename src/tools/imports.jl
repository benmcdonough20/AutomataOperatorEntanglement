using Pkg
Pkg.activate("automata")
using Plots
using RandomQuantum
using StatsBase
using LinearAlgebra
using Random
using LaTeXStrings
using ProgressBars
using SpecialFunctions
using TensorKit
using ProgressBars
using RandomMatrix
using LsqFit
using CurveFit
using Arrow
using Kronecker

include("../tools/pauli_tools.jl")
include("../tools/permutation_tools.jl")
include("../tools/ES_tools.jl")
include("../tools/distribution_tools.jl")
include("../tools/plotting_defaults.jl")