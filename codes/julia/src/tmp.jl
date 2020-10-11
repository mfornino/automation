# Temporary module
module tmp

# Include HamiltonJacobiBellmanEquations as if it were a package.
# include("hjb.jl")
# using .HamiltonJacobiBellmanEquations

# # Dependencies specific to this code:
using MarkovChains
using DynamicProgramming
using Distributions # used to get boundaries of zGrid from theoretical Gamma
using Parameters
# # Temporary Modules
# using Debugger


#-----------------------------------------------------------------------------#
# 							Run RCK Model Test 								  #
#-----------------------------------------------------------------------------#
# NOTE: This is a deterministic model, but can be made stochastic by setting Nz > 1.

# include("RamseyCassKoopmansTest.jl")
# parameters = RCK_params()
# display(parameters)
# problem = HJBEquation(parameters)
# solution = solve(problem)
# statDist = stationary_distribution(solution) # just for the spike :)


#-----------------------------------------------------------------------------#
# 							Run Aiyagari Model Test 						  #
#-----------------------------------------------------------------------------#
# include("AiyagariTest.jl")
# parameters = Aiyagari_params()
# display(parameters)
# # V0 = reshape(1:40000, parameters.Na, parameters.Nz)
# # shock = ContinuousTimeMarkovChain(parameters)
# problem = HJBEquation(parameters)
# solution = solve(problem)
# statDist = stationary_distribution(solution)
# samplePath = random_sample(solution)


#-----------------------------------------------------------------------------#
# 							Run FM Base Model Test 							  #
#-----------------------------------------------------------------------------#
include("FM19base.jl")
parameters = FM19base_params()
display(parameters)
# shock = ContinuousTimeMarkovChain(parameters)
problem = HJBEquation(parameters)
solution = solve(problem)
statDist = stationary_distribution(solution)
samplePath = random_sample(solution)


#-----------------------------------------------------------------------------#
# 						Run FM Rigid Labor Model Test 						  #
#-----------------------------------------------------------------------------#
# include("FM19rigidlabor.jl")
# parameters = FM19rigidlabor_params()
# shock = ContinuousTimeMarkovChain(parameters)
# problem = HJBEquation(parameters)
# sol = solve(problem)


end # Module





