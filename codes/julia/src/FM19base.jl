# Automation and the Future of Work
# Michele Fornino and Andrea Manera

#---------------------------------------------------#
#					  BASE model 					#
#---------------------------------------------------#


"""
	FM19base_params{F<:Real, I<:Integer}

	Composite type containing all parameters of the Base model in Fornino Manera (2019).


"""
@with_kw struct FM19base_params{F<:Real, I<:Integer} @deftype F

	# Parameters
	θ = 0.5
	w = 0.5
	pR = 1.5
	Γ = 0.5
	m = 0.2
	δ = log(1 + 1/12)
	ρ = log(1 + 0.04)
	ψR = 25.0
	θz = 0.1
	σz = 0.4
	zd = 1.0
	p = 1.0

	# Derived Parameters
	Ω = (1.0 - Γ)/Γ * w - m
	Rmax =  1.0 ./ (δ .* ψR) .* (Ω ./ (ρ .+ δ) .- pR)

	# Grids Parameters
	NR::I = 200
	Nz::I = 200
	zGrid_lower_quantile = 1e-4
	zGrid_upper_quantile = 1 - 1e-4
	RGrid_lower = 0.0
	RGrid_upper = Rmax

	# Assertions
	@assert (θ > 0) & (θ < 1)
	@assert w > 0
	@assert pR > 0
	@assert (Γ > 0) & (Γ < 1)
	@assert m > 0
	@assert δ > 0
	@assert ρ > 0
	@assert ψR > 0
	@assert θz > 0 
	@assert σz > 0
	@assert zd > 0
	@assert p > 0
	@assert Ω > 0 "Labor savings Ω cannot be negative."

	# Check values for grids make sense
	@assert NR > 1
	@assert Nz > 0
	@assert zGrid_upper_quantile > zGrid_lower_quantile
	@assert zGrid_lower_quantile > 0
	@assert zGrid_upper_quantile < 1
	@assert RGrid_lower >= 0

	# Enforce definitions (prevents accidental creation of object with wrong derived parameters)
	@assert isequal(Ω, (1.0 - Γ)/Γ * w - m) "Ω is a derived parameter. Please modify the underlying parameters if you wish to change it."
	@assert isequal(Rmax, 1.0 ./ (δ .* ψR) .* (Ω ./ (ρ .+ δ) .- pR)) "Rmax is a derived parameter. Please modify the underlying parameters if you wish to change it."

end


"""
	function HamiltonJacobiBellmanEquations.HJBEquation(params::FM19base_params)


Overloaded constructor for input of type FM19base_params.

"""
function HamiltonJacobiBellmanEquations.HJBEquation(params::FM19base_params)

	# Use UnPack.jl to localize all members of params
	@unpack_FM19base_params params

	# Profit Function
	function Π(R::T, z::T)::T where T<:Real
		Rb = 1.0/ (1.0 - Γ) * ((p * z * θ * Γ) / w) ^ (1.0 / (1.0 - θ))
		Rh = 1.0/ (1.0 - Γ) * ((p * z * θ * (1.0 - Γ)) / m) ^ (1.0 / (1.0 - θ))
		if R ≤ Rb
			static_profit = (1.0 - θ) * p * z * ((1.0 - Γ) * Rb) ^ θ + Ω * R
		elseif Rb ≤ R ≤ Rh
			static_profit = p * z * ((1.0 - Γ) * R) ^ θ - m * R
		else
			static_profit = p * z * ((1.0 - Γ) * Rh) ^ θ - m * Rh
		end
		return static_profit
	end

	# Returns Function
	function flow_returns(state::AbstractVector{T}, I::Union{T, AbstractVector{T}})::T where T<:Real 
		return Π(state[1], state[2]) - pR * I[1] - ψR / 2 * I[1] ^ 2
	end

	# Policy Function in term of marginal value function
	function optimal_controls(state::AbstractVector{T}, ∂V∂R::Union{T, AbstractVector{T}})::T where T<:Real
		return 1 / ψR * (∂V∂R[1] - pR)
	end

	# State equations
	function state_equations(state::AbstractVector{T}, I::Union{T, AbstractVector{T}})::T where T<:Real
		return I[1] - δ * state[1]
	end

	# State constraint
	function state_constraints(state::AbstractVector{T})::T where T<:Real
		return ψR * δ * state[1] + pR
	end

	# Tuple of endogenous state grids
	RGrid = Vector(range(RGrid_lower, RGrid_upper, length = NR))
	state_grids_endog = (RGrid, )

	# Exogenous shock
	exog_shock = ContinuousTimeMarkovChain(params)
	
	# Allocate hjb object
	return HJBEquation(flow_returns, 
					   state_equations, 
					   optimal_controls, 
					   state_constraints,
					   state_grids_endog,
					   exog_shock,
					   ρ)

end



"""
	function MarkovChains.ContinuousTimeMarkovChain(params::FM19base_params)


Overloaded constructor for input of type FM19base_params.

"""
function MarkovChains.ContinuousTimeMarkovChain(params::FM19base_params)

	@unpack_FM19base_params params

	# Construct exog_shock ContinuousTimeMarkovChain object (if problem is stochastic).
	if Nz > 1

		# Obtain discretize diffusion process object as a MarkovChains.ContinuousTimeMarkovChain.
		μ = z -> -θz .* (z .- zd)
		σ = z -> sqrt(2 * θz * σz^2 / zd) .* sqrt.(z)
		dp = ItoDiffusionProcess(μ, σ)
		
		# Theoretical benchmark for CIR stationary distribution, used for zGrid.
		shape = zd^2 / σz^2
		rate = zd / shape
		Fz_theory = Gamma(shape, rate)
		zGrid_lower = quantile(Fz_theory, zGrid_lower_quantile)
		zGrid_upper = quantile(Fz_theory, zGrid_upper_quantile)

		# Allocate Grid
		zGrid = Vector(range(zGrid_lower, zGrid_upper, length = Nz))

		# Exogenous Shock: ContinuousTimeMarkovChain object 
		exog_shock = ContinuousTimeMarkovChain(dp, zGrid)

	elseif Nz == 1

		exog_shock = ContinuousTimeMarkovChain(reshape([0], 1, 1), [zd])

	end

	return exog_shock

end


function HamiltonJacobiBellmanEquations.solve(params::FM19base_params)

	return solve(HJBEquation(params))

end



# # IDEA?
# function plot(params::FM19base_params, sol::HJBSolution)

# 	# PLOTTING FOR 2-dim problem (1 end + 1 exo)
	
# 	# IDEA IS TO MAKE MULTIPLE FIGURES
# 	# FIG 1: 
# 	# [subplot surf/contour of value function ]
# 	# [subplot slice x] ||| [subplot slice z  ]

# 	# FIG 2:
# 	# same as above but one for each policy c_j

# 	# FIG 3:
# 	# stationary distribution (countour)


# 	# PLOTTING FOR N-dim problem (N-1 end + 1 exo)
# 	# THINK ABOUT HOW TO MAKE IT EFFECTIVELY
# end

# # IDEA?
# function summarize(sol::HJBSolution)

# 	# OBTAIN STATISTICS REGARDING SOLUTION OBJECTS (MUST BE GENERIC FOR ANY PROBLEM)

# end
