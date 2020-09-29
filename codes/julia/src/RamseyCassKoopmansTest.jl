# Aiyagari Model Test
# Michele Fornino
#---------------------------------------------------#
#			Ramsey-Cass-Koopmans MODEL 				#
#---------------------------------------------------#


"""
	RCK_params{F<:Real, I<:Integer}

	Composite type containing all parameters of the standard Ramsey-Cass-Koopmans model.


"""
@with_kw struct RCK_params{F<:Real, I<:Integer} @deftype F

	# Parameters
	σ = 1.0
	n = 0.0
	α = 0.36
	δ = log(1 + 0.1)
	ρ = log(1 + 0.04)
	# σz = 0.2
	# θz = 0.1
	zd = 1.0

	# Derived Parameters
	γ = 1.0 / σ

	# Grids Parameters
	Nk::I = 2000
	Nz::I = 1
	# zGrid_lower_quantile = 1e-2
	# zGrid_upper_quantile = 1 - 1e-2
	kGrid_lower = 0.01
	# kGrid_upper = 20.
	kGrid_upper = 1.5 * (zd * α / (n + δ + ρ)) ^ (1 / (1 - α))

	# Assertions
	@assert σ > 0.0
	@assert n >= 0.0
	@assert δ > 0.0
	@assert ρ > 0.0
	@assert (α > 0.0) & (α < 1.0) 
	# @assert θz > 0 
	# @assert σz > 0
	@assert zd > 0

	# Check values for grids make sense
	@assert Nk > 1
	@assert Nz > 0
	# @assert zGrid_upper_quantile > zGrid_lower_quantile
	# @assert zGrid_lower_quantile > 0.0
	# @assert zGrid_upper_quantile < 1.0
	@assert kGrid_upper > kGrid_lower
	@assert kGrid_lower >= 0.

	# Enforce definitions (prevents accidental creation of object with wrong derived parameters)
	@assert γ == 1 / σ

end


"""
	function HamiltonJacobiBellmanEquations.HJBEquation(params::RCK_params)


Overloaded constructor for input of type RCK_params.

"""
function HamiltonJacobiBellmanEquations.HJBEquation(params::RCK_params)

	# Use UnPack.jl to localize all members of params
	@unpack_RCK_params params

	# Utility Function and derivatives.
	function u(c::T)::T where T<:Real
		if σ == 1.0
			return log(c)
		else
			return (c ^ (1 - γ) - 1) / (1 - γ)
		end
	end

	function u′(c::T)::T where T<:Real
		return c ^ -γ
	end

	function u′⁻¹(x::T)::T where T<:Real
		return x ^ (- 1 / γ)
	end

	function f(k::T, z::T)::T where T<:Real
		return z * (k ^ α)
	end

	# function u(c::T)::T where T<:Real
	# 	return - 1 / γ * exp(-γ * c)
	# end

	# function u′(c::T)::T where T<:Real
	# 	return exp( - γ * c)
	# end

	# function u′⁻¹(x::T)::T where T<:Real
	# 	return - 1 / γ * log(x)
	# end

	# Returns Function
	function flow_returns(state::AbstractVector{T}, c::Union{T, AbstractVector{T}})::T where T<:Real 
		return u(c[1])
	end

	# Policy Function in term of marginal value function
	function optimal_controls(state::AbstractVector{T}, ∂V∂k::Union{T, AbstractVector{T}})::T where T<:Real
		return u′⁻¹(∂V∂k[1])
	end

	# State equations
	function state_equations(state::AbstractVector{T}, c::Union{T, AbstractVector{T}})::T where T<:Real
		return f(state[1], state[2]) - (n + δ) * state[1] - c[1]
	end

	# State constraint
	function state_constraints(state::AbstractVector{T})::T where T<:Real
		return u′(f(state[1], state[2]) - (n + δ) * state[1])
	end

	# Tuple of endogenous state grids
	kGrid = Vector(range(kGrid_lower, kGrid_upper, length = Nk))
	state_grids_endog = (kGrid, )

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
	function MarkovChains.ContinuousTimeMarkovChain(params::RCK_params)


Overloaded constructor for input of type RCK_params.

"""
function MarkovChains.ContinuousTimeMarkovChain(params::RCK_params)

	@unpack_RCK_params params

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

