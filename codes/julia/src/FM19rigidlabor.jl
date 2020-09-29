# Automation and the Future of Work
# Michele Fornino and Andrea Manera

#---------------------------------------------------#
#			     RIGID LABOR model 					#
#---------------------------------------------------#


"""
	FM19rigidlabor_params{F<:Real, I<:Integer}

	Composite type containing all parameters of the model with rigid labor input detailed in Fornino Manera (2019).


"""
@with_kw struct FM19rigidlabor_params{F<:Real, I<:Integer} @deftype F

	# Parameters
	θ = 0.5
	w = 0.5
	pR = 1.5
	Γ = 0.5
	m = 0.2
	δ = log(1 + 1/12)
	s = δ
	ρ = log(1 + 0.04)
	ψR = 25.0
	ψL = 1.0
	θz = 0.4
	σz = 0.1
	zd = 1.0
	p = 1.0

	# Derived Parameters
	# Ω = (1.0 - Γ)/Γ * w - m
	Rmax =  1.0 ./ (δ .* ψR) .* (Ω ./ (ρ .+ δ) .- pR)

	# Grids Parameters
	NR::I = 50
	NL::I = 50
	Nz::I = 25
	zGrid_lower_quantile = 1e-4
	zGrid_upper_quantile = 1 - 1e-4
	RGrid_lower = 0.0
	RGrid_upper = 5 * Rmax
	LGrid_lower = 0.0
	LGrid_upper = ((3.0 * p * θ * Γ) / w)^(1.0 / (1.0 - θ)) / 3

	# Assertions
	@assert (θ > 0) & (θ < 1)
	@assert w > 0
	@assert pR > 0
	@assert (Γ > 0) & (Γ < 1)
	@assert m >= 0
	@assert δ > 0
	@assert ρ > 0
	@assert ψR > 0
	@assert θz > 0 "If you intended to use zero variance, just set Nz = 1."
	@assert σz > 0
	@assert zd > 0
	@assert p > 0

	# Check values for grids make sense
	@assert NR > 1
	@assert NL > 1
	@assert Nz > 0
	@assert zGrid_upper_quantile > zGrid_lower_quantile
	@assert zGrid_lower_quantile > 0
	@assert zGrid_upper_quantile < 1
	@assert RGrid_lower >= 0
	@assert LGrid_lower >= 0
	@assert RGrid_upper > RGrid_lower
	@assert LGrid_upper > LGrid_lower

	# Enforce definitions (prevents accidental creation of object with wrong derived parameters)
	# @assert isequal(Ω, (1.0 - Γ)/Γ * w - m) "Ω is a derived parameter. Please modify the underlying parameters if you wish to change it."

end


"""
	function HamiltonJacobiBellmanEquations.HJBEquation(params::FM19rigidlabor_params)


Overloaded constructor for input of type FM19rigidlabor_params.

"""
function HamiltonJacobiBellmanEquations.HJBEquation(params::FM19rigidlabor_params)

	# Use UnPack.jl to localize all members of params
	@unpack_FM19rigidlabor_params params

	# Profit Function
	function Π(R::T, L::T, z::T)::T where T<:Real
		Rh = 1.0/ (1.0 - Γ) * ((p * z * θ * (1.0 - Γ)) / m) ^ (1.0 / (1.0 - θ))
		if Γ * L ≤ (1.0 - Γ) * (Rh - R)
			static_profit = p * z * (Γ * L + (1.0 - Γ) * R) ^ θ - m * R - w * L
		elseif (1.0 - Γ) * (Rh - R) ≤ Γ * L ≤ (1.0 - Γ) * Rh
			static_profit = p * z * ((1.0 - Γ) * Rh) ^ θ - m * Rh - Γ / (1 - Γ) * Ω * L
		else
			static_profit = p * z * (Γ * L) ^ θ - w * L
		end
		return static_profit
	end

	# Returns Function
	function flow_returns(state::AbstractVector{T}, I::AbstractVector{T})::T where T<:Real 
		return Π(state[1], state[2], state[3]) - pR * I[1] - ψR / 2 * I[1] ^ 2 - ψL / 2 * I[2] ^ 2
	end

	# Policy Function in term of marginal value function
	function optimal_controls(state::AbstractVector{T}, 
							  ∇ₓV::AbstractVector{T})::AbstractVector{T} where T<:Real
		Iᵣ = 1 / ψR * (∇ₓV[1] - pR)
		Iₗ = ∇ₓV[2] / ψL 
		return [Iᵣ, Iₗ]
	end

	# State equations
	function state_equations(state::AbstractVector{T},
							I::AbstractVector{T})::AbstractVector{T} where T<:Real
		∂ₜR = I[1] - δ * state[1]
		∂ₜL = I[2] - s * state[2]
		return [∂ₜR, ∂ₜL]
	end

	# State constraint
	function state_constraints(state::AbstractVector{T})::AbstractVector{T} where T<:Real
		∂ᵣV = pR + ψR * δ * state[1]
		∂ₗV = ψL * s * state[2]
		return [∂ᵣV, ∂ₗV]
	end

	# Tuple of endogenous state grids
	RGrid = Vector(range(RGrid_lower, RGrid_upper, length = NR))
	LGrid = Vector(range(LGrid_lower, LGrid_upper, length = NL))
	state_grids_endog = (RGrid, LGrid)

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
	function MarkovChains.ContinuousTimeMarkovChain(params::FM19rigidlabor_params)


Overloaded constructor for input of type FM19rigidlabor_params.

"""
function MarkovChains.ContinuousTimeMarkovChain(params::FM19rigidlabor_params)

	@unpack_FM19rigidlabor_params params

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


function HamiltonJacobiBellmanEquations.solve(params::FM19rigidlabor_params)

	return solve(HJBEquation(params))

end
