# Solution of Hamilton-Jacobi-Bellman Equations
# Michele Fornino

module HamiltonJacobiBellmanEquations


# Dependencies
using LinearAlgebra
using SparseArrays
using MarkovChains
using ProgressMeter

import Base: 
	show

# Temporary Modules
# using Distributed
# using BenchmarkTools
# using Debugger
using PyPlot

export HJBEquation
export solve


#---------------------------------------------------------------------
# 						HJB EQUATION TYPE
# --------------------------------------------------------------------

struct HJBEquation
	
	# Basic objects
	discount_rate::Real									# Discount factor
	flow_returns::Function										# Returns function
	state_equations::Function 		 							# State Equations
	optimal_controls::Function 									# Optimal choice of controls from FOC(s)
	state_constraints::Function 								# Function to impose state constraints
	exog_infinitesimal_generator::AbstractMatrix{<:Real}	# Jump-Diffusion infinitesimal generator
	state_grids::Tuple{Vararg{AbstractVector{<:Real}}}					# Grids state variables
	
	# Useful objects
	N::Integer									# Number of endogenou states
	M::Integer									# Number of control variables
	nstates::Integer
	dx::Tuple{Vararg{<:Real}}	# diff(state_grids)
	state_grids_lengths::Tuple{Vararg{<:Integer}}
	state_grids_strides::Tuple{Vararg{<:Integer}}
	mc_state_grid::AbstractArray{<:AbstractVector{<:Real}}

# KEEP GOING

	# Constructor
	function HJBEquation(flow_returns::Function, 
						 state_equations::Function,
						 optimal_controls::Function,
						 state_constraints::Function,
						 state_grids_endog::Tuple{Vararg{AbstractVector{<:Real}}},
						 exog_shock::ContinuousTimeMarkovChain, 
						 discount_rate::Real)

		# Implement checks as needed
		
		state_grids = (state_grids_endog..., exog_shock.states)
		state_grids_lengths = length.(state_grids)
		state_grids_strides = Int.(cumprod(state_grids_lengths) ./ state_grids_lengths)
		N = length(state_grids_endog)
		M = length(optimal_controls([state_i[1] for state_i in state_grids], 
							 		zeros(Float64, N)))
		nstates = prod([n for n in state_grids_lengths])
		nstates_endog = Int.(nstates / state_grids_lengths[end])
		exog_infinitesimal_generator = kron(exog_shock.infinitesimal_generator,
											sparse(I, nstates_endog, nstates_endog))
		dx = map(x -> x[2] - x[1], state_grids[1:N])
		mc_state_grid = collect.(Iterators.product(state_grids...))

		new(discount_rate,
			flow_returns, 
			state_equations, 
			optimal_controls, 
			state_constraints,
			exog_infinitesimal_generator, 
			state_grids,
			N,
			M,
			nstates,
			dx,
			state_grids_lengths,
			state_grids_strides,
			mc_state_grid)
	end
end


# Overloading the show operator to get pretty-printing behavior.
function Base.show(io::IO, ::MIME"text/plain", hjb::HJBEquation) 
	print(io, "Hamilton-Jacobi-Bellman Equation.\n")
	print(io, "\nUse .field notation to get access to members.\n")
end


#---------------------------------------------------------------------
# 						HJB SOLUTION TYPE
# --------------------------------------------------------------------
struct HJBSolution
	infinitesimal_generator::AbstractMatrix{<:Real}
	mc_state_grid::AbstractArray{<:AbstractVector{<:Real}}
	value_function::AbstractArray{<:Real}
	policy_functions::AbstractArray{<:Real}
end

# Overloading the show operator to get pretty-printing behavior.
function Base.show(io::IO, ::MIME"text/plain", hjb::HJBSolution) 
	print(io, "Solution of Hamilton-Jacobi-Bellman Equation.\n")
	print(io, "\nUse .field notation to get access to members.\n")
end

# Wrapper for stationary_distribution
function MarkovChains.stationary_distribution(sol::HJBSolution)
	return stationary_distribution(ContinuousTimeMarkovChain(sol.infinitesimal_generator, sol.mc_state_grid))
end

# Wrapper for random_sample with random guess according to statdist if not provided.
function MarkovChains.random_sample(sol::HJBSolution)
	idx_start = rand(stationary_distribution(sol))
	return random_sample(sol, idx_start)
end

function MarkovChains.random_sample(sol::HJBSolution, idx_start::Integer)
	return random_sample(ContinuousTimeMarkovChain(sol.infinitesimal_generator, sol.mc_state_grid), idx_start)
end


	
#---------------------------------------------------------------------
# 						HJB SOLVE ROUTINE
# --------------------------------------------------------------------

# Wrapper for case in which guess is not provided. Use naive guess of controller who sets drifts to zero.
function solve(hjb::HJBEquation)
	
	# Construct standard guess based on "stay put" policy.
	∇V₀ = hjb.state_constraints.(hjb.mc_state_grid[:])
	c₀ = hjb.optimal_controls.(hjb.mc_state_grid[:], ∇V₀)
	guess = hjb.flow_returns.(hjb.mc_state_grid[:], c₀) ./ hjb.discount_rate
	
	return solve(hjb, guess)
end


function solve(hjb::HJBEquation, 
			   guess::AbstractArray;  # THIS NEEDS TO BE OPTIONAL PARAMETER
			   step_size::Tr=1e3, 
			   tol::Tr=1e-8, 
			   maxit::Ti where Ti<:Integer=1000) where Tr<:Real

	# Parameters for loop
	Δ = step_size
	ρ = hjb.discount_rate

	# GUESS OF THE VALUE FUNCTION
	V_next = guess

	# PREALLOCATE OBJECTS
	infinitesimal_generator = sparse([], [], [], hjb.nstates, hjb.nstates)
	policy_functions = Matrix{Float64}(undef, hjb.nstates, hjb.M)
	
	c₊ = Matrix{Float64}(undef, hjb.nstates, hjb.M)
	c₋ = similar(c₊)
	
	ẋ₊ = Matrix{Float64}(undef, hjb.nstates, hjb.N)
	ẋ₋ = similar(ẋ₊)

	up = Matrix{Bool}(undef, hjb.nstates, hjb.N)
	dwn = similar(up)
	stay = similar(up)

	R = Vector{Float64}(undef, hjb.nstates)

	lin_idx = Vector{Int64}(1:hjb.nstates)	

	function controls(states::Vector{<:Int64}, ∇V::Matrix{Float64})::Matrix{Float64}
		c = Matrix{Float64}(undef, length(states), hjb.M)
		for state in states
			c[state, :] .= hjb.optimal_controls(hjb.mc_state_grid[state], ∇V[state,:])
		end
		return c
	end

	function drifts(states::Vector{<:Int64}, c::Matrix{Float64})::Matrix{Float64}
		ẋ = Matrix{Float64}(undef, length(states), hjb.N)
		dx = [dx for dx in hjb.dx[1:hjb.N]]
		for state in states
			ẋ[state, :] .= hjb.state_equations(hjb.mc_state_grid[state], c[state, :]) ./ dx
		end
		# ẋ[abs.(ẋ) .< 1e-16] .= 0.0
		return ẋ
	end

	function returns(states::Vector{<:Int64}, c::Matrix{Float64})::Vector{Float64}
		R = Vector{Float64}(undef, length(states))
		for state in states
			R[state] = hjb.flow_returns(hjb.mc_state_grid[state], c[state, :])
		end
		return R
	end

	function constraints(states::Vector{<:Int64})
		∇R = Matrix{Float64}(undef, length(states), hjb.N)
		for state in states
			∇R[state, :] .= hjb.state_constraints(hjb.mc_state_grid[state])
		end
		return ∇R
	end

	# Enforce state constraints here!
	∇V₀ = constraints(lin_idx)
	∇V₀_cart = reshape(∇V₀, (hjb.state_grids_lengths..., hjb.N)...)
	∇V₊_cart = copy(∇V₀_cart)
	∇V₋_cart = copy(∇V₀_cart)
	c₀ = controls(lin_idx, ∇V₀)

	# Create sparse matrix rows, cols, vals triplet. Start with shock
	rows = Vector{Int64}(undef, 3 * (hjb.N+1) * hjb.nstates)
	cols = similar(rows)
	vals = Vector{Float64}(undef, 3 * (hjb.N+1) * hjb.nstates)
	r_mc , c_mc, v_mc = findnz(hjb.exog_infinitesimal_generator)
	len_mc = length(r_mc)
	ran = 1:len_mc
	rows[ran] = view(r_mc, :)
	cols[ran] = view(c_mc, :)
	vals[ran] = view(v_mc, :)
	
	# Initialize iterators
	p = ProgressThresh(log10(tol), desc="Iterating until convergence. Showing log₁₀(supnorm(Vₙ₊₁, Vₙ)): ")
	sup = 1e8
	it = 0
	while it < maxit

		# Update		
		it += 1
		V = reshape(V_next, hjb.state_grids_lengths...)

		# Compute marginal value for each endogenous state, imposing the state constraints
		for dim in 1:hjb.N

			view_∇V₊_cart = selectdim(∇V₊_cart[axes(V)..., dim], dim, 1:hjb.state_grids_lengths[dim]-1)
			view_∇V₋_cart = selectdim(∇V₋_cart[axes(V)..., dim], dim, 2:hjb.state_grids_lengths[dim])

			∇V₊_cart[parentindices(view_∇V₊_cart)..., dim] = (diff(V, dims = dim) ./ hjb.dx[dim])[:]
			∇V₋_cart[parentindices(view_∇V₋_cart)..., dim] = (diff(V, dims = dim) ./ hjb.dx[dim])[:]
			
			view_∇V₊_cart = selectdim(∇V₊_cart[axes(V)..., dim], dim, hjb.state_grids_lengths[dim])
			view_∇V₋_cart = selectdim(∇V₋_cart[axes(V)..., dim], dim, 1)
		end

		c₊ = controls(lin_idx, reshape(∇V₊_cart, (hjb.nstates, hjb.N)...))
		c₋ = controls(lin_idx, reshape(∇V₋_cart, (hjb.nstates, hjb.N)...))

		ẋ₊ = drifts(lin_idx, c₊)
		ẋ₋ = drifts(lin_idx, c₋)

		# Upwinding 
		up .= ẋ₊ .> 0.0
		dwn .= ẋ₋ .< 0.0
		stay .= (.~ up) .& (.~ dwn)

		# Construct upwind scheme
		# ∇V = reshape(∇V₊_cart, (hjb.nstates, hjb.N)...) .* up + 
		# 	 reshape(∇V₋_cart, (hjb.nstates, hjb.N)...) .* dwn + 
		# 	 ∇V₀ .* stay

		# Optimal choice of controls
		# c = controls(lin_idx, ∇V)
		c = c₊ .* up .+ c₋ .* dwn .+ c₀ .* stay

		# Maximized returns function
		R = returns(lin_idx, c)
		
		accumul = len_mc
		for dim in 1:hjb.N

			# Going UP
			up_dim = view(up, :, dim)
			len = sum(up_dim)
			ran = accumul .+ (1:len)
			rows[ran] = view(lin_idx, up_dim)
			cols[ran] = view(lin_idx, up_dim) .+ hjb.state_grids_strides[dim]
			vals[ran] = view(ẋ₊, up_dim, dim)
			accumul += len
			ran = accumul .+ (1:len)
			rows[ran] = view(lin_idx, up_dim)
			cols[ran] = view(lin_idx, up_dim)
			vals[ran] = .- view(ẋ₊, up_dim, dim)
			accumul += len

			# Going DOWN
			dwn_dim = view(dwn, :, dim)
			len = sum(dwn_dim)
			ran = accumul .+ (1:len)
			rows[ran] = view(lin_idx, dwn_dim)
			cols[ran] = view(lin_idx, dwn_dim) .- hjb.state_grids_strides[dim]
			vals[ran] = .- view(ẋ₋, dwn_dim, dim)
			accumul += len
			ran = accumul .+ (1:len)
			rows[ran] = view(lin_idx, dwn_dim)
			cols[ran] = view(lin_idx, dwn_dim)
			vals[ran] = view(ẋ₋, dwn_dim, dim)
			accumul += len
		
		end

		# Build Sparse transition matrix A (infinitesimal generator)
		A = sparse(view(rows, 1:accumul), 
				   view(cols, 1:accumul),
				   view(vals, 1:accumul),
				   hjb.nstates, 
				   hjb.nstates)

		# Compute value function update
		V_next = (sparse( (1/Δ + ρ)I, hjb.nstates, hjb.nstates) .- A) \ 
				 (R + reshape(V, hjb.nstates) ./ Δ)

		# Compute supnorm
		sup = maximum(abs.(V[:] - V_next[:]))
		update!(p, log10(sup))

		# Decision point for convergence
		if sup <= tol
			infinitesimal_generator = A
			policy_functions = c
			break
		end
	end

	if (it == maxit) & (sup > tol)
		error("Reached iteration $it/$maxit without converging.")
	end

	return HJBSolution(infinitesimal_generator, 
					   hjb.mc_state_grid[:],
					   V_next, 
					   policy_functions)
end


end # module




