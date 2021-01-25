% Fornino Manera (2019)
% 
% DATE: January 24, 2021
%
% Project: Automation and the Future of Work: Assessing the Role of Labor
%          Flexibility

function out_eqm = partialEquilibrium(variables,params)
% Routine that computes partial equilibrium outcomes for given prices.
%
% Input ("settings" struct object)
% - variables: vector with prices
% - params: struct with parameters
%
% Outputs("out" struct object)
% - out_eqm: cell object with outputs of LaborDemand_traps.m

% Convention is that variables are the prices + value for AF
prices = variables(1:14);

% Preallocate useful objects
out_eqm = cell(13,1);


% Cycle through the 13 sectors and compute the solution of the model
% for given prices and parameters.
parfor ii = 1:13

    % Trick to construct appropriate params structure.
    params_eqm = params{ii};

    if isfield(params_eqm, 'A_sec_mrts')
        params_eqm.A_prod = prices(ii) * params_eqm.A_sec_mrts;
    else
        % Fix the price of the good in sector ii
        params_eqm.A_prod = prices(ii);         % NOTE: A_prod acts as a price
    end

    % Set Prices
    params_eqm.W = prices(end);

    % Compute R_max^star to set a suitable grid for the Robot input
    kmax = 1 ./ params_eqm.psi_k ./ params_eqm.delta * (( params_eqm.W * ...
        (1 - params_eqm.Gamma) / params_eqm.Gamma - params_eqm.E) / ...
        (params_eqm.rho + params_eqm.delta) - params_eqm.pR );

    % Takes care of the case with negative PDV of robots
    if kmax > 0
        params_eqm.kmax = kmax;
    end

    % Run Labor Demand
    [~, out] = LaborDemand_trapz(params_eqm.W, params_eqm.pR, params_eqm);

    out.RL = out.RL / 1e3;

    out_eqm{ii, 1} = out;
end


end