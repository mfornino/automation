function out = generalEquilibrium(params)

% Fornino Manera (2019)
%
% This function solves for GE prices and quantities.
%
% Input ("params" struct object, as constructed from RunFigures.m)
%
% Outputs("out" struct object, same as LaborDemand_trapz.m)


    % Define function handle for fsolve to find equilibrium at given
    % Gammas.
    fun = @(xCon) mktClearingDeviations( exp(xCon), params);
    % Read the guess from the disk to save on computation time.
    try
        x0 = readmatrix('x0.csv');
%         x0 = x0(1:14);
    catch
        x0 = [.1*ones(13,1); 2.7e-5; 1];
    end
    opts = optimoptions('fsolve',...
                        'display', 'off', ...
                        'useparallel', false, ...
                        'MaxFunctionEvaluations', 15 * 10000,...
                        'MaxIterations', 1e5,...
                        'algorithm', 'trust-region-dogleg');
    % Solve for the equilibrium
    [x, ~, FLAG] = fsolve(fun,log(x0),opts);
    x = exp(x);

    % Return NaN if the solver did not converge
    if FLAG < 1
        
        out = NaN;
    
    % Else return the deviations from the R/L targets
    else
            
        % Re-Compute equilibrium at solution to obtain the R/L in the model
        [~, out] = mktClearingDeviations(x, params);

    end
end


%% FUNCTION TO COMPUTE MARKET CLEARING
function [eq, out_eqm] = mktClearingDeviations(variables, params)

% Fornino Manera (2019)
%
% This subroutine computes the deviations from Market Clearing given the
% prices. Used in the call to fsolve in the main routine at the top.
%
% Input
% - variables: vector containing the prices and the shifter A_F
% - params: as constructed in RunFigures.m
%
% Outputs
% - eq: vector of deviations from Market Clearing 
% - out_eqm: structure as constructed by LaborDemand_trapz.m

    % Convention is that variables are the prices + value for AF
    prices = variables(1:14);
%     AF = params{1}.AF;
    
    chi = params{1}.chi;
    varphi = params{1}.varphi;
    
%     AF = 1;
    AF = variables(15);
    
    % Preallocate useful objects
    out_eqm = cell(13,1);
    labordemand_sec = zeros(13,1);
    output_sec = labordemand_sec;
    income_sec = labordemand_sec;
    xi_vec = labordemand_sec;
    % Cycle through the 13 sectors and compute the solution of the model
    % for given prices and parameters.
    parfor ii = 1:13
        
        xi_vec(ii) = params{ii}.xi;
        
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
        [labordemand_sec(ii), out] = LaborDemand_trapz(params_eqm.W, params_eqm.pR, params_eqm);
        

        if isfield(params_eqm, 'A_sec_mrts')
            output_sec(ii) = out.GE.q * params_eqm.A_sec_mrts;
        else
            output_sec(ii) = out.GE.q;
        end
        income_sec(ii) = prices(ii) * output_sec(ii);
        out_eqm{ii, 1} = out;
        
        % Substitute Robots per 1000 Employees with Robot per Employee
        out_eqm{ii, 1}.RL = out_eqm{ii, 1}.RL / 1000;
    end
  
    % MARKET CLEARING DEVIATION FOR INTERMEDIATE GOODS
    eq(1:13,1) = (output_sec - xi_vec .* sum(income_sec) ./ prices(1:13,1) ./ AF) ./ output_sec;
    
    % MARKET CLEARING DEVIATION FOR LABOR
    eq(14,1) = (sum(labordemand_sec) - (prices(end)./chi) .^ varphi) ./ sum(labordemand_sec);

    % NUMERAIRE (NOT USED)
    eq(15) = prod((variables(1:13) ./xi_vec).^ xi_vec) - 1;
end