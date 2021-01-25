% Fornino Manera (2019)
% 
% DATE: January 24, 2021
%
% Project: Automation and the Future of Work: Assessing the Role of Labor
%          Flexibility

function out = generalEquilibrium_onesector(params, calibrate_chi)
% This function solves for GE prices and quantities for the case of one
% representative sector.
%
% Input ("params" struct object, as constructed from RunFigures.m)
%
% Outputs("out" struct object, same as LaborDemand_trapz.m)


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% VERY IMPORTANT: WE HAVE HARD CODED L_target here

L_target = 1;

%%%%%%%%%%%%%%%%%%%%%%%%%%

% Define function handle for fsolve to find equilibrium at given
% Gammas.
fun = @(xCon) mktClearingDeviations( exp(xCon), params, calibrate_chi, L_target);
% Read the guess from the disk to save on computation time.
%     try
%         x0 = readmatrix('x0.csv');
% %         x0 = x0(1:14);
%     catch
    x0 = .4298;
%     end
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
    if calibrate_chi
        params.chi = x ./ ( L_target .^ (1/params.varphi) );
    end
    [~, out] = mktClearingDeviations(x, params, false, L_target);
end
end


%% FUNCTION TO COMPUTE MARKET CLEARING
function [eq, out_eqm] = mktClearingDeviations(wage, params, calibrate_chi, L_target)
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

% Labor Supply parameters
chi = params.chi;
varphi = params.varphi;

% Trick to construct appropriate params structure.
params_eqm = params;

% Here there is only one good market
params_eqm.A_prod = 1;

% Set Prices
params_eqm.W = wage;
%     if relative_pR
    params_eqm.pR = params.pR .* wage;
%     end

% Run Labor Demand
[labordemand_sec, out] = LaborDemand_trapz(params_eqm.W, params_eqm.pR, params_eqm);

out_eqm = out;

% Substitute Robots per 1000 Employees with Robot per Employee
out_eqm.RL = out_eqm.RL / 1000;

% MARKET CLEARING DEVIATION FOR LABOR
if calibrate_chi
    eq(1) = (labordemand_sec - L_target) ./ labordemand_sec;
else
    eq(1) = (labordemand_sec - (wage ./ chi) .^ varphi) ./ labordemand_sec;
end
end