% Fornino Manera (2019)
% 
% DATE: January 24, 2021
%
% Project: Automation and the Future of Work: Assessing the Role of Labor
%          Flexibility

function params = SetParametersGE(pR_rel_pre, pR_rel_post, varargin)
% This function is called by RunFigures.m to calibrate the model in General
% Equilibrium. The objective is to match the aggregate change in robot 
% penetration between 2010 and 2014, as well as the sectoral robot
% penetration in 2014.
%
% Input:
% - variable: check below
%
% Outputs
% - params structure for the GE case

%% Current path
current_path = cd;

%% Check inputs
narginchk(2,5)

%% Load Targets
% *** Import variables to create calibration targets
opts = delimitedTextImportOptions("NumVariables", 12);
opts.DataLines = [2, Inf];
opts.Delimiter = ",";
opts.VariableNames = ["ifrCode", ...
                      "ifrString", ...
                      "apr_lv04", ...
                      "apr_lv07", ...
                      "apr_lv10", ...
                      "apr_lv14", ...
                      "thetaProd", ...
                      "theta_p", ...
                      "sigma_p", ...
                      "emp_us_", ...
                      "vadd_1989", ...
                      "shares"];
opts.VariableTypes = ["string", ...
                      "string", ...
                      "double", ...
                      "double", ...
                      "double", ...
                      "double", ...
                      "double", ...
                      "double", ...
                      "double", ...
                      "double", ...
                      "double", ...
                      "double"];
opts = setvaropts(opts, 1, "WhitespaceRule", "preserve");
opts = setvaropts(opts, 1, "EmptyFieldRule", "auto");
opts.ExtraColumnsRule = "ignore";
opts.EmptyLineRule = "read";
Statistics = readtable([current_path '/data/OUStatistics.csv'], opts);
Statistics.shares = Statistics.shares / sum(Statistics.shares);
clear opts

% *** Data on average relative price of robots in years 2004, 2007, 2010,
%     2014. Source: see appendix.

pR_rel10 = pR_rel_pre;
pR_rel14 = pR_rel_post;

% *** Create Targets
employment_target = Statistics.emp_us_ * 1000;
vadd_target = Statistics.vadd_1989 * 1e6;
rl_target04 = Statistics.apr_lv04 / 1000;
rl_target07 = Statistics.apr_lv07 / 1000;
rl_target10 = Statistics.apr_lv10 / 1000;
rl_target14 = Statistics.apr_lv14 / 1000;
theta_prod = Statistics.thetaProd;
% theta_p = Statistics.theta_p_detrend;
% sigma_p = Statistics.sigma_p_detrend;
% theta_p = Statistics.theta_p_HPF;
% sigma_p = Statistics.sigma_p_HPF;
theta_p = Statistics.theta_p;
sigma_p = Statistics.sigma_p;
shares = Statistics.shares;

% Set theta_i and xi_i
params_cal.xi = shares;
params_cal.theta_P = theta_p;
params_cal.sigma_P = sigma_p;
params_cal.theta = theta_prod;

params_cal.m_rel = 0;
params_cal.rho = log(1 + .04);
params_cal.delta = log(1 + 1/12);
params_cal.L_target = 1;
params_cal.rl_target = rl_target14;
params_cal.rl_target_previous = rl_target10;
params_cal.pR_rel = pR_rel14;
params_cal.pR_rel_previous = pR_rel10;
params_cal.Utilization = 'Yes';

params_cal.N_k = 100;
params_cal.N_p = 25;



% Stochastic process switch
if nargin > 4 && ~isempty(varargin{3})
    params_cal.ShockType = varargin{3};
else
    params_cal.ShockType = 'Diffusion';
end


%% FIND GUESSES FOR GAMMA AND PSI_R_REL

if nargin > 2 && ~isempty(varargin{1})
    
    disp('STEP 0: VALUE OF \psi_R PROVIDED BY USER... SKIPPING')
    
    x = varargin{1};
    
else
    disp('STEP 0: FIND VALUE OF \psi_R TO MATCH CHANGE IN AGGREGATE ROBOT PENETRATION')

    opts = optimset('display','off', 'TolX', 1e-6);
%     opts = optimoptions('fsolve',...
%                         'display','iter');
    fun = @(log_x) deltarlDeviation( exp(log_x), params_cal);
    log_x0 = log(8900);
    if strcmp(params_cal.ShockType, 'GBM')
        log_x0 = 9.5389;
    end

    log_x = fzero(fun, log_x0,opts);
%     log_x = fsolve(fun, log_x0, opts);
    x = exp(log_x);
end

params_cal.psi_R_rel = x;

if ~isreal(params_cal.psi_R_rel)
    error('Imaginary number for psi_R_rel')
end

gamma_upperbar = 1 ./ (1 + params_cal.m_rel + ...
    (params_cal.rho + params_cal.delta) .* params_cal.pR_rel);

R_max_i_guess = params_cal.rl_target .* params_cal.xi;

gamma_guess = 1 ./ (1 + params_cal.m_rel + ...
        (params_cal.rho + params_cal.delta) .*(params_cal.delta .* ...
        params_cal.psi_R_rel .* R_max_i_guess + params_cal.pR_rel));
    
if strcmp(params_cal.ShockType, 'GBM') % discount factor varies with ii
    StatisticsGBM = readtable([current_path '/data/GBMStatistics']);
    gamma_upperbar = 1 ./ (1 + params_cal.m_rel + ...
        (params_cal.rho + params_cal.delta + StatisticsGBM.reset_rate ...
        ) .* params_cal.pR_rel);
    gamma_guess = 1 ./ (1 + params_cal.m_rel + ...
        (params_cal.rho + params_cal.delta + StatisticsGBM.reset_rate) .* ...
        (params_cal.delta .* ...
        params_cal.psi_R_rel .* R_max_i_guess + params_cal.pR_rel));
    
end
    
    
% gamma_guess = min(gamma_guess, ones(size(gamma_guess)) .* gamma_upperbar *.99);
    
psi_R_rel_base = params_cal.psi_R_rel;


disp('DONE! ')
disp(' ')

%% RUN NESTED CALIBRATION EXERCISE IN GENERAL EQUILIBRIUM

if nargin > 3 && ~isempty(varargin{2})
    
    disp('STEP 1: VALUES OF \Gamma_i PROVIDED BY USER... SKIPPING')
    
    x = varargin{2};
    
else


    disp('STEP 1: FIND VALUES OF \Gamma_i MATCHING THE ROBOT PENETRATION ACROSS SECTORS IN 2014')

    % Function handle to compute the deviations from R/L targets in Gen Eqm.
    fun = @(x) rlDeviations( gamma_upperbar ./ (1 + exp(-x)), params_cal);

    % Read the guess from the disk to save on computation time.
    if strcmp(params_cal.ShockType, 'GBM')
        mat = 'x0_gamma_GBM.csv';
    else
        mat = 'x0_gamma.csv';
    end
    try 
        x0 = readmatrix([current_path '/guesses/' mat]);
    catch
        x0 = gamma_guess;
    end
    
    if ~isempty(find((x0 - gamma_upperbar) > 0, 1))
        error('guess exceeds gamma_upperbar (%0.2f)', gamma_upperbar);
    end
%     x0 = min(gamma_guess, ones(size(gamma_guess)) * gamma_upperbar *.99);
    
    opts = optimoptions('fsolve',...
                        'display', 'off', ...
                        'useparallel', false, ...
                        'MaxFunctionEvaluations', 13 * 10000,...
                        'MaxIterations', 1e5,...
                        'algorithm', 'trust-region-dogleg');
    [x, ~, FLAG] = fsolve(fun,log(x0./ (gamma_upperbar -x0)),opts);
    x = gamma_upperbar./(1 + exp(-x));

    % Write guess to disk if successful for future use.
    % if FLAG >= 1
    %     writematrix(x, [current_path '/guesses/' mat])
    % end
end

disp('DONE!')
disp(' ')

%% OBTAIN EQUILIBRIUM OBJECTS AT CALIBRATED VALUES
    
disp('STEP 2: GATHER FINAL CALIBRATION IN GENERAL EQUILIBRIUM.')

[~, out] = rlDeviations(x, params_cal);


% Unpack prices and wage
prices = zeros(13,1);
labor = prices;
rl = prices;
q = prices;
for ii = 1:13
    prices(ii) = out{ii}.params.A_prod;
    labor(ii) = out{ii}.L;
    rl(ii) = out{ii}.RL;
    q(ii) = out{ii}.GE.q;
end
wage = out{ii}.params.W;



%% POPULATE params STRUCTURE FOR OUTPUT

params = cell(13,1);

for ii = 1:13
    params{ii} = out{ii}.params;
    params{ii}.AF = out{ii}.params.AF;
    params{ii}.xi = params_cal.xi(ii);
    params{ii}.ReducedOutput = false;
    params{ii}.varphi = 1; % FRISCH ELASTICITY OF L.S.
    params{ii}.chi = wage ./ ( params_cal.L_target .^ (1/params{ii}.varphi) );    
    params{ii}.outCalibration = out{ii};
end

disp('DONE!')

end


%% FUNCTION TO FIND BALLPARK VALUE FOR PSI_R_REL

function eq = deltarlDeviation(psi_R_rel, params_cal)

    params_cal.psi_R_rel = psi_R_rel;

    R_max_i_guess = params_cal.rl_target .* params_cal.xi;

    gamma_guess = 1 ./ (1 + params_cal.m_rel + ...
        (params_cal.rho + params_cal.delta) .*(params_cal.delta .* ...
        psi_R_rel .* R_max_i_guess + params_cal.pR_rel));
%     
%     if strcmp(params_cal.ShockType, 'GBM') % discount factor varies with ii
%         StatisticsGBM = readtable('GBMStatistics');
%         gamma_upperbar = 1 ./ (1 + params_cal.m_rel + ...
%             (params_cal.rho + params_cal.delta + StatisticsGBM.reset_rate ...
%             ) .* params_cal.pR_rel);
%         gamma_guess = 1 ./ (1 + params_cal.m_rel + ...
%             (params_cal.rho + params_cal.delta + StatisticsGBM.reset_rate) .* ...
%             (params_cal.delta .* ...
%             params_cal.psi_R_rel .* R_max_i_guess + params_cal.pR_rel));
%         gamma_guess = min(gamma_guess, .9999* gamma_upperbar);
%     
%     end

    % OBTAIN EQUILIBRIUM OBJECTS AT CALIBRATED VALUES
    params_previous = params_cal;
    params_previous.pR_rel = params_previous.pR_rel_previous;

       
    [~, out_previous] = rlDeviations(gamma_guess, params_previous);
    [~, out] = rlDeviations(gamma_guess, params_cal);

    rl_change_model = 0;
    rl_change_data = 0;
    rl_base_model = 0;
    rl_base_data = 0;

    for ii = 1:13
        rl_change_data = rl_change_data + params_cal.xi(ii) .* (  ...
            params_cal.rl_target(ii) - params_cal.rl_target_previous(ii) );
        rl_change_model = rl_change_model + params_cal.xi(ii) .* (  ...
             (out{ii}.RL - out_previous{ii}.RL) );
        rl_base_model = rl_base_model + params_cal.xi(ii) .* out_previous{ii}.RL;
        rl_base_data = rl_base_data + params_cal.xi(ii) .* params_cal.rl_target_previous(ii); 
    end
    
    eq = rl_change_data/rl_base_data - rl_change_model/rl_base_model;
    
end




%% FUNCTION TO COMPUTE DEVIATIONS FROM R/L IN GENERAL EQUILIBRIUM

function [eq, out] = rlDeviations(gamma,params_cal)

    % Current path
    current_path = cd;

    % Impose provisional values of Gamma for each sector
    params_cal.Gamma = gamma;

    % Define function handle for fsolve to find equilibrium at given
    % Gammas.
    fun = @(xCon) mktClearingDeviations_cal( exp(xCon), params_cal);
    % Read the guess from the disk to save on computation time.
    try
        x0 = readmatrix([current_path '/guesses/x0.csv']);
        if ~isreal(x0)
          x0 = [.1*ones(13,1); .01; 1]; 
        end
    catch
        x0 = [.1*ones(13,1); .01; 1];
    end
%     opts = optimoptions('fsolve',...
%                         'display', 'iter', ...
%                         'useparallel', false, ...
%                         'MaxFunctionEvaluations', 15 * 10000,...
%                         'MaxIterations', 1e5,...
%                         'algorithm', 'trust-region-dogleg');
    opts = optimoptions('fsolve',...
                        'display', 'off', ...
                        'useparallel', false, ...
                        'algorithm', 'trust-region-dogleg');
    % Solve for the equilibrium
    [x, ~, FLAG] = fsolve(fun,log(x0),opts);
    x = exp(x);

    % Return NaN if the solver did not converge
    if FLAG < 1
        
        if nargout > 1
            out = NaN;
        end
        eq = NaN(13,1);
    
    % Else return the deviations from the R/L targets
    else
        
%         disp('Equilibrium Found!')
        % Write the equilibrium for current guess to disk to save on
        % computation time.
         writematrix(x, [current_path '/guesses/x0.csv']);
        
        % Re-Compute equilibrium at solution to obtain the R/L in the model
        [~, out] = mktClearingDeviations_cal(x, params_cal);

        % Compute deviations from R/L targets (in %)
        eq = zeros(13,1);
        for ii = 1:13
            eq(ii,1) =  (out{ii}.RL - params_cal.rl_target(ii))./ params_cal.rl_target(ii);
        end
    end
end


%% FUNCTION TO COMPUTE DEVIATIONS FROM MARKET CLEARING AND NUMERAIRE

function [eq, out_eqm] = mktClearingDeviations_cal(variables, params_cal)

    % Convention is that variables are the prices + value for AF
    prices = variables(1:14);
    AF = variables(15);
%     AF = 1;
    
    % Preallocate useful objects
    out_eqm = cell(13,1);
    labordemand_sec = zeros(13,1);
    output_sec = labordemand_sec;
    income_sec = labordemand_sec;
    
    % Cycle through the 13 sectors and compute the solution of the model
    % for given prices and parameters.
    parfor ii = 1:13
        
        % Trick to construct appropriate params structure.
        params = SetParameters(ii, params_cal.ShockType);
        
        params.ReducedOutput = true;
        
        params.AF = AF;
        
        % Fix the parameters to needed values
        params.N_k = params_cal.N_k;
        params.N_p = params_cal.N_p;
        params.A_prod = prices(ii);         % NOTE: A_prod acts as a price
        params.Gamma = params_cal.Gamma(ii);
        params.psi_k = params_cal.psi_R_rel * prices(end);
%         params.psi_k = params_cal.psi_R_rel;
        params.t = params_cal.theta(ii);
        
        % Set Prices
        params.W = prices(end);
        params.pR = params_cal.pR_rel * prices(end);
        params.E = params_cal.m_rel * prices(end);
        
        params = rmfield(params, 'k_max');
        params = rmfield(params, 'Omega');
        
        kmax = 1 ./ params.psi_k ./ params.delta * (( params.W * ...
            (1 - params.Gamma) / params.Gamma - params.E) / ...
            (params.rho + params.delta) - params.pR );

        if kmax > 0
            params.kmax = kmax;
        end
        
        % Run Labor Demand
        [labordemand_sec(ii), out] = LaborDemand_trapz(params.W, params.pR, params);
        output_sec(ii) = out.GE.q;
%         disp(['sector' num2str(ii)])
        income_sec(ii) = prices(ii)*out.GE.q ;
        
        out_eqm{ii, 1} = out;
        
        out_eqm{ii, 1}.RL = out_eqm{ii, 1}.RL ./ 1e3;
        
        % Save value of AF
        out_eqm{ii,1}.params.AF = AF;
        
    end
  
    % IMPOSE MARKET CLEARING FOR INTERMEDIATE GOODS
    eq(1:13,1) = (output_sec - params_cal.xi(1:13,1) .* sum(income_sec) ./ prices(1:13,1) ./ AF) ./ output_sec;
    
    % IMPOSE MARKET CLEARING FOR LABOR
    eq(14,1) = (sum(labordemand_sec) - params_cal.L_target) ./ params_cal.L_target;

    % IMPOSE THE NUMERAIRE
    eq(15) = prod((variables(1:13) ./params_cal.xi).^ params_cal.xi) - 1;

    

end

%}


