% Fornino Manera (2019)
% 
% DATE: January 24, 2021
%
% Project: Automation and the Future of Work: Assessing the Role of Labor
%          Flexibility

function params = SetParameters(varargin)
% Routine that sets all parameters needed for calling LaborDemand.m. It
% creates a structure "params" with the following fields:
%
% Gamma             % Relative MRTS of Labor and Robots
% E                 % Flow Cost of Robots m
% R                 % Robot investment price pR
% P                 % Average Price of task output
% N_p               % Dimension of idiosyncratic state realization set
% N_k               % Dimension of approximating grid for Robots
% rho_P             % Persistence of Underlying AR(1)
% sigma_P           % Variance of Productivity Process
% time_scale        % Meaning of one unit time ( 1/4 for quarterly )
% t                 % DRS parameter
% delta             % Depreciation rate
% rho               % Time discount factor for firm
% psi_k             % Scale of Adjustment Cost function
% tolG              % Tolerance for Stationary Distribution
% DT                % DT for HJB Implicit Method
% print_conv        % Option to print HJB Convergence
% Utilization       % Option to activate Robot Utilization margin

% Process input
narginchk(0,2)


% Current path
current_path = cd;


if nargin > 1 && ~isempty(varargin{2})
    shockProcess =  varargin{2};
else
    shockProcess = 'Diffusion';
end

% Wage and Price of Robots
params.W = .5;
params.pR_rel = 3;
params.pR = params.pR_rel * params.W;

% Time Parameters
params.time_scale = 1;

% Approximation Parameters
params.N_k = 1000;
params.N_p = 100;

% Model parameters
params.Gamma = 0.5;
params.E = 0.2;
params.t = 0.5;
params.delta = log(1 + params.time_scale * 1 / 12);
params.rho = log(1 + params.time_scale * 0.04);
params.psi_k = 25;
params.A_prod = 1;          % productivity shifter of intermediate good producer
params.P = 1;
params.rho_P = .9;
params.theta_P = 0.1;
params.sigma_P = 0.4;
% Exogenous Shock
if strcmp(shockProcess,'Diffusion')
    params.ShockType = 'Diffusion';
    % parameters of the lognormal distribution of the stationary distribution.
    params.sigma_hat = sqrt(log(params.sigma_P.^2 ./ params.P.^2 + 1));
    params.mu_hat = log(params.P) - 1/2 .* log(params.sigma_P.^2 ./ params.P.^2 + 1);
    % definition of drift and s.d.
    params.settings.mu = @(p) - params.theta_P .* (log(p) - params.mu_hat - (params.sigma_hat).^2) .* p;
    params.settings.sigma = @(p) sqrt(2 * params.theta_P) .* params.sigma_hat .* p;
    % Exogenous shock
    opts = optimoptions('fsolve','Display','off', 'UseParallel', false);
    params.pmin = fsolve(@(p) logncdf(p, params.mu_hat, params.sigma_hat) - 1e-6, 1,opts);
    params.pmax = fsolve(@(p) 1-logncdf(p, params.mu_hat, params.sigma_hat) - 1e-6, 1,opts);
elseif strcmp(shockProcess,'GBM')
    % Approximation Parameters
    params.N_k = 100;
    params.N_p = 150;
    % calibration for representative sector
    params.ShockType = 'GBM';
    params.settings.ResetType = 'LogNormal';
    params.settings.reset_rate = 0.11302;
    params.settings.mu = @(x) 0.011091 *x;
    params.settings.sigma = @(x) 0.074446 *x; 
    params.settings.mu_reset = -0.106 ;
    params.settings.sigma_reset =  0.034357;  
    params.pmin = 0.001;
    params.pmax = 3;
    % discount factor has to be augmented with death rate
    params.rho = params.rho + params.settings.reset_rate;
    params.pR = 1.4 * params.W;
else
    error('specify second arg. as Diffusion or GBM or leave empty')
end

params.settings.GammaPGrid = 1;    

% HJB and KFE Parameters
params.tolG = 1e-10;
params.DT =1e6;
params.print_conv = 0;

% GE
params.A_robots = 0.01;    % Productivity of the Robot-Producing sector
params.alpha = 0.2;         % 1- share of capital
params.sigma_c = 2;         % IES
params.eta = .01;           % THE 1%
params.chi = 1;             % Scale Labor Supply in steady state
params.varphi = 1;         % Frisch El of Labor Supply


%     params.p_grid = exp(linspace(log(params.pmin), log(params.pmax), params.N_p)); 

% Compute Labor Savings and R_max
params.Omega = (1 - params.Gamma) ./ params.Gamma * params.W - params.E;
params.k_max = 1 / params.delta / params.psi_k * (params.Omega / (params.rho + params.delta) - params.pR);

% Reversible Investment 
params.ReversibleInvestment = 'Yes';

if nargin > 0 && ~isempty(varargin{1})

    if isnumeric(varargin{1})
        ifrCode = varargin{1};

        % Load Targets
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
        clear opts

        % *** Data on average relative price of robots in years 2004, 2007, 2010,
        %     2014. Source: see appendix.
        pR_rel04 = 2.3224;
        pR_rel07 = 1.7865;
        pR_rel10 = 1.4348;
        pR_rel14 = 1.0209;

        % *** Create Targets
        employment_target = Statistics.emp_us_ * 1000;
        vadd_target = Statistics.vadd_1989 * 1e6;
        rl_target04 = Statistics.apr_lv04 / 1000;
        rl_target07 = Statistics.apr_lv07 / 1000;
        rl_target10 = Statistics.apr_lv10 / 1000;
        rl_target14 = Statistics.apr_lv14 / 1000;
        theta_prod = Statistics.thetaProd;
        theta_p = Statistics.theta_p;
        sigma_p = Statistics.sigma_p;

        if strcmp(shockProcess, 'GBM') % import statistics
            StatisticsGBM = readtable([current_path '/data/GBMStatistics']);
        end


        % --------------------------------------------------------------- %

        % *** Calibrate Parameters

        % Wage and Maintenance Cost of Robots
        params.W = .35;
        params.E = 0;
        params.pR_rel = pR_rel10;
        params.pR = params.pR_rel * params.W;

        if strcmp(shockProcess, 'Diffusion')
            % Exogenous shock
            params.sigma_P = sigma_p(ifrCode);
            params.theta_P = theta_p(ifrCode);

            % parameters of the lognormal distribution of the stationary distribution.
            params.sigma_hat = sqrt(log(params.sigma_P.^2 ./ params.P.^2 + 1));
            params.mu_hat = log(params.P) - 1/2 .* log(params.sigma_P.^2 ./ params.P.^2 + 1);

            params.settings.mu = @(p) - params.theta_P .* (log(p) - params.mu_hat - (params.sigma_hat).^2) .* p;
            params.settings.sigma = @(p) sqrt(2 * params.theta_P) .* params.sigma_hat .* p;
            params.settings.GammaPGrid = 1;

            % Choose pmax following the variance of the stationary distribution
            opts = optimoptions('fsolve', 'Display', 'off', 'UseParallel', false);
            params.pmin = fsolve(@(p) logncdf(p, params.mu_hat, params.sigma_hat) - 1e-6, 1,opts);
            params.pmax = fsolve(@(p) 1-logncdf(p, params.mu_hat, params.sigma_hat) - 1e-6, 1,opts); 
                        % Set parameters Gamma, psi_R, A_prod
            gamma_tilde = ( params.E .* ( rl_target14(ifrCode) - rl_target10(ifrCode) ) + ...
                           ( params.rho + params.delta) * (params.W * pR_rel10 * rl_target14(ifrCode) - params.W * pR_rel14 * rl_target10(ifrCode))) ./ ...
                          ( params.W * (rl_target14(ifrCode) - rl_target10(ifrCode) ) - ...
                           ( params.rho + params.delta) .* rl_target14(ifrCode) .* rl_target10(ifrCode) .* (params.W * pR_rel10 - params.W * pR_rel14)  );
            %integral_A_prod = integral(@(p) p.^(1./(1 - params.t)) .* lognpdf(p, params.mu_hat, params.sigma_hat), 0, params.pmax);

            params.Gamma = 1 ./ (1 + gamma_tilde);
            params.psi_k = (1 + gamma_tilde * rl_target10((ifrCode))) / (employment_target((ifrCode)) * rl_target10((ifrCode)) * params.delta) * ...
                           ( (gamma_tilde * params.W - params.E) / (params.rho + params.delta) - params.W * pR_rel10 );
            params.A_prod = 1;
%             params.A_prod = (integral_A_prod.^(1 - params.t) .* params.t .* ...
%                         params.Gamma.^params.t ./ params.W ./ employment_target((ifrCode)).^(1 - params.t)).^-1;

        elseif strcmp(shockProcess, 'GBM')
            % collect all guess values
            params.settings.mu = @(x) StatisticsGBM.mu(ifrCode) .* x;
            params.settings.sigma = @(x) StatisticsGBM.sigma(ifrCode) * x;
            params.settings.reset_rate = StatisticsGBM.reset_rate(ifrCode);
            params.settings.sigma_reset = StatisticsGBM.sigma_reset(ifrCode);
            params.settings.mu_reset = StatisticsGBM.mu_reset(ifrCode);
            params.pmin = 0.001;
            params.pmax = StatisticsGBM.pmax(ifrCode);
            params.rho = params.rho + params.settings.reset_rate;
        else
        end

        % Set parameter theta
        params.t = theta_prod(ifrCode);



        % Compute Labor Savings and R_max
        params.Omega = (1 - params.Gamma) ./ params.Gamma * params.W - params.E;
        params.kmax = 1 / params.delta / params.psi_k * (params.Omega / (params.rho + params.delta) - params.pR);

        % Store Targets
        params.targets.ifrCode = ifrCode;
        params.targets.ifrString = Statistics.ifrString{ifrCode};
        params.targets.pR_rel04 = pR_rel04;
        params.targets.pR_rel07 = pR_rel07;
        params.targets.pR_rel10 = pR_rel10;
        params.targets.pR_rel14 = pR_rel14;
        params.targets.employment_target = employment_target(ifrCode);
        params.targets.vadd_target = vadd_target(ifrCode);
        params.targets.rl_target04 = rl_target04(ifrCode) * 1000;
        params.targets.rl_target07 = rl_target07(ifrCode) * 1000;
        params.targets.rl_target10 = rl_target10(ifrCode) * 1000;
        params.targets.rl_target14 = rl_target14(ifrCode) * 1000;

    else
        error('Need to specify an integer between 1 and 13 to pick the IFR Sector.')
    end
end

end

