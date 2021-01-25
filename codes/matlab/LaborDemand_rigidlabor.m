% Fornino Manera (2019)
% 
% DATE: January 24, 2021
%
% Project: Automation and the Future of Work: Assessing the Role of Labor
%          Flexibility

function [L, out] = LaborDemand_rigidlabor(W, R, params)
% This function finds the solution to the firm problem as well as the
% stationary distribution for given parameters. The model is the one with
% labor adjustment costs instead of perfectly flexible labor.
%
% Input ("settings" struct object)
% - W: wage for workers
% - R: pR, or purchase price of robots
% - params: struct object created by SetParameters()
%
% Outputs("out" struct object)
% - L: aggregate sectoral labor demand at the wage W
% - out: struct object with a host of information on the model solution.

opts = optimset('display','off');

out.params = params;
out.params.pR = R;
out.params.pR_rel = R/W;
out.params.W = W;

%% Shock Type
if isfield(params,'ShockType')
    ShockType = params.ShockType;
else
    ShockType = 'Lognormal';
    out.params.ShockType = ShockType;
end


%% Choice for Putty-Clay v. Utilization
if isfield(params,'Utilization')
    Utilization = params.Utilization;
    if ~strcmp(Utilization, 'Yes')
        warning('Current code only for Robot Utilization = Yes, resetting')
        Utilization = 'Yes';
        out.params.Utilization = Utilization;
    end
else
    Utilization = 'Yes';
    out.params.Utilization = Utilization;
end

if isfield(params,'LaborUtilization')
    LaborUtilization = params.LaborUtilization;
else
    LaborUtilization = 'No';
    out.params.LaborUtilization = LaborUtilization;
end

%% Irreversibility v. costly reversibility
if isfield(params,'ReversibleInvestment')
    ReversibleInvestment = params.ReversibleInvestment;
else
    ReversibleInvestment = 'No';
    out.params.ReversibleInvestment = ReversibleInvestment;
end


%% Modelling choices for GE:
if isfield(params,'RobotsProducers')
    %Use either 'Workers', 'Capitalists', 'None', def to none
    RobotsProducers = params.RobotsProducers;
    if ~strcmp(params.RobotsProducers, 'None')
        if isfield(params,'A_robots')
            A_robots = params.A_robots;
            R = W/A_robots; %(see end of code for explanation of prod fun)
            out.params.R = R;
        else
            error('Specify robot producers productivity A_robots in params')
        end
    end     
else
    RobotsProducers = 'None';
    out.params.RobotsProducers = RobotsProducers;
end

if isfield(params,'InvGoodCost')
    %Use 'Final', 'External' or 'Intermediate', latter is default
    %This means: produced by final good sector, produced by a mysterious
    %external sector, or SELF-produced by the intermediate goods firm
    %internally
    InvGoodCost = params.InvGoodCost;
    if strcmp(InvGoodCost, 'Final')
        if isfield(params,'p_y')
            p_y = params.p_y;
        else
            if isfield(params,'r')
                r = params.r;
                alpha = params.alpha;
                p_y = r^alpha / ( (1-alpha) ^ (1-alpha) * alpha ^ alpha);
                out.params.p_y = p_y;
            else
                error('either p_y or r have to provided for investment costs in terms of the final good!')
            end
        end
    end           
else
    InvGoodCost = 'Intermediate';
    out.params.InvGoodCost = InvGoodCost;
end


if isfield(params, 'EnergyCosts')
    %Use 'Final', 'External' or 'Intermediate', latter is default
    %This means: produced by final good sector, produced by a mysterious
    %external sector, or SELF-produced by the intermediate goods firm
    %internally
    EnergyCosts = params.EnergyCosts;
else
    EnergyCosts = 'Final';
    out.params.EnergyCosts = 'Final';
end

%% Set  Parameter Values: 

if isfield(params, 'A_prod')
    A_prod = params.A_prod;
else
    A_prod = 1;
    out.params.A_prod = 1;
end
    
% set some value for R, fundamental for solution of unconstrained problem

% Price Process Parameters
P = params.P; % average price
N_p = params.N_p; %number of grid points
rho_P = params.rho_P; %autoregressive coefficient
sigma_P = params.sigma_P; %price standard deviation
sigma_p_eps = sqrt((1 - rho_P^2) * log(sigma_P / P^2 + 1)); % standard deviation of the price error
% sigma_p_eps = sqrt( log(sigma_P / P^2 + 1)); % standard deviation of the price error
mu_p = (1 - rho_P) * (log(P) - 0.5 * (sigma_p_eps^2 / (1 - rho_P^2))); % mean of the process
% mu_p = (log(P) - 0.5 * sigma_p_eps^2); % mean of the process
time_scale = params.time_scale;
Gamma = params.Gamma;
E = params.E;
% decreasing return parameter (theta in my notation)
t = params.t;
% depreciation
delta = params.delta;
sep = params.sep;
% discount rate 
rho = params.rho;
% Adjustment cost parameters
psi_k = params.psi_k; % capital
psi_l = params.psi_l; % labor
% Volatility
if isfield(params, 'volatility')
    volatility = params.volatility;
else
    volatility = 1;
    out.params.volatility = volatility;
end

% shrink k-grid? 1 if no
if isfield(params, 'scale_k')
    scale_k = params.scale_k;
else
    scale_k =1;
end

if isfield(params, 'print_conv')
    print_conv = params.print_conv;
else
    print_conv = 0;
    out.params.print_conv = 0;
end

if (1-Gamma)/Gamma * W - E < 0
    error('Change E, W to get positive labor savings!')
end
%% Cost options
% Pay for delta?

if isfield(params, 'pay_for_delta')
    pay_for_delta = params.pay_for_delta;
else
    pay_for_delta = 1;
end

% Stochastic?
if isfield(params, 'stochastic')
    stochastic = params.stochastic;
else
    stochastic = 1;
end

% Investment in percent? (NOTE: need to exclude alpha,k ==0 if in percent)
if isfield(params, 'inv_in_percent')
    inv_in_percent = params.stochastic.inv_in_percent;
else
    inv_in_percent = 1;
end
%% Functions

%OLD GAMMA

% % Operating Profit function; WITHOUT INVESTMENT-RELATED COSTS
% k_bar = @(P) ( (P * A_prod * Gamma * t) / W ).^(1/(1 - t) );
% 
% Pi = @ (k , P) (k <= k_bar(P)) .* (...
%     (1 - t) * A_prod * P  .* (k_bar(P)).^(t) + (W / Gamma - E) * k ) + ...
%     (k > k_bar(P)) .* (...
%     A_prod * P .* k .^t - E * k);
% 
% % Optimal unconstrained capital (FOC of above)
% K_star = @(P) max( (t * A_prod * P/ E) .^(1/(1 - t)), (Gamma * t * A_prod * P/ W) .^(1/(1 - t)));
% %K_star = @(P) 2 * k_bar(P);

%NEW GAMMA
% Operating Profit function; WITHOUT INVESTMENT-RELATED COSTS

if strcmp(Utilization,'No')
    k_bar = @(P) 1/(1-Gamma)*( W./(A_prod*P*Gamma*t)).^(1/(t-1)) ;

    Pi = @ (k ,l, P) ((Gamma * l + (1-Gamma) * k) <= (1-Gamma)*k_bar(P)).* (...
         A_prod * P (Gamma * l + (1-Gamma) * k).^t - W * l - E * k ) +...
         ((Gamma * l + (1-Gamma) * k) > (1-Gamma) * k_bar(P) ).* (...
    (k <= k_bar(P)) .* (...
        (1 - t) * (1 - Gamma)^t * A_prod * P  .* (k_bar(P)).^(t) ...
        + (W * (1-Gamma) / Gamma - E) * k ) + ...
        (k > k_bar(P)) .* (...
        A_prod * (1 - Gamma)^t * P .* k .^t - E * k) );
else
   
    if ~strcmp(LaborUtilization,'No')

        % AMENDED PROFIT AND KBAR if there is a utilization margin
        k_bar = @(P) 1/(1-Gamma)*( W./(A_prod*P*Gamma*t)).^(1/(t-1)) ;
        k_hat = @(P) ( (A_prod*P*(1-Gamma)^t*t)./E).^(1/(1-t)) ;


        p_hat = @(R) E / (A_prod * t * (1 - Gamma)^t) .* R .^ (1 - t);
        p_bar = @(R) W / (A_prod * t * Gamma) .* (R .* (1 - Gamma)).^ (1 - t);


        Pi = @ (k ,l, P) ((Gamma * l + (1-Gamma) * k) <= (1-Gamma)*k_bar(P)).* (...
             A_prod .* P .* (Gamma * l + (1-Gamma) * k).^t - W * l - E * k ) +...
             ((Gamma * l + (1-Gamma) * k) > (1-Gamma) * k_bar(P) ).* (...
        (k <= k_bar(P)) .* (...
            (1 - t) * (1 - Gamma)^t * A_prod * P  .* (k_bar(P)).^(t) ...
            + (W * (1-Gamma) / Gamma - E) * k ) + ...
            (k > k_bar(P)) .* (k < k_hat(P)) .* (...
            A_prod * (1 - Gamma)^t * P .* k .^t - E * k) +...
            (k >= k_hat(P)) .* (...
            A_prod * (1 - Gamma)^t * P .* k_hat(P) .^t - E * k_hat(P)));

        % derivative wrt k
        Pi_prime = @ (k ,l, P) ((Gamma * l + (1-Gamma) * k) <= (1-Gamma)*k_bar(P)).* (...
             A_prod .* P .* t .* (1-Gamma) * (Gamma * l + (1-Gamma) * k).^(t - 1)  - E  ) +...
             ((Gamma * l + (1-Gamma) * k) > (1-Gamma) * k_bar(P) ).* (...
        (k <= k_bar(P)) .* (...
             (W * (1-Gamma) / Gamma - E)  ) + ...
            (k > k_bar(P)) .* (k < k_hat(P)) .* (...
            A_prod * t (1 - Gamma)^t * P .* k .^(t - 1) - E ) );

        Pi_prime_l = @ (k ,l, P) ((Gamma * l + (1-Gamma) * k) <= (1-Gamma)*k_bar(P)).* (...
             A_prod * P * t * Gamma * (Gamma * l + (1-Gamma) * k).^(t - 1)  - W  );


        if E==0
            % Here k_hat(P) is Inf, so we are never in the highest region.
            % Labor savings cannot turn negative here!
            % Utilization always set to u=1 because of no flow cost.
            Pi = @ (k ,l, P) ((Gamma * l + (1-Gamma) * k) <= (1-Gamma)*k_bar(P)).* (...
             A_prod * P .* (Gamma * l + (1-Gamma) * k).^t - W * l - E * k ) +...
             ((Gamma * l + (1-Gamma) * k) > (1-Gamma) * k_bar(P) ).* (...
            (k <= k_bar(P)) .* (...
            (1 - t) * (1 - Gamma)^t * A_prod * P  .* (k_bar(P)).^(t) ...
            + (W * (1-Gamma) / Gamma - E) * k ) + ...
            (k > k_bar(P)) .* (...
            A_prod * (1 - Gamma)^t * P .* k .^t - E * k) );

        end
    else
        % AMENDED PROFIT AND KBAR if there is a utilization margin
        k_bar = @(P) 1/(1-Gamma)*( W./(A_prod*P*Gamma*t)).^(1/(t-1)) ;
        k_hat = @(P) ( (A_prod*P*(1-Gamma)^t*t)./E).^(1/(1-t)) ;


        p_hat = @(R) E / (A_prod * t * (1 - Gamma)^t) .* R .^ (1 - t);
        p_bar = @(R) W / (A_prod * t * Gamma) .* (R .* (1 - Gamma)).^ (1 - t);

        Pi = @ (k ,l, P) ((Gamma * l) < (1-Gamma)*(k_hat(P) - k)).* (...
             A_prod .* P .* (Gamma * l + (1-Gamma) * k).^t - W * l - E * k ) +...
             ((Gamma * l) >= (1-Gamma)*(k_hat(P) - k)).*... 
         ((Gamma * l)<= (1-Gamma)*(k_hat(P) )).* (...
            A_prod .* P .* ((1-Gamma) * k_hat(P)).^t - (W - Gamma/(1-Gamma) * E) * l...
            -  E * k_hat(P)) + ...
            ((Gamma * l) > (1-Gamma)*(k_hat(P) )) .* (...
             A_prod .* P .* (Gamma * l ).^t - W * l);
         
         Pi_prime = @ (k ,l, P) ((Gamma * l) < (1-Gamma)*(k_hat(P) - k)).* (...
             A_prod .* P .* t .* (1-Gamma) * (Gamma * l + (1-Gamma) * k).^(t - 1)  - E  );
        
    end
        if E==0
            % Here k_hat(P) is Inf, so we are never in the highest region.
            % Labor savings cannot turn negative here!
            % Utilization always set to u=1 because of no flow cost.
            Pi = @ (k ,l, P)... 
             A_prod .* P .* (Gamma * l + (1-Gamma) * k).^t - W * l - E * k;
        end

end
    






% Adjustment cost functions
if pay_for_delta == 0
    phi_k = @(I_k) psi_k/2 .* (I_k) .^2;
%     phi_k = @(I_k) psi_k/2 .* (I_k - delta) .^2;
else
    phi_k = @(I_k) psi_k/2 .* (I_k) .^2;
end
% Derivatives

if pay_for_delta == 0
    dphi_k = @(I_k) psi_k .* (I_k) ;
%     dphi_k = @(I_k) psi_k .* (I_k - delta) ;
else
    dphi_k = @(I_k) psi_k .* (I_k) ;
end
% Inverse derivatives
if pay_for_delta == 0
    inv_dphi_k = @(y) 1/psi_k .* y;
%     inv_dphi_k = @(y) 1/psi_k .* y + delta ;
else
    inv_dphi_k = @(y) 1/psi_k .* y;
end


% Optimal unconstrained capital (FOC of above)
% K_star = @(P) max( (t * (1 - Gamma)^t * A_prod * P/ E) .^(1/(1 - t)), k_bar(P));

if E==0
    K_star_guess = @(P) k_bar(P);
else
    K_star_guess = @(P) (t * (1 - Gamma)^t * A_prod * P/ E) .^(1/(1 - t));
end


if pay_for_delta == 0
   solve_max_k = @(k,P)  t * A_prod * (1 - Gamma)^t * P .* k .^(t-1) - E - delta * R;
else 
   solve_max_k = @(k,P) -delta * dphi_k(delta * k) + t * A_prod * (1 - Gamma)^t * P .* k .^(t-1) - E - delta * R; 
end

K_star = @(P) fsolve(@(k) solve_max_k(k,P), K_star_guess(P), opts);
% K_star = @(P) fsolve(@(k) solve_max_k(k,P), 1, opts);
% K_star = @(P) K_star_guess(P);

% DEBUG
% opts_deb = optimoptions('fsolve', 'display', 'iter');
% K_star = @(P) fsolve(@(k) solve_max_k(k,P), K_star_guess(P), opts_deb);


% Adjustment costs labor

% Adjustment cost functions
if pay_for_delta == 0
    phi_l = @(I_l) psi_l/2 .* (I_l) .^2;
%     phi_l = @(I_l) psi_k/2 .* (I_l - sep) .^2;
else
    phi_l = @(I_l) psi_l/2 .* (I_l) .^2;
end
% Derivatives

if pay_for_delta == 0
    dphi_l = @(I_l) psi_l .* (I_l) ;
%     dphi_l = @(I_l) psi_l .* (I_l - sep) ;
else
    dphi_l = @(I_l) psi_l .* (I_l) ;
end
% Inverse derivatives
if pay_for_delta == 0
    inv_dphi_l = @(y) 1/psi_l .* y;
%     inv_dphi_l = @(y) 1/psi_l .* y + sep ;
else
    inv_dphi_l = @(y) 1/psi_l .* y;
end



%% Grid parameters

% Size
N_k = params.N_k; % capital
N_l = params.N_l; % labor

if stochastic == 1
    N_p = N_p; % Poisson states for the price
else
    N_p = 1;
end


if strcmp(ShockType, 'Lognormal')
    
    %% POISSON PROCESS
    settings.nz = N_p;              % number of grid points
    % settings.nz = 201;              % number of grid points
    settings.rho = rho_P;           % autoregressive coefficient
    settings.sigma = sigma_p_eps;   % standard deviation of the price error
    settings.mu = mu_p;             % offset of the process
    settings.criterion = 1e-6;      % convergence for stationary distribution
    settings.nostat = true;         % flag not to compute stationary distribution
    if isfield(params.settings, 'm')
        settings.m = params.settings.m; % Standard deviations of interest (Tauchen)
    else
        settings.m = 3;
    end


    % % Rouwenhorst Method for Discretizing Highly Persistent AR Processes
    % [TM_large, states_large, ~] = discretizeAR(settings);
    % states_large = states_large - mean(states_large);
    % [TM_large, states_large] = discretizeAR_tauchen(settings);
    % states_large = states_large - mean(states_large);

    [TM, states] = discretizeAR_tauchen(settings);
elseif strcmp(ShockType, 'Markovian')
    TM = params.TM;
    states = log(params.states);
elseif  strcmp(ShockType, 'Diffusion')    
    settings = params.settings;
    settings.N_p = N_p;
    if isfield(params, 'pmin')
        settings.pmin = params.pmin;
    else
        settings.pmin = 0.01;
    end
    if isfield(params, 'pmax')
        settings.pmax = params.pmax;
    else
        settings.pmax = 10;
    end
    if isfield(params, 'p_grid')
        settings.p_grid = params.p_grid;
    end
    
    if settings.sigma(settings.pmax) == 0
        N_p = 1;
        p_grid = P;
        pmin = min(p_grid);
        pmax = max(p_grid);
        pi0 = 1;
        Lambda = 0;
        stochastic = 0;
        outDiffusion.dp_mean = 1;
        dp_mean = 1;
    else
        outDiffusion = DiscretizeDiffusion(settings);
        p_grid = outDiffusion.p_grid;
        pmin = min(p_grid);
        pmax = max(p_grid);
        % standardize so that trapz(p_grid,pi0) = 1
        pi0 = outDiffusion.g_0./outDiffusion.dp_mean';
        Lambda = outDiffusion.Lambda;
    end
    
else
    error('Specify ShockType, accepts Markovian, Lognormal or Diffusion')    
end   

if ~strcmp(ShockType, 'Diffusion')
    % Compute Stat Dist of Exogenous Shock
    criterion_statdist = 1e-6;
    distance_statdist = criterion_statdist + 1;
    pi0 = 1/N_p * ones(1,N_p); 

    while distance_statdist > criterion_statdist
        temp = pi0;
        pi0 = pi0 * TM;
        distance_statdist = max(abs(temp-pi0));
    end

    pi0 = pi0';

    out.params.TM = TM;

    % states = log(params.p_grid);

    % Matrix containing Poisson intensities from (row_state) to (col_state)
    Lambda = volatility .* (TM - diag(diag(TM)) - diag(sum(TM - diag(diag(TM)), 2))) .* time_scale;

    % Input Validation
    if size(Lambda,1) ~= size(Lambda,2)
        error('Lambda is not square');
    end
    if size(Lambda,1) < N_p
        error('Lambda is too small');
    end
    if min(sum(Lambda,2)) > 0
        error('Flows not summing to 0!');
    end

   %% Boundaries Price


    if stochastic
        p_grid = exp(states); % for three states
        p_grid = p_grid ./ (p_grid * pi0); % Eliminate numerical errors
        pmin = min(p_grid);
        pmax = max(p_grid);
    else
        p_grid = P;
        pmin = min(p_grid);
        pmax = max(p_grid);
    end

end


 %% Boundaries Capital
kmin = 0;

if isfield(params, 'kmax')
    kmax = params.kmax;
else
    [kmax, ~, EXITFLAG] = K_star(pmax); % set to the unconstrained solution with minimal labor share and max price
    if EXITFLAG <= 0
        kmax = fsolve(@(k) solve_max_k(k,pmax), 1, opts);
    end
    if kmax < 0 || ~isreal(kmax)
        kmax = K_star_guess(pmax);
    end
end

out.params.kmax = kmax;

 %% Boundaries Labor
 
 lmin = 0;
 
if isfield(params, 'lmax')
    lmax = params.lmax;
else
    lmax = (1 - Gamma)/(Gamma) * k_bar(pmax);
end

out.params.lmax = lmax;

%% Convergence parameters
Delta = 1000;
tol = 1e-6;
maxit = 50; % maximum value function iterations
%% grids

%option to have unevenly spaced grid Gamma > 1 gives more points for lower
%K
if isfield(params, 'GridKGamma')
    GridKGamma = params.GridKGamma;
else
    GridKGamma = 1;
    out.params.GridKGamma = GridKGamma;
end

if isfield(params, 'GridLGamma')
    GridLGamma = params.GridLGamma;
else
    GridLGamma = 1;
    out.params.GridLGamma = GridLGamma;
end


% grids
k_grid = kmin + (kmax - kmin) * (linspace(0, 1, N_k).^GridKGamma)';
l_grid = lmin + (lmax - lmin) * (linspace(0, 1, N_l).^GridLGamma)';

% in single vector form with p outer state, then a, then k

kk = repmat(k_grid, N_l, 1);
ll = kron(l_grid, ones(N_k,1));

kk = repmat(kk, N_p, 1);
ll = repmat(ll, N_p, 1);
pp = kron(p_grid' , ones(N_k*N_l , 1));


%DELTAS for k
dk_vec_up = diff(k_grid);
dk_vec_up = [dk_vec_up; dk_vec_up(end)];
dk_vec_dwn = diff(k_grid) ;
dk_vec_dwn = [dk_vec_dwn(1); dk_vec_dwn];

dk_vec_up = repmat(dk_vec_up, N_l, 1);
dk_vec_up = repmat(dk_vec_up, N_p, 1);

dk_vec_dwn = repmat(dk_vec_dwn, N_l, 1);
dk_vec_dwn = repmat(dk_vec_dwn, N_p, 1);

dk_mean = (dk_vec_up + dk_vec_dwn)/2;
dk_mean(1) = dk_mean(1)/2;
dk_mean(end) = dk_mean(end)/2;

%DELTAS for l
dl_vec_up = diff(l_grid);
dl_vec_up = [dl_vec_up; dl_vec_up(end)];
dl_vec_dwn = diff(l_grid) ;
dl_vec_dwn = [dl_vec_dwn(1); dl_vec_dwn];

dl_vec_up = kron(dl_vec_up, ones(N_k,1));
dl_vec_up = repmat(dl_vec_up, N_p, 1);

dl_vec_dwn = kron(dl_vec_dwn, ones(N_k,1));
dl_vec_dwn = repmat(dl_vec_dwn, N_p, 1);

dl_mean = (dl_vec_up + dl_vec_dwn)/2;
dl_mean(1) = dl_mean(1)/2;
dl_mean(end) = dl_mean(end)/2;


%DELTAS for p

if stochastic
    dp_vec_up = diff(p_grid');
    dp_vec_up = [dp_vec_up; dp_vec_up(end)];
    dp_vec_dwn = diff(p_grid') ;
    dp_vec_dwn = [dp_vec_dwn(1); dp_vec_dwn] ;
    
    dp_vec_up = kron(dp_vec_up , ones(N_k*N_l , 1));
    dp_vec_dwn = kron(dp_vec_dwn , ones(N_k*N_l , 1));

    dp_mean = (dp_vec_up + dp_vec_dwn)/2;
    dp_mean(1) = dp_mean(1)/2;
    dp_mean(end) = dp_mean(end)/2;
end


%% DEFINE TRANSITIONS

% CAPITAL GRID
% create indices (and implicitly set boundary conditions)
ind = (1:1: N_k * N_l * N_p )'; % starting state indices

vector_ind = ind2sub_vec([N_k N_l N_p], ind); %matrix with indices for all states
k_up = (N_k  - vector_ind(:,1)) > 0 ; %positions where it is possible to go up
k_dwn = ( vector_ind(:,1) -1) > 0 ; %positions where it is possible to go down


%transition indices if drift is positive (up) or negative (dwn)
% set indices so that if we cannot move up, the "up" index is the same as
% the original index, similar for down
ind_up_k = sub2ind_vec([N_k  N_l  N_p], vector_ind + [k_up zeros(length(ind),2)]);
ind_dwn_k = sub2ind_vec([N_k N_l  N_p], vector_ind - [k_dwn zeros(length(ind),2)]);


% LABOR GRID

l_up = (N_l  - vector_ind(:,2)) > 0 ; %positions where it is possible to go up
l_dwn = ( vector_ind(:,2) -1) > 0 ; %positions where it is possible to go down


%transition indices if drift is positive (up) or negative (dwn)
% set indices so that if we cannot move up, the "up" index is the same as
% the original index, similar for down
ind_up_l = sub2ind_vec([N_k  N_l  N_p], vector_ind + [zeros(length(ind),1) l_up zeros(length(ind),1)]);
ind_dwn_l = sub2ind_vec([N_k N_l  N_p], vector_ind - [zeros(length(ind),1) l_dwn zeros(length(ind),1)]);


if ~strcmp(ShockType, 'Diffusion')

     % Definition of ind_up selecting only the instances where the state can
     % ACTUALLY jump up

     ind_up_p = []; % for the colums the state jumps to
     ind_dwn_p = []; % for the colums the state jumps to
     ind_p = []; % for the rows the state jumps from
     F_p = [];
     B_p = [];
     DIAG_p = [];


     %% JUMPS UP
     for shift_state = 1 : N_p - 1 %over max number of possible jumps up 
         shift_possible = vector_ind(: , 3) + shift_state <= N_p; % indicator u-bound
         ind_col_up = sub2ind_vec([N_k  N_l  N_p],...
             vector_ind + [zeros(length(ind),2) shift_possible * shift_state ]);
         % select only states where it actually goes up
         ind_col_up = ind_col_up(logical(shift_possible));
         ind_row_up = ind(logical(shift_possible));

         ind_up_p = [ind_up_p;  ind_col_up];
         ind_p = [ind_p;  ind_row_up];

         for state = 1: (N_p - shift_state) %loop for states where shift_state is possible
             % get off-diag. elms that are above current state
             F_p = [F_p; repmat( Lambda(state, state + shift_state),N_k * N_l, 1)];
         end
     end
     %% JUMPS DOWN
    for shift_state = 1 : N_p - 1 %over max number of possible jumps down
         shift_possible = vector_ind(: , 3) - shift_state >= 1; % indicator l-bound
         ind_col_dwn = sub2ind_vec([N_k  N_l  N_p],...
             vector_ind - [zeros(length(ind),2) shift_possible * shift_state ]);
         % select only states where it actually goes up
         ind_col_dwn = ind_col_dwn(logical(shift_possible));
         ind_row_dwn = ind(logical(shift_possible));

         ind_dwn_p = [ind_dwn_p;  ind_col_dwn];
         ind_p = [ind_p;  ind_row_dwn];

         for state = (shift_state + 1): N_p %loop for states where jumps possible
             % get off-diag. elms that are below current state
             B_p = [B_p; repmat( Lambda(state, state - shift_state), N_k * N_l, 1)];
         end
     end 
    %% DIAGONAL
    if stochastic
         for state = 1: N_p
             DIAG_p = [DIAG_p; repmat( Lambda(state, state), N_k * N_l, 1)];
         end
    else
        DIAG_p = zeros( N_k * N_l, 1);
    end


    %  %% GRAPH for states only
    %  sparse_rows = [ind; ind_p];
    %  sparse_cols = [ind; ind_up_p; ind_dwn_p];
    %  sparse_vals = [DIAG_p; F_p; B_p]; % will set  inside loop
    %  A_new = sparse(sparse_rows, sparse_cols, sparse_vals);
    % 
    %  figure()	
    %  spy(A_new)
    %  title('new method')
else
    
    if stochastic
        F_p = kron(outDiffusion.F_p' , ones(N_k * N_l , 1));
        B_p = kron(outDiffusion.B_p' , ones(N_k * N_l , 1));
    
        DIAG_p = - F_p - B_p;
    else
         DIAG_p = zeros(N_k*N_l, 1);
         F_p = DIAG_p;
         B_p = DIAG_p;
    end
    
    
    p_up = (N_p  - vector_ind(: , 3)) > 0 ; %positions where it is possible to go up
    p_dwn = ( vector_ind(: , 3) -1) > 0 ; %positions where it is possible to go down


    %transition indices if drift is positive (up) or negative (dwn)
    % set indices so that if we cannot move up, the "up" index is the same as
    % the original index, similar for down
    ind_up_p = sub2ind_vec([N_k  N_l  N_p], vector_ind + [zeros(length(ind),2) p_up]);
    ind_dwn_p = sub2ind_vec([N_k  N_l  N_p], vector_ind - [zeros(length(ind),2) p_dwn]);
    ind_p = repmat(ind, 2, 1);
    


end
%% ALLOCATION
% set up the sparse rows and columns for the matrices
sparse_rows = [repmat(ind, 5 ,1); ind_p];
sparse_cols = [ind; ind_up_k;  ind_dwn_k;...
    ind_up_l;  ind_dwn_l; ...
     ind_up_p;  ind_dwn_p];
 

% CHECK the approach makes sense for forward approximation and jump up of
% Poisson

% sparse_rows = repmat(ind,4,1);
% sparse_cols = [ind; ind_up_k;  ind_up_a; ind_up_p];
% DIAG = - ones(length(ind),1);
% F_k = -DIAG;
% F_a = -DIAG;
% F_p = -DIAG;
% sparse_vals = [DIAG; F_k;F_a; F_p];
% 
% A = sparse(sparse_rows, sparse_cols, sparse_vals);
% figure()
% spy(A)
%% SOLVE WITH IMPLICIT METHOD

Pi_vec = Pi(kk, ll,  pp);
% Initialize the value function: assume no positive drift in initial states,
% value is just current profit minus upkeep of stocks


v_0 = (Pi_vec - phi_k(delta.*kk) - phi_l(sep.*ll))./rho;

V_next = v_0;


supnorm = zeros(maxit, 1);

for J = 1:maxit
    
    %% FORWARD DIFFERENCES
    V = V_next;
    % Forward differences: get a zero when the step forward is not possible
    dVF_k = (V(ind_up_k) - V(ind))./dk_vec_up;
    dVB_k = - (V(ind_dwn_k) - V(ind))./dk_vec_dwn;
    
    dVF_l = (V(ind_up_l) - V(ind))./dl_vec_up;
    dVB_l = - (V(ind_dwn_l) - V(ind))./dl_vec_dwn;

    % Check it is correct
    % [kk V (V(ind_up_k) - V(ind))]
    % [aa V (V(ind_up_a) - V(ind))]
    % [pp V (V(ind_up_p) - V(ind))]

    %% STATE CONSTRAINTS : HERE

    % Impose state constraints to correctly compute the optimal controls.
    % Set the up_drift at 0 at upper boundary, so that 
    % k_dot = (I - delta)*k -> I = delta
    % ind_up_k == ind shows where upward movement is not possible
    
    if pay_for_delta ==1
        
       
        % Forward derivative
        dVF_k(ind_up_k == ind) = dphi_k(delta * kmax) + R;

        % Backward derivative
        dVB_k(ind_dwn_k == ind) = dphi_k(delta * kmin) + R;

        % Optimal controls
        i_F_k = inv_dphi_k(dVF_k - R);
        i_B_k = inv_dphi_k(dVB_k - R);
        
        % Forward derivative
        dVF_l(ind_up_l == ind) = dphi_l(delta * lmax);

        % Backward derivative
        dVB_l(ind_dwn_l == ind) = dphi_l(delta * lmin);

        % Optimal controls
        i_F_l = inv_dphi_l(dVF_l);
        i_B_l = inv_dphi_l(dVB_l);

        % Impose positive investment constraint
        if strcmp(ReversibleInvestment, 'No')
            i_F_k = (i_F_k >= 0) .* i_F_k;
            i_B_k = (i_B_k >= 0) .* i_B_k;
        end
    else
        % Forward derivative
        dVF_k(ind_up_k == ind) = dphi_k(delta * kmax) - psi_k * delta * kmax + R;

        % Backward derivative
        dVB_k(ind_dwn_k == ind) = dphi_k(delta * kmin) - psi_k * delta * kmin + R;

        % Optimal controls
        i_F_k = inv_dphi_k(dVF_k + psi_k * delta * kk -  R);
        i_B_k = inv_dphi_k(dVB_k + psi_k * delta * kk - R);
%         i_0_k = inv_dphi_k(.5*(dVF_k + dVB_k) + psi_k * delta * kk - R);
        i_0_k = delta * kk;
        % Impose positive investment constraint
        if strcmp(ReversibleInvestment, 'No')
            i_F_k = (i_F_k >= 0) .* i_F_k;
            i_B_k = (i_B_k >= 0) .* i_B_k;
        end
        
        % Forward derivative
        dVF_l(ind_up_l == ind) = dphi_l(sep * lmax) - psi_l * sep * lmax;

        % Backward derivative
        dVB_l(ind_dwn_l == ind) = dphi_l(sep * lmin) - psi_l * sep * lmin;
        
        i_F_l = inv_dphi_l(dVF_l + psi_l * sep * ll);
        i_B_l = inv_dphi_l(dVB_l + psi_l * sep * ll);
%         i_0_k = inv_dphi_k(.5*(dVF_k + dVB_k) + psi_k * delta * kk - R);
        i_0_l = sep * ll;
        
    end

    % Drifts
    s_F_k = i_F_k - delta .* kk;
    s_B_k = i_B_k - delta .* kk;
    
    s_F_l = i_F_l - sep .* ll;
    s_B_l = i_B_l - sep .* ll;

    % instantaneous costs of investment
    if pay_for_delta == 1
        Phi_F_k =  phi_k(i_F_k) + R * i_F_k;
        Phi_B_k =  phi_k(i_B_k) + R * i_B_k;
        Phi_0_k =  phi_k(delta * kk) + R * (delta * kk);
        
        Phi_F_l =  phi_l(i_F_l);
        Phi_B_l =  phi_l(i_B_l);
        Phi_0_l =  phi_l(sep * ll);
        
    else
        Phi_F_k =  phi_k(i_F_k - delta .* kk) + R * i_F_k;
        Phi_B_k =  phi_k(i_B_k - delta .* kk) + R * i_B_k;
        Phi_0_k =  phi_k(i_0_k - delta .* kk) + R * i_0_k;
        
        Phi_F_l =  phi_l(i_F_l - sep * ll);
        Phi_B_l =  phi_l(i_B_l - sep * ll);
        Phi_0_l =  phi_l(i_0_l - sep * ll);
    end

    % indicators
    I_F_k = s_F_k > 0;
    I_B_k = s_B_k < 0;
    I_0_k = (1 - I_F_k - I_B_k);
    
    I_F_l = s_F_l > 0;
    I_B_l = s_B_l < 0;
    I_0_l = (1 - I_F_l - I_B_l);

%     % Correction when there is ambiguity on where to move, use maximum drift of
%     % the value function
%     ind_amb_k = ((I_F_k == 1) & (I_B_k == 1)); %indices where there is ambiguitiy
%     decision_k = abs(dVF_k) > abs(dVB_k);
%     I_F_k(ind_amb_k) = decision_k(ind_amb_k);
%     I_B_k(ind_amb_k) = abs(1 - decision_k(ind_amb_k));
% 
%     I_0_k = (1 - I_F_k - I_B_k);

    % instantaneous utility
    if strcmp(InvGoodCost, 'Final')
        u = Pi_vec + p_y * ( - I_F_k .* Phi_F_k - I_B_k .* Phi_B_k - I_0_k .* Phi_0_k - I_F_l .* Phi_F_l - I_B_l .* Phi_B_l - I_0_l .* Phi_0_l);
    else
        u = Pi_vec + ...
            ( - I_F_k .* Phi_F_k - I_B_k .* Phi_B_k - I_0_k .* Phi_0_k  ...
            - I_F_l .* Phi_F_l - I_B_l .* Phi_B_l - I_0_l .* Phi_0_l );
    end
        


    % Construct the sparse matrix
    F_k = I_F_k .* s_F_k ./ dk_vec_up; % positive drifts that are actually positive
    B_k = I_B_k .* s_B_k ./ dk_vec_dwn; 
    
    F_l = I_F_l .* s_F_l ./ dl_vec_up; % positive drifts that are actually positive
    B_l = I_B_l .* s_B_l ./ dl_vec_dwn; 

    DIAG = DIAG_p - F_k + B_k - F_l + B_l ;
    % Populate A

    A = sparse(sparse_rows, sparse_cols, [DIAG ;...
        F_k; -B_k; F_l; -B_l; F_p; B_p]);
    
    %% UPDATE VALUE
    B = (1/Delta + rho)*speye(N_p * N_k * N_l) - A;
    b = u + V/Delta;
    
%     Bgpu = gpuArray(B);
%     bgpu = gpuArray(b);
%     V_next = gather(mldivide(Bgpu,bgpu));
    V_next = B\b;

    
%     V_next = cgs(B,b);
    supnorm(J) = max(abs(V - V_next) ./ V_next * 100);
    if print_conv
        disp(['HJB Iteration ' num2str(J) '. Supnorm = ' num2str(supnorm(J)) '.'])
    end
    if supnorm(J) < tol   
        break
    end

end




%% STATIONARY DISTRIBUTION STARTING FROM g_0;
AT = A';
if isfield(params, 'DT')
    DT = params.DT;
else
    DT = 1000;
    out.params.DT = DT;
end

if isfield(params, 'g_0')
    g_0 = params.g_0;
else
    g_0 = zeros(N_k*N_l*N_p, 1);
    if strcmp(ShockType, 'Diffusion') 
%         g_0(1:N_k*N_l:N_k*N_l*N_p) = pi0 ./ sum(pi0);
        g_0(1:N_k*N_l:N_k*N_l*N_p) = pi0.*outDiffusion.dp_mean';% start everyone at 0 capital on statdist of p
    else
        g_0(1:N_k*N_l:N_k*N_l*N_p) = pi0;
    end
    
    %standardize for non-uniform grids
%     g0_mat = reshape(g_0, [N_k N_p]);
    
%     if strcmp(ShockType, 'Diffusion')
%         g_0(1:N_k:N_k*N_p) = g_0(1:N_k:N_k*N_p)./(trapz(k_grid, trapz(p_grid, g0_mat, 2)));
%     else
%         g_0(1:N_k:N_k*N_p) = g_0(1:N_k:N_k*N_p)./(trapz(k_grid, sum(g0_mat, 2)));
%     end
    
    out.params.g_0 = g_0;
end

if isfield(params, 'tolG')
    tol_kfe = params.tolG;
else
    tol_kfe = 1e-6;
    out.params.tolG = tol_kfe;
end
    

II = speye(size(AT));
TransMat = plus(II, - DT*AT);

gg = ones(size(g_0));
converged = false;
it_kfe = 0;
while ~converged
    
    it_kfe = it_kfe + 1;
    gg = TransMat \  g_0;

    supnorm_kfe = max(abs(gg-g_0));
    g_0 = gg;
    
    if print_conv
        disp(['KFE Iteration: ' num2str(it_kfe) '. Supnorm = ' num2str(supnorm_kfe) '.'])
    end
    
    converged = supnorm_kfe < tol_kfe;
    
end

out.statDist_pmf = gg;

%rescale to obtain pdf
if strcmp(ShockType, 'Diffusion')
    
    gg = gg./(dp_mean.*dk_mean.*dl_mean);
%     gg = gg./(dp_mean.*dk_mean.*dl_mean); 
else
    gg = gg./(dk_mean.*dl_mean);
end
 
 
%  %Check that the sum is 1
% if strcmp(ShockType, 'Diffusion')
%     if abs(sum(gg.*dp_mean.*dk_mean) - 1) > 1e-6
%         warning(['StatDist is summing to ' num2str(sum(gg.*dp_mean.*dk_mean))]);
%     end
% else
%     if abs(sum(gg.*dk_mean) - 1) > 1e-6
%         warning(['StatDist is summing to ' num2str(sum(gg.*dk_mean))]);
%     end
% end 


%% DEFINE INTEGRALS
if strcmp(ShockType, 'Diffusion')
    integrate_kp = @(fun) sum(fun.*dk_mean.*dp_mean.*dl_mean);
else
    integrate_kp = @(fun) sum(fun.*dk_mean.*dl_mean);
end

%% POLICIES and STATISTICS

%Investment policy for each firm
i_k_pol = I_F_k .* i_F_k + I_B_k .* i_B_k + I_0_k * delta .* kk;
i_l_pol = I_F_l .* i_F_l + I_B_l .* i_B_l + I_0_l * delta .* ll;
% i_k_pol = I_F_k .* i_F_k + I_B_k .* i_B_k + I_0_k * i_0_k;

% Labor demand
% for each firm
indL = ((k_bar(pp) - kk)>=0) .*  ( ((k_bar(pp) - kk)  *(1 - Gamma)/ Gamma <=1) .* ...
    (k_bar(pp) - kk)  *(1 - Gamma)/ Gamma  +...
    ((k_bar(pp) - kk)  *(1 - Gamma)/ Gamma >1) .* 1) .* ll;
% overall
L = integrate_kp(gg .* indL);
% stock of labor
L_contracts = integrate_kp(gg .* ll);
out.L_contracts = L_contracts;

% If workers produce robots, add robot purchases to labor demand
if strcmp(RobotsProducers, 'Workers')
    % Assume production function is R = A*L Therefore add I_k/A to L^d
    L = L + integrate_kp(gg .* i_k_pol/A_robots);
end

%compute the mass of firms fully automating
gg_auto = gg;
gg_auto(kk<k_bar(pp)) = 0;

out.Mass_auto = integrate_kp(gg_auto);
out.gg_auto = gg_auto;


% compute measure of automation in eqm
out.auto_extent = integrate_kp(gg .* min(kk ./ k_bar(pp), 1));

out.kmaxgg = max(kk(gg>tol_kfe));

% Optional output
out.statDist = gg;
out.A = A;
out.V = V_next;


% Decision to automate or not: derivative of value function at R=0 in
% maximum productivity state (using only Forward Difference because at the
% boundary) must be positive! 
% out.INV_TRADEOFF_MAX_P = dVF_k(N_k*(N_p - 1) + 1) - R;

%% INV CUTOFFS COMMENTED

% Compute Investment Cutoffs by state
% Note: correct approximation of V' must use UPWIND SCHEME!!!
dV = dVF_k .* I_F_k + I_B_k .* dVB_k + .5 .* I_0_k .* (dVF_k + dVB_k);
%increment to use for Inv_Cutoffs
%dk_dV = dk_vec_up .* I_F_k + dk_vec_dwn .* dVB_k + .5 .* I_0_k .* (dk_vec_up + dk_vec_dwn);

%{
if pay_for_delta == 1 %find the R_star where the inv. covers depreciation
    dV_reshaped = reshape(inv_dphi_k(dV - R) - delta * kk, [N_k, N_p]);
else %find where the drift equals 0 without paying for depreciation
    dV_reshaped = reshape(inv_dphi_k(dV - R), [N_k, N_p]);
end

out.INV_CUTOFFS = zeros(N_p,1);
% Linearly interpolate the investment condition to find the robots level
% such that the firm keeps a constant size. (perfectly offset depreciation)
for idx = 1:N_p
    idx_cutoff = find(dV_reshaped(:,idx) < 0, 1);
    if idx_cutoff > 1
        idx_cutoff_lin = sub2ind([N_k, N_p], idx_cutoff, idx);
        out.INV_CUTOFFS(idx) = kk(idx_cutoff_lin - 1) - ...
           (kk(idx_cutoff_lin) - kk(idx_cutoff_lin - 1)) * dV_reshaped(idx_cutoff-1,idx) / ...
                          (dV_reshaped(idx_cutoff,idx) - dV_reshaped(idx_cutoff-1,idx));
    end
end
out.dV = dV_reshaped;

% Alternative computation of Investment Cutoffs

netInv = reshape(i_k_pol - delta*kk, N_k, N_p);
try
    out.INV_CUTOFFS_ALTERNATIVE = zeros(N_p,1);
    for idx = 1:N_p
        out.INV_CUTOFFS_ALTERNATIVE(idx) = interp1(netInv(:,idx),k_grid,0);
    end

    out.INV_CUTOFFS_MEAN = .5*(out.INV_CUTOFFS + out.INV_CUTOFFS_ALTERNATIVE);
catch
    out.INV_CUTOFFS_MEAN = out.INV_CUTOFFS;
    out.INV_CUTOFFS_ALTERNATIVE = out.INV_CUTOFFS;
end
%}
%% OUT GE OBJECTS: COMMENTED

%{
%Revenues
if strcmp(Utilization,'Yes')
    Rev = @(k,P) (k <= k_bar(P)) .* (...
         (1 - Gamma)^t * A_prod * P  .* (k_bar(P)).^(t)  ) + ...
        (k > k_bar(P)) .* (k < k_hat(P)) .* (...
        A_prod * (1 - Gamma)^t * P .* k .^t ) + ...
        (k >= k_hat(P)) .* ( A_prod * (1 - Gamma)^t * P .* (k_hat(P)) .^t ) ;

    if E == 0
        
        Rev = @(k,P) (k <= k_bar(P)) .* (...
                    (1 - Gamma)^t * A_prod * P  .* (k_bar(P)).^(t)  ) + ...
                    (k > k_bar(P)) .* (...
                    A_prod * (1 - Gamma)^t * P .* k .^t );
        
    end
    
    IndRev = Rev(kk,pp);
end

% CHOICE OF INVESTMENT COSTS: if InvGoodCost = FINAL value in final good price, else in intermediate

% Individual adjustment costs at the optimal policy
if ~strcmp(InvGoodCost, 'Final')
    if pay_for_delta == 1
        IndAdjCosts = phi_k(i_k_pol);
    else
        IndAdjCosts = phi_k(i_k_pol - delta * kk);
    end
else
    if pay_for_delta == 1
        IndAdjCosts = p_y * phi_k(i_k_pol);
    else
        IndAdjCosts = p_y * phi_k(i_k_pol - delta * kk);
    end
end


% REVENUES NET OF INVESTMENT COSTS IF PAID IN INTERMEDIATES (DEFAULT)
if ~strcmp(InvGoodCost, 'Final') && ~strcmp(InvGoodCost, 'External')
    netRev = IndRev - IndAdjCosts;
else %costs in final good, so the convention on Q is that it is gross of costs
    netRev = IndRev;
end

if ~strcmp(EnergyCosts,'Final') && ~strcmp(InvGoodCost, 'External') %Energy costs are paid in intermediates

    %compute intermediate good quantities and profits in the various scenarios
    
    if strcmp(RobotsProducers,'None') % production of robots is paid in intermediates: Net output excludes robot purchases
        indQ = netRev - R * i_k_pol - E * kk; 
        indPi = indQ - W * indL;
    end

    if strcmp(RobotsProducers,'Workers') % production of robots is paid to robot producing sectors by workers
        indQ = netRev - E * kk;
        indPi = indQ - W * indL - R * i_k_pol;
    end
    
    if strcmp(RobotsProducers,'Capitalists') % production of robots is paid out to capitalists as income
        indQ = netRev - E * kk;
        indPi = indQ - W * indL;
    end
    
else % IF THERE IS AN ENERGY SECTOR, DO THE SAME AS BEFORE, NOW NET OUTPUT DOES NOT EXCLUDE ENERGY COSTS
    %compute intermediate good quantities and profits in the various scenarios
    if strcmp(RobotsProducers,'None') % production of robots is paid in intermediates: Net output excludes robot purchases
        indQ = netRev - R * i_k_pol;
        indPi = indQ - W * indL - E * kk;
    end

    if strcmp(RobotsProducers,'Workers') % production of robots is paid to robot producing sectors by workers
        indQ = netRev;
        indPi = indQ - W * indL - R * i_k_pol - E * kk;
    end
    
    if strcmp(RobotsProducers,'Capitalists')  % production of robots is paid out to capitalists as income
        indQ = netRev;
        indPi = indQ - W * indL - E * kk;
    end    
    
end

if strcmp(InvGoodCost, 'Final') || strcmp(InvGoodCost, 'External')
%costs in final good, so the profits bear the burden of adjustment costs
%(multiplied by p_y in the definition of IndAdjCosts above
    indPi = indPi - IndAdjCosts; 
end
    
% Compute Useful GE objects
out.GE.Pi = integrate_kp(gg .* indPi); % Aggregate Profits (less aggregate adjustment costs)
out.GE.q = integrate_kp(gg .* indQ); % Aggregate Output
out.GE.IndRev = IndRev;
out.GE.IndAdjCosts = IndAdjCosts;
out.GE.AggAdjCosts = integrate_kp(gg .*IndAdjCosts);
out.GE.netRev = netRev;
out.GE.AggnetRev = integrate_kp(gg .*netRev);
out.GE.AggRev = integrate_kp(gg .* IndRev);
out.GE.AggECosts = integrate_kp(gg .* E .* kk);
out.GE.indQ = indQ;
out.GE.indPi = indPi;
%}
%% OUT OTHER PARAMETERS
out.i_k_pol = i_k_pol;
out.i_l_pol = i_l_pol;
out.kk = kk;
out.pp = pp;
out.ll = ll;
out.k_grid = k_grid;
out.p_grid = p_grid;
out.l_grid = l_grid;
out.W = W;
out.L_pol = indL;
out.u = u;
out.params.pi0 = pi0;
out.params.Lambda = Lambda;
out.L = L;
if ~strcmp(RobotsProducers,'None')
    out.L_robots_sector = integrate_kp(gg .*i_k_pol/A_robots);
else
    out.L_robots_sector = 0;
end
out.L_intermediate_goods_sector = out.L - out.L_robots_sector;
out.GE.AggInv = integrate_kp(gg .*i_k_pol);
out.params.Omega = ( 1 - Gamma)/Gamma * W - E;

out.RL = integrate_kp(gg .* kk)./out.L * 1e3; % Robots per thousand employees
%out.LS_AggRev = W * L / out.GE.AggRev;
%out.LS_q = W * L / out.GE.q;
%out.YL = out.GE.AggRev / L;


%% OUT USEFUL FUNCTION HANDLES
out.handles.Pi = Pi;
out.handles.R_bar = k_bar;
if strcmp(Utilization,'Yes')
    out.handles.R_hat = k_hat;
    out.handles.Pi_prime = Pi_prime;
end
out.handles.R_hat_hat = K_star;
out.handles.Phi = phi_k;
out.handles.dPhi = dphi_k;
out.handles.dPhi_inv = inv_dphi_k;
out.handles.p_bar = p_bar;
out.handles.p_hat = p_hat;

out.handles.integrate_kp = integrate_kp;

if strcmp(ShockType, 'Diffusion')
    out.pi0_int = pi0.*outDiffusion.dp_mean';
else
    out.pi0_int = pi0;
end


out.params.sigma_P = sqrt(out.p_grid.^2*out.pi0_int-(out.p_grid*out.pi0_int)^2);

% 
out.dV_Fk = dVF_k;
out.dV_Bk = dVB_k;
out.dV_Fl = dVF_l;
out.dV_Bl = dVB_l;

% 
% dk_vec_up = diff(k_grid);
% dk_vec_up = [dk_vec_up; dk_vec_up(end)];
% dk_vec_dwn = diff(k_grid) ;
% dk_vec_dwn = [dk_vec_dwn(1); dk_vec_dwn];
% 
% 
% dk_mean = (dk_vec_up + dk_vec_dwn)/2;
% dk_mean(1) = dk_mean(1)/2;
% dk_mean(end) = dk_mean(end)/2;
% 
% %DELTAS for p
% dp_vec_up = diff(p_grid');
% dp_vec_up = [dp_vec_up; dp_vec_up(end)];
% dp_vec_dwn = diff(p_grid') ;
% dp_vec_dwn = [dp_vec_dwn(1); dp_vec_dwn] ;
% 
% 
% dp_mean = (dp_vec_up + dp_vec_dwn)/2;
% dp_mean(1) = dp_mean(1)/2;
% dp_mean(end) = dp_mean(end)/2;
% 
% 
% out.dp_mean = dp_mean;
% out.dk_mean = dk_mean;

if stochastic == 0
    out.params.N_p = 1;
    out.params.stochastic = stochastic;
end

end





    