% Fornino Manera (2019)
% 
% DATE: January 24, 2021
%
% Project: Automation and the Future of Work: Assessing the Role of Labor
%          Flexibility

function [L, out] = LaborDemand_trapz(W, R, params)
% This function finds the solution to the firm problem as well as the
% stationary distribution for given parameters.
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

% GBM reset value of capital option
if strcmp(ShockType, 'GBM')
    if isfield(params,'ReEntryK')
        ReEntryK = params.ReEntryK;
    else
        ReEntryK = 'Kstar';
    end
end
if isfield(params, 'AdjAtEntry')
    AdjAtEntry = params.AdjAtEntry;
else
    AdjAtEntry = true;
end   


%% Choice for Putty-Clay v. Utilization
if isfield(params,'Utilization')
    Utilization = params.Utilization;
else
    Utilization = 'Yes';
    out.params.Utilization = Utilization;
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
% discount rate 
rho = params.rho;
% Adjustment cost parameters
psi_k = params.psi_k; % capital

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
%% Cost options
% Pay for delta?

if isfield(params, 'pay_for_delta')
    pay_for_delta = params.pay_for_delta;
else
    pay_for_delta = 1;
end

% Kinked?
if isfield(params, 'kinked')
    kinked = params.kinked;
else
    kinked = false;
end
if isfield(params, 'chi_k')
    chi_k = params.chi_k;
else
    if kinked
        error('please specify parameter chi_k when kinked adj costs are assumed.')
    end
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

    Pi = @ (k , P) (k <= k_bar(P)) .* (...
        (1 - t) * (1 - Gamma)^t * A_prod * P  .* (k_bar(P)).^(t) ...
        + (W * (1-Gamma) / Gamma - E) * k ) + ...
        (k > k_bar(P)) .* (...
        A_prod * (1 - Gamma)^t * P .* k .^t - E * k);
else
    Utilization_low_region = ((1-Gamma)/Gamma * W - E)>=0;
% AMENDED PROFIT AND KBAR if there is a utilization margin
    k_bar = @(P) 1/(1-Gamma)*( W./(A_prod*P*Gamma*t)).^(1/(t-1)) ;
    k_hat = @(P) ( (A_prod*P*(1-Gamma)^t*t)./E).^(1/(1-t)) ;

    p_hat = @(R) E / (A_prod * t * (1 - Gamma)^t) .* R .^ (1 - t);
    p_bar = @(R) W / (A_prod * t * Gamma) .* (R .* (1 - Gamma)).^ (1 - t);
    
    Pi = @ (k , P) (k <= k_bar(P)) .* (...
        (1 - t) * (1 - Gamma)^t * A_prod * P  .* (k_bar(P)).^(t) ...
        + (W * (1-Gamma) / Gamma - E) * Utilization_low_region * k ) + ...
        Utilization_low_region .* ((k > k_bar(P)) .* (k <= k_hat(P)).*  (...
        A_prod * (1 - Gamma)^t * P .* k .^t - E * k) +...
        (k > k_hat(P)).* (A_prod * (1 - Gamma)^t * P .* k_hat(P) .^t - E * k_hat(P))) + ...
        (1 - Utilization_low_region) .* (k > k_bar(P)) .* ( ...
        (1 - t) * (1 - Gamma)^t * A_prod * P  .* (k_bar(P)).^(t) ...
        + (W * (1-Gamma) / Gamma - E) * Utilization_low_region * k );
    
    Pi_prime = @ (k , P) (k <= k_bar(P)) .* (...
         (W * (1-Gamma) / Gamma - E) * Utilization_low_region) + ...
        Utilization_low_region .* ((k > k_bar(P)) .* (k <= k_hat(P)).*  (...
        A_prod * (1 - Gamma)^t * t * P .* k .^(t - 1) - E) );
    
    if E==0
        % Here k_hat(P) is Inf, so we are never in the highest region.
        % Labor savings cannot turn negative here!
        % Utilization always set to u=1 because of no flow cost.
        Pi = @ (k , P) (k <= k_bar(P)) .* (...
        (1 - t) * (1 - Gamma)^t * A_prod * P  .* (k_bar(P)).^(t) ...
        + (W * (1-Gamma) / Gamma - E) * k ) + ...
        (k > k_bar(P)) .* (k <= k_hat(P)).*  (...
        A_prod * (1 - Gamma)^t * P .* k .^t - E * k);
    
    end
    
end


% Adjustment cost function
if pay_for_delta == 0
    phi_k = @(I_k) psi_k/2 .* (I_k) .^2;
else
    if kinked
        phi_k = @(I_k) psi_k/2 .* (I_k) .^2 + abs(I_k) .* chi_k;
    else
        phi_k = @(I_k) psi_k/2 .* (I_k) .^2;
    end
end

% Derivative
if pay_for_delta == 0
    dphi_k = @(I_k) psi_k .* (I_k) ;
else
    if kinked
        dphi_k = @(I_k) psi_k .* (I_k) + sign(I_k) .* chi_k;
    else
        dphi_k = @(I_k) psi_k .* (I_k) ;
    end
end

% Inverse derivative
if pay_for_delta == 0
    inv_dphi_k = @(y) 1/psi_k .* y;
else
    if kinked
        inv_dphi_k = @(y) 1/psi_k .* ( max(y - chi_k, 0) + min(y + chi_k, 0));
    else
        inv_dphi_k = @(y) 1/psi_k .* y;
    end
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


%% Grid parameters

% Size
N_k = params.N_k; % capital

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
        out.outDiffusion = outDiffusion;
        p_grid = outDiffusion.p_grid;
        pmin = min(p_grid);
        pmax = max(p_grid);
        % standardize so that trapz(p_grid,pi0) = 1
        pi0 = outDiffusion.g_0./outDiffusion.dp_mean';
        Lambda = outDiffusion.Lambda;
    end
    
elseif  strcmp(ShockType, 'GBM')    
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

    outDiffusion = DiscretizeGBM(settings);
    out.outDiffusion = outDiffusion;
    p_grid = outDiffusion.p_grid;
    pmin = min(p_grid);
    pmax = max(p_grid);
    % standardize so that trapz(p_grid,pi0) = 1
    pi0 = outDiffusion.g_0./outDiffusion.dp_mean';
    Lambda = outDiffusion.Lambda;

else
    error('Specify ShockType, accepts Markovian, Lognormal, Diffusion or GBM')    
end   

if ~(strcmp(ShockType, 'Diffusion')||strcmp(ShockType, 'GBM'))
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

    %% Boundaries
    % labor share


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


% Capital
if isfield(params, 'kmin')
    kmin = params.kmin;
else
    kmin = 0;
end

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
%% Convergence parameters
Delta = 1000;
tol = 1e-6;
maxit = 20; % maximum value function iterations
%% grids

% option to have unevenly spaced grid Gamma > 1 gives more points for lower
%K
if isfield(params, 'GridKGamma')
    GridKGamma = params.GridKGamma;
else
    GridKGamma = 1;
    out.params.GridKGamma = GridKGamma;
end


% grids
k_grid = kmin + (kmax - kmin) * (linspace(0, 1, N_k).^GridKGamma)';


% in single vector form with p outer state, then a, then k
kk = repmat(k_grid, N_p, 1);
pp = repmat (p_grid , N_k , 1);
pp = reshape( pp, [ N_k * N_p 1]);


%DELTAS for k
dk_vec_up = diff(k_grid);
dk_vec_up = [dk_vec_up; dk_vec_up(end)];
dk_vec_dwn = diff(k_grid) ;
dk_vec_dwn = [dk_vec_dwn(1); dk_vec_dwn];

dk_vec_up = kron(ones(length(p_grid),1), dk_vec_up);
dk_vec_dwn = kron(ones(length(p_grid),1), dk_vec_dwn);

dk_mean = (dk_vec_up + dk_vec_dwn)/2;
dk_mean(1) = dk_mean(1)/2;
dk_mean(end) = dk_mean(end)/2;


if stochastic
    dp_vec_up = diff(p_grid');
    dp_vec_up = [dp_vec_up; dp_vec_up(end)];
    dp_vec_dwn = diff(p_grid') ;
    dp_vec_dwn = [dp_vec_dwn(1); dp_vec_dwn] ;

    dp_vec_up = kron(dp_vec_up, ones(length(k_grid),1));
    dp_vec_dwn = kron(dp_vec_dwn, ones(length(k_grid),1));

    dp_mean = (dp_vec_up + dp_vec_dwn)/2;
    dp_mean(1) = dp_mean(1)/2;
    dp_mean(end) = dp_mean(end)/2;
end


%% DEFINE TRANSITIONS

% create indices (and implicitly set boundary conditions)
ind = (1:1:N_k * N_p)'; % starting state indices

vector_ind = ind2sub_vec([N_k N_p], ind); %matrix with indices for all states
k_up = (N_k  - vector_ind(:,1)) > 0 ; %positions where it is possible to go up
k_dwn = ( vector_ind(:,1) -1) > 0 ; %positions where it is possible to go down


%transition indices if drift is positive (up) or negative (dwn)
% set indices so that if we cannot move up, the "up" index is the same as
% the original index, similar for down
ind_up_k = sub2ind_vec([N_k  N_p], vector_ind + [k_up zeros(length(ind),1)]);
ind_dwn_k = sub2ind_vec([N_k  N_p], vector_ind - [k_dwn zeros(length(ind),1)]);


%% STOCHASTIC PROCESSES
if ~(strcmp(ShockType, 'Diffusion')||strcmp(ShockType, 'GBM'))

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
         shift_possible = vector_ind(: , 2) + shift_state <= N_p; % indicator u-bound
         ind_col_up = sub2ind_vec([N_k  N_p],...
             vector_ind + [zeros(length(ind),1) shift_possible * shift_state ]);
         % select only states where it actually goes up
         ind_col_up = ind_col_up(logical(shift_possible));
         ind_row_up = ind(logical(shift_possible));

         ind_up_p = [ind_up_p;  ind_col_up];
         ind_p = [ind_p;  ind_row_up];

         for state = 1: (N_p - shift_state) %loop for states where shift_state is possible
             % get off-diag. elms that are above current state
             F_p = [F_p; repmat( Lambda(state, state + shift_state),N_k, 1)];
         end
     end
     %% JUMPS DOWN
    for shift_state = 1 : N_p - 1 %over max number of possible jumps down
         shift_possible = vector_ind(: , 2) - shift_state >= 1; % indicator l-bound
         ind_col_dwn = sub2ind_vec([N_k  N_p],...
             vector_ind - [zeros(length(ind),1) shift_possible * shift_state ]);
         % select only states where it actually goes up
         ind_col_dwn = ind_col_dwn(logical(shift_possible));
         ind_row_dwn = ind(logical(shift_possible));

         ind_dwn_p = [ind_dwn_p;  ind_col_dwn];
         ind_p = [ind_p;  ind_row_dwn];

         for state = (shift_state + 1): N_p %loop for states where jumps possible
             % get off-diag. elms that are below current state
             B_p = [B_p; repmat( Lambda(state, state - shift_state), N_k, 1)];
         end
     end 
    %% DIAGONAL
    if stochastic
         for state = 1: N_p
             DIAG_p = [DIAG_p; repmat( Lambda(state, state), N_k, 1)];
         end
    else
        DIAG_p = zeros( N_k, 1);
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
    F_p = repmat(outDiffusion.F_p, N_k, 1);
    F_p = F_p(:);
    B_p = repmat(outDiffusion.B_p, N_k, 1);
    B_p = B_p(:);
    DIAG_p = - F_p - B_p;
    else
         DIAG_p = zeros( N_k, 1);
         F_p = DIAG_p;
         B_p = DIAG_p;
    end
    
    
    p_up = (N_p  - vector_ind(:,2)) > 0 ; %positions where it is possible to go up
    p_dwn = ( vector_ind(:,2) -1) > 0 ; %positions where it is possible to go down


    %transition indices if drift is positive (up) or negative (dwn)
    % set indices so that if we cannot move up, the "up" index is the same as
    % the original index, similar for down
    ind_up_p = sub2ind_vec([N_k  N_p], vector_ind + [zeros(length(ind),1) p_up]);
    ind_dwn_p = sub2ind_vec([N_k  N_p], vector_ind - [zeros(length(ind),1) p_dwn]);
    ind_p = repmat(ind, 2, 1);
    


end
%% ALLOCATION
% set up the sparse rows and columns for the matrices
sparse_rows = [repmat(ind, 3 ,1); ind_p];
sparse_cols = [ind; ind_up_k;  ind_dwn_k; ...
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

Pi_vec = Pi(kk, pp);
% Initialize the value function: assume no positive drift in initial states,
% value is just current profit minus upkeep of stocks


%% Value Function Guess

if isfield(params, 'ValueGuess')
    ValueGuess = params.ValueGuess;
else
    ValueGuess = (Pi_vec - phi_k(delta.*kk))./rho;
    out.params.ValueGuess = ValueGuess;
end

V_next = ValueGuess;


supnorm = zeros(maxit, 1);

for J = 1:maxit
    
    %% FORWARD DIFFERENCES
    V = V_next;
    % Forward differences: get a zero when the step forward is not possible
    dVF_k = (V(ind_up_k) - V(ind))./dk_vec_up;
    dVB_k = - (V(ind_dwn_k) - V(ind))./dk_vec_dwn;
    

    % Check it is correct
    % [kk V (V(ind_up_k) - V(ind))]
    % [aa V (V(ind_up_a) - V(ind))]
    % [pp V (V(ind_up_p) - V(ind))]

    %% STATE CONSTRAINTS

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
    end

    % Drifts
    s_F_k = i_F_k - delta .* kk;
    s_B_k = i_B_k - delta .* kk;

    % instantaneous costs of investment
    if pay_for_delta == 1
        Phi_F_k =  phi_k(i_F_k) + R * i_F_k;
        Phi_B_k =  phi_k(i_B_k) + R * i_B_k;
        Phi_0_k =  phi_k(delta * kk) + R * (delta * kk);
    else
        Phi_F_k =  phi_k(i_F_k - delta .* kk) + R * i_F_k;
        Phi_B_k =  phi_k(i_B_k - delta .* kk) + R * i_B_k;
        Phi_0_k =  phi_k(i_0_k - delta .* kk) + R * i_0_k;
    end

    % indicators
    I_F_k = s_F_k > 0;
    I_B_k = s_B_k < 0;
    I_0_k = (1 - I_F_k - I_B_k);
    %I_0_a = (1 - I_F_a - I_B_a);

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
        u = Pi_vec + p_y * ( - I_F_k .* Phi_F_k - I_B_k .* Phi_B_k - I_0_k .* Phi_0_k);
    else
        u = Pi_vec + ( - I_F_k .* Phi_F_k - I_B_k .* Phi_B_k - I_0_k .* Phi_0_k);
    end
        


    % Construct the sparse matrix
    F_k = I_F_k .* s_F_k ./ dk_vec_up; % positive drifts that are actually positive
    B_k = I_B_k .* s_B_k ./ dk_vec_dwn; 

    DIAG = DIAG_p - F_k + B_k ;
    % Populate A

    A = sparse(sparse_rows, sparse_cols, [DIAG ;...
        F_k; -B_k; F_p; B_p]);
    
    %% UPDATE VALUE
    B = (1/Delta + rho)*speye(N_p * N_k) - A;
    b = u + V/Delta;
    V_next = B\b;
    supnorm(J) = max(abs(V - V_next) ./ V_next * 100);
    if print_conv
        disp(['Iteration ' num2str(J) ' completed'])
        disp(['Convergence : supnorm ' num2str(supnorm(J))])
    end
    if supnorm(J) < tol   
        break
    end
    
%     disp(['ho fatto un iterazione della vf'])

end

% disp(['ho finito la value function.'])




%% Compute Investment Cutoffs by state
dV = dVF_k .* I_F_k + I_B_k .* dVB_k + .5 .* I_0_k .* (dVF_k + dVB_k);


%Investment policy for each firm
i_k_pol = I_F_k .* i_F_k + I_B_k .* i_B_k + I_0_k * delta .* kk;


%increment to use for Inv_Cutoffs
%dk_dV = dk_vec_up .* I_F_k + dk_vec_dwn .* dVB_k + .5 .* I_0_k .* (dk_vec_up + dk_vec_dwn);

if pay_for_delta == 1 %find the R_star where the inv. covers depreciation
    dV_reshaped = reshape(inv_dphi_k(dV - R) - delta * kk, [N_k, N_p]);
else %find where the drift equals 0 without paying for depreciation
    dV_reshaped = reshape(inv_dphi_k(dV - R), [N_k, N_p]);
end


if AdjAtEntry % reshape value function for entrants: They max V(R) - Costs(R)
    dV_entrant = reshape(inv_dphi_k(dV - R) -  kk, [N_k, N_p]);
else
    dV_entrant = reshape(dV - R, [N_k, N_p]);
end

out.INV_CUTOFFS = zeros(N_p,1);
out.INV_CUTOFF_MAX_VALUE = out.INV_CUTOFFS;
% Linearly interpolate the investment condition to find the robots level
% such that the firm keeps a constant size. (perfectly offset depreciation)
idx_entry_cutoffs = zeros(N_p,1);
idx_inv_cutoffs = zeros(N_p,1);
for idx = 1:N_p
    idx_max_value = find(dV_entrant(:,idx) < 0, 1);
    idx_cutoff = find(dV_reshaped(:,idx) < 0, 1);
    if idx_cutoff > 1
        idx_cutoff_lin = sub2ind([N_k, N_p], idx_cutoff, idx);
        out.INV_CUTOFFS(idx) = kk(idx_cutoff_lin - 1) - ...
           (kk(idx_cutoff_lin) - kk(idx_cutoff_lin - 1)) * dV_reshaped(idx_cutoff-1,idx) / ...
                          (dV_reshaped(idx_cutoff,idx) - dV_reshaped(idx_cutoff-1,idx));
    end
    if idx_max_value > 1
        idx_cutoff_lin = sub2ind([N_k, N_p], idx_max_value, idx);
        out.INV_CUTOFF_MAX_VALUE(idx) = kk(idx_cutoff_lin - 1) - ...
           (kk(idx_cutoff_lin) - kk(idx_cutoff_lin - 1)) * dV_entrant(idx_max_value-1,idx) / ...
                          (dV_entrant(idx_max_value,idx) - dV_entrant(idx_max_value-1,idx));
    end
    
    if isempty(idx_cutoff)
        idx_inv_cutoffs(idx) = N_k;
    else
        idx_inv_cutoffs(idx) = idx_cutoff;
    end
    if isempty(idx_max_value)
        idx_entry_cutoffs(idx) = N_k;
    else
        idx_entry_cutoffs(idx) = idx_max_value;
    end
end
out.idx_inv_cutoffs = idx_inv_cutoffs;
out.idx_entry_cutoffs = idx_entry_cutoffs;
% if ~ReducedOutput
%     out.dV = dV_reshaped;
% end

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

%% STATIONARY DISTRIBUTION STARTING FROM g_0;
if ~strcmp(ShockType, 'GBM')

    AT = A';

else
   % step 1: preserve existing entries of the A matrix encoding the
   % policies and the revenue drifts conditional on survival
   AT_survive = A';
   [row, col, val] = find(AT_survive);
   
   if strcmp(settings.ResetType,'Deterministic')
       % step 2: find the index for the reset value of p and associated k_star
       p_reset_idx = find(p_grid >= params.settings.reset_val,1);
       if strcmp(ReEntryK, 'Optimal')
           k_reset_idx = idx_entry_cutoffs(p_reset_idx);
       elseif strcmp(ReEntryK, 'Kstar')
           % Step 2: extract reset indices
            k_reset_idx = idx_inv_cutoffs(p_reset_idx);
       elseif strcmp(ReEntryK, 'Zero')
            k_reset_idx = ones(length(p_reset_idx),1);
       else
           error('Supported ReEntryK: Optimal, Kstar OR Zero')
       end
       % reset index in vectorized kk, pp
       reset_idx = sub2ind_vec([N_k N_p], [k_reset_idx p_reset_idx]);
        row = [row; (1:N_k*N_p)'; repmat(reset_idx, N_k*N_p, 1)];
        col = [col; (1:N_k*N_p)'; (1:N_k*N_p)'];
        val = [val; -repmat(params.settings.reset_rate, N_k*N_p, 1);...
            repmat(params.settings.reset_rate, N_k*N_p, 1)];

   elseif strcmp(settings.ResetType,'LogNormal')
       reset_rates_logn = outDiffusion.reset_rates_logn;
       p_reset_idx = outDiffusion.in_indexes';
       if strcmp(ReEntryK, 'Optimal')
            k_reset_idx = idx_entry_cutoffs(p_reset_idx);
       elseif strcmp(ReEntryK, 'Kstar')
            k_reset_idx = idx_inv_cutoffs(p_reset_idx);
       elseif strcmp(ReEntryK, 'Zero')
            k_reset_idx = ones(length(p_reset_idx),1);
       else
           error('Supported ReEntryK: Optimal, Kstar OR Zero')
       end
        reset_idx = sub2ind_vec([N_k N_p], [k_reset_idx p_reset_idx]);
        row = [row; (1:N_k*N_p)'; repmat(reset_idx, N_k*N_p, 1)];
        col = [col; (1:N_k*N_p)'; kron((1:N_k*N_p)', ones(length(reset_idx), 1))];
        val = [val; -repmat(params.settings.reset_rate, N_k*N_p, 1); ...
        kron(ones(N_k*N_p, 1), params.settings.reset_rate * reset_rates_logn') ];
       
   else
       error('Wrong settings.ResetType')
   end
    AT = sparse(row, col, val, N_k*N_p, N_k*N_p);
    out.AT = AT;
end
 

if isfield(params, 'DT')
    DT = params.DT;
else
    DT = 1000;
    out.params.DT = DT;
end

if isfield(params, 'g_0')
    g_0 = params.g_0;
else
    g_0 = zeros(N_k*N_p, 1);
    if (strcmp(ShockType, 'Diffusion') ||strcmp(ShockType, 'GBM')) 
        g_0(1:N_k:N_k*N_p) = pi0.*outDiffusion.dp_mean'; % start everyone at 0 capital on statdist of p
    else
        g_0(1:N_k:N_k*N_p) = pi0;
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
    tolG = params.tolG;
else
    tolG = 1e-6;
    out.params.tolG = tolG;
end


II = speye(size(AT));
TransMat = plus(II, - DT*AT);

gg = ones(size(g_0));
supnorm = 1;
while supnorm > tolG
    gg = TransMat \  g_0;
    supnorm = max(abs(gg - g_0));
    g_0 = gg;
end

% disp('ho finito la KFE')

%rescale to obtain pdf
if (strcmp(ShockType, 'Diffusion')||strcmp(ShockType, 'GBM'))
    gg = gg./(dp_mean.*dk_mean); 
else
    gg = gg./(dk_mean);
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


%% POLICIES and STATISTICS

%% DEFINE INTEGRALS
if (strcmp(ShockType, 'Diffusion')||strcmp(ShockType, 'GBM'))
    integrate_kp = @(fun) sum(sum(fun.*dk_mean.*dp_mean));
else
    integrate_kp = @(fun) sum(sum(fun.*dk_mean));
end

%% LABOR DEMAND



if isfield(params, 'ReducedOutput')
    ReducedOutput = params.ReducedOutput;
else
    ReducedOutput = false;
end


% Labor demand
% for each firm
indL = ((k_bar(pp) - kk)>=0) .*  ((k_bar(pp) - kk)) *(1 - Gamma)/ Gamma;
% overall
L = integrate_kp(gg .* indL);

% If workers produce robots, add robot purchases to labor demand
if strcmp(RobotsProducers, 'Workers')
    % Assume production function is R = A*L Therefore add I_k/A to L^d
    L = L + integrate_kp(gg .* i_k_pol/A_robots);
end

%compute the mass of firms fully automating
gg_auto = gg;
gg_auto(kk<k_bar(pp)) = 0;

out.Mass_auto = integrate_kp(gg_auto);
if ~ReducedOutput
    out.gg_auto = gg_auto;
end

% compute measure of automation in eqm
out.auto_extent = integrate_kp(gg .* min(kk ./ k_bar(pp), 1));

out.kmaxgg = max(kk(gg>tolG));

% Optional output
if ~ReducedOutput
    out.statDist = gg;
    out.A = A;
    out.V = V_next;
end
% Decision to automate or not: derivative of value function at R=0 in
% maximum productivity state (using only Forward Difference because at the
% boundary) must be positive! 
% out.INV_TRADEOFF_MAX_P = dVF_k(N_k*(N_p - 1) + 1) - R;

%% OUT GE OBJECTS

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
%% OLD GE
%{

% REVENUES NET OF INVESTMENT COSTS IF PAID IN INTERMEDIATES (DEFAULT)
if ~strcmp(InvGoodCost, 'Final') && ~strcmp(InvGoodCost, 'External')
    netRev = IndRev - IndAdjCosts; % DEFAULT IS THISSSSSSSS!!!!!!!!!!!!!!!!!!!!
%     netRev = IndRev / A_prod - IndAdjCosts./A_prod; % DEFAULT IS THISSSSSSSS!!!!!!!!!!!!!!!!!!!!
else %costs in final good, so the convention on Q is that it is gross of costs
    netRev = InqdRev;
end

if ~strcmp(EnergyCosts,'Final') && ~strcmp(InvGoodCost, 'External') %Energy costs are paid in intermediates

    %compute intermediate good quantities and profits in the various scenarios
    
    % DEFAULT IS THIS!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    if strcmp(RobotsProducers,'None') % production of robots is paid in intermediates: Net output excludes robot purchases
        indQ = netRev - R * i_k_pol - E * kk; 
%         indQ = netRev - (R * i_k_pol + E * kk) ./ A_prod; 
        indPi = indQ - W * indL;
%         indPi = A_prod * indQ - W * indL;
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

%}
%% GE    
% NEW DEFAULT: Robots are produced using the final good (HENCE DIVIDE Q by A_PROD!!!!)
indQ = IndRev;
indPi = IndRev - IndAdjCosts - R * i_k_pol - E * kk - W * indL;
netRev = IndRev - IndAdjCosts;

% Compute Useful GE objects
out.GE.Pi = integrate_kp(gg .* indPi); % Aggregate Profits (less aggregate adjustment costs)
out.GE.q = integrate_kp(gg .* indQ)./A_prod; % Aggregate Output
% NOTE FOR GE.q: note that A_prod is treated as the price of the good here,
% and adjustment costs are paid in unit of the good produced.

if ~ReducedOutput
    out.GE.IndRev = IndRev;
    out.GE.IndAdjCosts = IndAdjCosts;
    out.GE.AggAdjCosts = integrate_kp(gg .*IndAdjCosts);
    out.GE.netRev = netRev;
end
out.GE.AggnetRev = integrate_kp(gg .*netRev);
out.GE.AggRev = integrate_kp(gg .* IndRev);
out.GE.AggECosts = integrate_kp(gg .* E .* kk);
if ~ReducedOutput
    out.GE.indQ = indQ;
    out.GE.indPi = indPi;
end
%% OUT OTHER PARAMETERS
if ~ReducedOutput
    out.i_k_pol = i_k_pol;
    out.kk = kk;
    out.pp = pp;
    out.L_pol = indL;
end
out.k_grid = k_grid;
out.p_grid = p_grid;
out.W = W;



if ~ReducedOutput
    out.u = u;
end
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
out.LS_AggRev = W * L / out.GE.AggRev;
out.LS_q = W * L / out.GE.q;
out.YL = out.GE.AggRev / L;


if ReducedOutput
    out.params = rmfield(out.params, 'ValueGuess');
    out.params = rmfield(out.params, 'g_0');
    out.params = rmfield(out.params, 'Lambda');
end

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

if (strcmp(ShockType, 'Diffusion')||strcmp(ShockType, 'GBM'))
    out.pi0_int = pi0.*outDiffusion.dp_mean';
else
    out.pi0_int = pi0;
end


out.params.sigma_P = sqrt(out.p_grid.^2*out.pi0_int-(out.p_grid*out.pi0_int)^2);

% 
if ~ReducedOutput
    out.dV_Fk = dVF_k;
    out.dV_Bk = dVB_k;
end
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





    