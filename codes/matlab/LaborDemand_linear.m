% Fornino Manera (2019)
% 
% DATE: January 24, 2021
%
% Project: Automation and the Future of Work: Assessing the Role of Labor
%          Flexibility

function [L, out] = LaborDemand_linear(W, R, params)
% This function finds the solution to the firm problem as well as the
% stationary distribution for given parameters. The model is the one with
% linear adjustment costs.
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
psi_plus = params.psi_plus; % linear adjustment cost up
psi_minus = params.psi_minus; % linear adjustment cost down

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
kmin = 0;

if isfield(params, 'kmax')
    kmax = params.kmax;
else
    kmax = ((t * A_prod * (1 - Gamma)^t * pmax)/(E + delta * (psi_plus + R))).^(1/(1-t));
end

out.params.kmax = kmax;
%% Convergence parameters
Delta = 1e3;
tol = 1e-4;
maxit = 2000; % maximum value function iterations
%% grids

%option to have unevenly spaced grid Gamma > 1 gives more points for lower
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

%DELTAS for p

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
%     ValueGuess = (Pi_vec - (psi_plus + R) * delta* kk)./(rho);
    ValueGuess = Pi_vec/rho;
%     ValueGuess = sqrt(kk.*pp);
    out.params.ValueGuess = ValueGuess;
end

V_next = ValueGuess;


supnorm = zeros(maxit, 1);
i_inv = zeros(N_p,1);
i_disinv = i_inv;
R_inv = i_inv;
R_disinv = R_inv;



supnorms = nan(50,1);



for J = 1:maxit
    
    %% FORWARD DIFFERENCES
    V = V_next;
    % Forward differences: get a zero when the step forward is not possible
    dVF_k = (V(ind_up_k) - V(ind))./dk_vec_up;
    dVB_k = - (V(ind_dwn_k) - V(ind))./dk_vec_dwn;
    dVF_k(ind(ind_up_k == ind)) = dVF_k(ind(ind_up_k == ind) - 1);
    dVB_k(ind(ind_dwn_k == ind)) = dVB_k(ind(ind_dwn_k == ind) + 1);
    dVC_k = (dVF_k + dVB_k)/2;
    F_p = repmat(outDiffusion.F_p, N_k, 1);
    F_p = F_p(:);
    B_p = repmat(outDiffusion.B_p, N_k, 1);
    B_p = B_p(:);


    %% STATE CONSTRAINTS
    % preallocation of sparse rows and sparse cols
    sparse_rows = [repmat(ind, 3, 1); ind_p];
    sparse_cols = [repmat(ind, 3, 1);  ind_up_p;  ind_dwn_p] ;
    vals = [delta * kk; - delta * kk; zeros(3*N_k*N_p , 1)];
    vals_transpose = -ones(N_k*N_p,1);
    u = Pi_vec;
    rho_delta_vec = ones(N_k*N_p,1);
    delta_vec = zeros(N_k*N_p,1);
    inaction_ind = zeros(N_k*N_p,1);
    i_inv_vec = zeros(N_p, 2);
    % Impose state constraints to correctly compute the optimal controls.
    % Set the up_drift at 0 at upper boundary, so that 
    % k_dot = (I - delta)*k -> I = delta
    % ind_up_k == ind shows where upward movement is not possible
    
    for jp = 1:N_p
        
%         i_inv = (jp - 1) * N_k + find(dVC_k((1 + (jp - 1) * N_k): jp * N_k ) <= (R + psi_plus), 1);
%         i_disinv = (jp - 1) * N_k + find(dVC_k((1 + (jp - 1) * N_k): jp * N_k ) <= (R - psi_minus), 1);
        i_inv = (jp - 1) * N_k + find(dVF_k((1 + (jp - 1) * N_k): jp * N_k ) < (R + psi_plus) - (psi_plus + psi_minus)/1000, 1);
        i_disinv = (jp - 1) * N_k + find(dVB_k((1 + (jp - 1) * N_k): jp * N_k ) < (R - psi_minus) , 1);
         
        
        if isempty(i_inv)
            error(['No investment cutoff found at jp = ' num2str(jp)])
            %i_inv = (jp - 1) * N_k + 1;
        end
        
        if isempty(i_disinv)
            subplot(1,2,1)
            dVC_k((1 + (jp - 1) * N_k): jp * N_k )
            plot(dVC_k((1 + (jp - 1) * N_k): jp * N_k ))
            hold on
            hline((R - psi_minus))
            
            dVB_k((1 + (jp - 1) * N_k): jp * N_k )
            plot(dVB_k((1 + (jp - 1) * N_k): jp * N_k ))
            hold off
            
            subplot(1,2,2)
            plot(V((1 + (jp - 1) * N_k): jp * N_k ))
            error(['No disinvestment cutoff found at jp = ' num2str(jp)])
        end
        
            
        % impose that R_inv = R(i_inv - 1) and R_disinv = R(i_disinv)
        % now i_inv - 1 and i_disinv are both inside the inaction area
        
        if i_inv == (1 + (jp - 1) * N_k) % first index is returned
            R_inv(jp) = kk(1 + (jp - 1) * N_k);
        else
            R_inv(jp) = kk(i_inv - 1);
        end
        
        R_disinv(jp) = kk(i_disinv);
        
        
        % Fill in the "left" column
        
        % investment
        if (i_inv - 1) > (1 + (jp - 1) * N_k)
            sparse_cols((1 + (jp - 1) * N_k) : i_inv - 2) = i_inv - 1;
            vals((1 + (jp - 1) * N_k) : i_inv - 2) = 1;
        elseif (i_inv - 1) == (1 + (jp - 1) * N_k)
            vals((1 + (jp - 1) * N_k):(i_inv - 1)) = 0;
        else
            vals((1 + (jp - 1) * N_k):(i_inv)) = 0;
        end
        
        % Inaction Region : vals is - delta * kk
        sparse_cols(max(i_inv - 1, (1 + (jp - 1) * N_k))  : i_disinv) =...
            ind_dwn_k(max(i_inv - 1, (1 + (jp - 1) * N_k))  : i_disinv);
        
        % Disinvestment
        sparse_cols(i_disinv + 1 : jp * N_k) = i_disinv;
        vals(i_disinv + 1 : jp * N_k) = 1;
        
        

        % Fill in the "right" column
        
        vals( N_k * N_p +  (1 + (jp - 1) * N_k): ...
             N_k * N_p  + max(i_inv - 2, (1 + (jp - 1) * N_k))) = 0;
         
        vals( N_k * N_p + min(i_disinv + 1, jp * N_k): ...
             N_k * N_p  + jp * N_k) = 0;
         
         % No inaction is needed as sparse_cols has diag as default and
         % vals has delta * kk
        
         
         
        % Manage Ito Terms: they should be zero outside of inaction region,
        % i.e. below i_inv - 1 and above i_disinv
        if jp < N_p
            F_p( (jp - 1) * N_k + 1 : i_inv - 2) = 0;
            F_p(i_disinv + 1 : N_k * jp) = 0;
        end
        
        if jp > 1
            B_p( (jp - 1) * N_k + 1 : i_inv - 2) = 0;
            B_p(i_disinv + 1 : N_k * jp) = 0; 
        end
        
        
        % instantaneous utility is Pi by default, and resale value of
        % capital when outside inaction region, below inv - 1 and above
        % i_disinv
        if i_inv - 2 >= (jp - 1) * N_k + 1
            u((jp - 1) * N_k + 1: i_inv - 2) = ...
            - (psi_plus + R) * (R_inv(jp) - kk((jp - 1) * N_k + 1: i_inv - 2));
        end
        
        if i_disinv + 1 <= jp * N_k
            u(i_disinv + 1 : jp * N_k) = ...
                - (psi_minus - R) * (kk(i_disinv + 1: jp * N_k) -   R_disinv(jp) );
        end

        
        rho_delta_vec(max(i_inv - 1, (jp - 1) * N_k + 1):i_disinv) = 1/Delta + rho;
        delta_vec(max(i_inv - 1, (jp - 1) * N_k + 1):i_disinv) = 1/Delta;
        inaction_ind(max(i_inv - 1, (jp - 1) * N_k + 1):i_disinv) = 1;
        
        i_inv_vec(jp, 1) = i_inv - 1 - (jp - 1) * N_k ;
        i_inv_vec(jp, 2) = i_disinv - (jp - 1) * N_k ;

        vals_transpose(max(1 + (jp-1)*N_k, i_inv - 1): i_disinv) = 0;
        
    end
    
    DIAG_p = - F_p - B_p;
    vals(2 * N_k * N_p + 1:end) = [DIAG_p ; F_p ; B_p];

    % Populate A

    A = sparse(sparse_rows, sparse_cols, vals);

    
    %% UPDATE VALUE
    B = sparse(ind, ind, rho_delta_vec) - A;
    b = u + V.* delta_vec;
     V_next = B\b;
%      V_next =  inaction_ind .* (V + Delta .* ( u + A * V - rho * V)) + (1 - inaction_ind) .* (A * V + u);
    supnorm(J) = max(abs(V - V_next));
    
%     supnorms = [supnorms(2:50); supnorm(J)];
% 
%     plot(supnorms)
%     hold on
%     hline(tol)
%     hold off
%     drawnow()
%     mesh(reshape(V_next,N_k,N_p))
%     pause(.01);
    if any(isnan(V_next))
        error('NaN in V')
    end
%     if print_conv
%         disp(['Iteration ' num2str(J) ' completed'])
%         disp(['Convergence : supnorm ' num2str(supnorm(J))])
%     end
    if supnorm(J) < tol   
        break
    end

end


    
    AT = transpose(sparse(vertcat(sparse_rows, ind),...
                          vertcat(sparse_cols, ind),...
                          vertcat(vals, vals_transpose)));



%% STATIONARY DISTRIBUTION STARTING FROM g_0;
% AT = A';
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
    if strcmp(ShockType, 'Diffusion') 
        g_0(1:N_k:N_k*N_p) = pi0.*outDiffusion.dp_mean'; % start everyone at 0 capital on statdist of p
        g_0(find(inaction_ind == true, 1)) = 1;
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

 while max(abs(gg-g_0)) > tolG
gg = TransMat \  g_0;
g_0 = gg;
 end

%rescale to obtain pdf
if strcmp(ShockType, 'Diffusion')
    gg = gg./(dp_mean.*dk_mean); 
else
    gg = gg./(dk_mean);
end
 
out.V = V_next;
out.kk = kk;
out.k_grid = k_grid;
out.p_grid = p_grid;
out.pp = pp;
out.inaction_indices = i_inv_vec;
out.inaction_region = inaction_ind;
L = 1;
out.gg = gg; 



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
    integrate_kp = @(fun) sum(sum(fun.*dk_mean.*dp_mean));
else
    integrate_kp = @(fun) sum(sum(fun.*dk_mean));
end



%% POLICIES and STATISTICS


if isfield(params, 'ReducedOutput')
    ReducedOutput = params.ReducedOutput;
else
    ReducedOutput = false;
end

%Investment policy for each firm
% i_k_pol = I_F_k .* i_F_k + I_B_k .* i_B_k + I_0_k * delta .* kk;
% i_k_pol = I_F_k .* i_F_k + I_B_k .* i_B_k + I_0_k * i_0_k;

% Labor demand
% for each firm
indL = ((k_bar(pp) - kk)>=0) .*  ((k_bar(pp) - kk)) *(1 - Gamma)/ Gamma;
% overall
L = integrate_kp(gg .* indL);
% 
% % If workers produce robots, add robot purchases to labor demand
% if strcmp(RobotsProducers, 'Workers')
%     % Assume production function is R = A*L Therefore add I_k/A to L^d
%     L = L + integrate_kp(gg .* i_k_pol/A_robots);
% end

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
% 
% % Optional output
% if ~ReducedOutput
%     out.statDist = gg;
%     out.A = A;
%     out.V = V_next;
% end
% Decision to automate or not: derivative of value function at R=0 in
% maximum productivity state (using only Forward Difference because at the
% boundary) must be positive! 
% out.INV_TRADEOFF_MAX_P = dVF_k(N_k*(N_p - 1) + 1) - R;

% Compute Investment Cutoffs by state
% Note: correct approximation of V' must use UPWIND SCHEME!!!

%{

if ~ReducedOutput
    dV = dVF_k .* I_F_k + I_B_k .* dVB_k + .5 .* I_0_k .* (dVF_k + dVB_k);
    %increment to use for Inv_Cutoffs
    %dk_dV = dk_vec_up .* I_F_k + dk_vec_dwn .* dVB_k + .5 .* I_0_k .* (dk_vec_up + dk_vec_dwn);

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
end
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
%}

%% OUT USEFUL FUNCTION HANDLES
out.handles.Pi = Pi;
out.handles.R_bar = k_bar;
if strcmp(Utilization,'Yes')
    out.handles.R_hat = k_hat;
    out.handles.Pi_prime = Pi_prime;
end
out.handles.p_bar = p_bar;
out.handles.p_hat = p_hat;

% out.handles.integrate_kp = integrate_kp;

% if strcmp(ShockType, 'Diffusion')
%     out.pi0_int = pi0.*outDiffusion.dp_mean';
% else
%     out.pi0_int = pi0;
% end


% out.params.sigma_P = sqrt(out.p_grid.^2*out.pi0_int-(out.p_grid*out.pi0_int)^2);

% 
% if ~ReducedOutput
%     out.dV_Fk = dVF_k;
%     out.dV_Bk = dVB_k;
% end
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
%}



end





    