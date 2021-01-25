% Fornino Manera (2019)
% 
% DATE: January 24, 2021
%
% Project: Automation and the Future of Work: Assessing the Role of Labor
%          Flexibility

function out = GetGrids_rigidlabor(params, ReducedOutput)
% This subroutine is called by SolveTransition.m to compute the grid given
% the initial and the final steady states. Then, the upper envelope of the
% two is chosen for the actual grid used along the transition. This function
% applies to the case of labor adjustment costs.
%
% Input
% - V: value function
% - params: parameters of main model
%
% Outputs("out" struct object, same as LaborDemand_trapz.m)
% - out: struct object containing a host of information on the transitional
%        dynamics

%% Set  Parameter Values: 

if isfield(params, 'A_prod')
    A_prod = params.A_prod;
else
    A_prod = 1;
end
    
% set some value for R, fundamental for solution of unconstrained problem
R = params.pR;
W = params.W;
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
% Adjustment cost parameters
psi_k = params.psi_k; % capital
psi_l = params.psi_l; % labor
% Volatility
if isfield(params, 'volatility')
    volatility = params.volatility;
else
    volatility = 1;
end
ShockType = params.ShockType;


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



if isfield(params, 'Utilization')
    Utilization = params.Utilization;
else
    Utilization = 'Yes';
end


if isfield(params, 'LaborUtilization')
    LaborUtilization = params.LaborUtilization;
else
    LaborUtilization = 'No';
end

%% Functions

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


opts = optimset('display', 'off');

K_star = @(P) fsolve(@(k) solve_max_k(k,P), K_star_guess(P), opts);


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

 %% Boundaries Labor
if isfield(params, 'lmin') 
    lmin = params.lmin;
else
    lmin = 0;
end
 
if isfield(params, 'lmax')
    lmax = params.lmax;
else
    lmax = (1 - Gamma)/(Gamma) * k_bar(pmax);
end


if ~ ReducedOutput
    %% grids

    %option to have unevenly spaced grid Gamma > 1 gives more points for lower
    %K
    if isfield(params, 'GridKGamma')
        GridKGamma = params.GridKGamma;
    else
        GridKGamma = 1;
    end

    if isfield(params, 'GridLGamma')
        GridLGamma = params.GridLGamma;
    else
        GridLGamma = 1;
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
end



%% EXTRACT GRIDS and kmax

if ~ReducedOutput
    out.ll = ll;
    out.pp = pp;
    out.kk = kk;
    out.ind_up_k = ind_up_k;
    out.ind_dwn_k = ind_dwn_k;
    out.ind_up_p = ind_up_p;
    out.ind_dwn_p = ind_dwn_p;
    out.ind_up_l = ind_up_l;
    out.ind_dwn_l = ind_dwn_l;
    out.dk_vec_up = dk_vec_up;
    out.dl_vec_up = dl_vec_up;
    out.dp_vec_up = dp_vec_up;
    out.dk_vec_dwn = dk_vec_dwn;
    out.dl_vec_dwn = dl_vec_dwn;
    out.dp_vec_dwn = dp_vec_dwn;
    out.dp_mean = dp_mean;
    out.dk_mean = dk_mean;
    out.dl_mean = dl_mean;
    out.ind_p = ind_p;
    out.F_p = F_p;
    out.B_p = B_p;
    out.pmin = pmin;
    out.pmax = pmax;
    out.kmin = kmin;
    out.kmax = kmax;
    out.lmin = lmin;
    out.lmax = lmax;
    %% DEFINE INTEGRALS
    if strcmp(ShockType, 'Diffusion')
        out.integrate_kp = @(fun) sum(fun.*dk_mean.*dp_mean.*dl_mean);
    else
        out.integrate_kp = @(fun) sum(fun.*dk_mean.*dl_mean);
    end
else
    out.pmin = pmin;
    out.pmax = pmax;
    out.kmin = kmin;
    out.kmax = kmax;
    out.lmin = lmin;
    out.lmax = lmax;
end










