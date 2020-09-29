function out = SolvePolicy_rigidlabor(V, params) 

% Fornino Manera (2019)
%
% This subroutine is called by SolveTransition.m to find policies and 
% returns at time t given the value function at time t. This is used when 
% solving the transitional dynamics iterating the HJB backwards.
%
% Input
% - V: value function
% - params: parameters of main model
%
% Outputs("out" struct object, same as LaborDemand_trapz.m)
% - out: struct object containing a host of information on the transitional
%        dynamics


if isfield(params,'ReversibleInvestment')
    ReversibleInvestment = params.ReversibleInvestment;
else
    ReversibleInvestment = 'Yes';
end


%% Set  Parameter Values: 

if isfield(params, 'A_prod')
    A_prod = params.A_prod;
else
    A_prod = 1;
    out.params.A_prod = 1;
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

sep = params.sep;
% Volatility
if isfield(params, 'volatility')
    volatility = params.volatility;
else
    volatility = 1;
    out.params.volatility = volatility;
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

%% EXTRACT GRIDS
ind = (1:1:(params.N_k * params.N_l * params.N_p))';
ind_up_k = params.outGrids.ind_up_k;
ind_dwn_k  = params.outGrids.ind_dwn_k;
ind_up_p  = params.outGrids.ind_up_p;
ind_dwn_p  = params.outGrids.ind_dwn_p;
ind_up_l  = params.outGrids.ind_up_l;
ind_dwn_l  = params.outGrids.ind_dwn_l;
dk_vec_up  = params.outGrids.dk_vec_up;
dl_vec_up  = params.outGrids.dl_vec_up;
dk_vec_dwn = params.outGrids.dk_vec_dwn;
dl_vec_dwn = params.outGrids.dl_vec_dwn;
ind_p = params.outGrids.ind_p;
F_p = params.outGrids.F_p;
B_p = params.outGrids.B_p;
kk = params.outGrids.kk;
ll = params.outGrids.ll;
pp = params.outGrids.pp;
kmax = params.outGrids.kmax;
kmin = params.outGrids.kmin;
lmax = params.outGrids.lmax;
lmin = params.outGrids.lmin;
DIAG_p = - F_p - B_p;

%%

% set up the sparse rows and columns for the matrices
sparse_rows = [repmat(ind, 5 ,1); ind_p];
sparse_cols = [ind; ind_up_k;  ind_dwn_k;...
    ind_up_l;  ind_dwn_l; ...
     ind_up_p;  ind_dwn_p];

Pi_vec = Pi(kk, ll,  pp);

%% FORWARD DIFFERENCES
% Forward differences: get a zero when the step forward is not possible
dVF_k = (V(ind_up_k) - V(ind))./dk_vec_up;
dVB_k = - (V(ind_dwn_k) - V(ind))./dk_vec_dwn;

dVF_l = (V(ind_up_l) - V(ind))./dl_vec_up;
dVB_l = - (V(ind_dwn_l) - V(ind))./dl_vec_dwn;

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
%% GET DRIFTS AND POPULATE A
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

out.u = Pi_vec + ...
    ( - I_F_k .* Phi_F_k - I_B_k .* Phi_B_k - I_0_k .* Phi_0_k  ...
    - I_F_l .* Phi_F_l - I_B_l .* Phi_B_l - I_0_l .* Phi_0_l );


% Construct the sparse matrix
F_k = I_F_k .* s_F_k ./ dk_vec_up; % positive drifts that are actually positive
B_k = I_B_k .* s_B_k ./ dk_vec_dwn; 

F_l = I_F_l .* s_F_l ./ dl_vec_up; % positive drifts that are actually positive
B_l = I_B_l .* s_B_l ./ dl_vec_dwn; 

DIAG = DIAG_p - F_k + B_k - F_l + B_l ;
% Populate A

out.A = sparse(sparse_rows, sparse_cols, [DIAG ;...
    F_k; -B_k; F_l; -B_l; F_p; B_p]);


out.handles.Pi = Pi;
out.handles.R_bar = k_bar;
if strcmp(Utilization,'Yes')
    out.handles.R_hat = k_hat;
end



end
