% Fornino Manera (2019)
% 
% DATE: January 24, 2021
%
% Project: Automation and the Future of Work: Assessing the Role of Labor
%          Flexibility

function out_transition = SolveTransition_rigidlabor(transition, params)
% This subroutine is called by RunFigures.m. It solves for the transitional
% dynamics of the model with rigid labor. 
%
% Input
% - transition: struct object containing the specifications of the
%               transition of interest (e.g., path of exogenous variable) 
% - params: parameters of main model
%
% Outputs("out" struct object, same as LaborDemand_trapz.m)
% - out_transition: struct object containing a host of information on the
%                   transitional dynamics

% Check arguments

% par_model = transition.model;
par_path = transition.par_path;
par_str = transition.par_str;
D = transition.D;
NT = length(par_path);
Ns = params.N_k * params.N_l * params.N_p;
% T = NT * D;
print_conv = params.print_conv;

% STEP 1: compute initial & final steady states

params_initial = params;
params_initial.(par_str) = par_path(1);
params_final = params;
params_final.(par_str) = par_path(end);

% find right kmax and lmax
outGrids_initial = GetGrids_rigidlabor(params_initial, true);
outGrids_final = GetGrids_rigidlabor(params_final, true);

params_initial.lmin = min(outGrids_initial.lmin, outGrids_final.lmin);
params_initial.kmin = min(outGrids_initial.kmin, outGrids_final.kmin);
params_initial.pmin = min(outGrids_initial.pmin, outGrids_final.pmin);
params_initial.lmax = max(outGrids_initial.lmax, outGrids_final.lmax);
params_initial.kmax = max(outGrids_initial.kmax, outGrids_final.kmax);
params_initial.pmax = max(outGrids_initial.pmax, outGrids_final.pmax);

params_final.lmin = min(outGrids_initial.lmin, outGrids_final.lmin);
params_final.kmin = min(outGrids_initial.kmin, outGrids_final.kmin);
params_final.pmin = min(outGrids_initial.pmin, outGrids_final.pmin);
params_final.lmax = max(outGrids_initial.lmax, outGrids_final.lmax);
params_final.kmax = max(outGrids_initial.kmax, outGrids_final.kmax);
params_final.pmax = max(outGrids_initial.pmax, outGrids_final.pmax);

params.lmin = min(outGrids_initial.lmin, outGrids_final.lmin);
params.kmin = min(outGrids_initial.kmin, outGrids_final.kmin);
params.pmin = min(outGrids_initial.pmin, outGrids_final.pmin);
params.lmax = max(outGrids_initial.lmax, outGrids_final.lmax);
params.kmax = max(outGrids_initial.kmax, outGrids_final.kmax);
params.pmax = max(outGrids_initial.pmax, outGrids_final.pmax);

outGrids = GetGrids_rigidlabor(params, false);

disp('Solving for starting point...')
[~, out_initial] = LaborDemand_rigidlabor(params_initial.W, params_initial.pR, params_initial);
disp('Done!')

disp('Solving for final point...')
[~, out_final] = LaborDemand_rigidlabor(params_final.W, params_final.pR, params_final);
disp('Done!')

disp(['Labor at initial point = ' num2str(out_initial.L_contracts)])
disp(['Labor at final point = ' num2str(out_final.L_contracts)])

% STEP 2: solve HJB backwards and KFE forwards.

% Preallocate necessary objects
params_trans = params;
params_trans.outGrids = outGrids;

g_seq = cell(NT, 1);
v_seq = cell(NT, 1);
handles_seq = cell(NT, 1);
g_seq{1} = out_initial.statDist_pmf;
v_seq{end} = out_final.V;
A_next = out_final.A;
u_next = out_final.u;


% SOLVE HJB BACKWARDS
disp('Solving HJB Equation backward...')
tic
t0 = toc;
rem_factor = floor(NT / 20);
for jt = NT-1:-1:1
    
    % Solve for value function at t as function of V_{t + Delta} using
    % HJB (A_{t + 1}, u_{t + 1} encode policies and values)
    B_it = (params_trans.rho + 1/D) * speye(Ns) - A_next;
    b_it = u_next + 1/D * v_seq{jt + 1};
    v_seq{jt} = B_it \ b_it;

    % Set right prices    
    params_trans.(par_str) = par_path(jt);
    
    % Solve for policy given value function
    out_it = SolvePolicy_rigidlabor(v_seq{jt}, params_trans);

    % Update sequence of Infinitesimal Generators encoding policies
    A_next = out_it.A;
    u_next = out_it.u;
    handles_seq{jt} = out_it.handles;
    
    if rem(jt + 1, rem_factor) == 0 && print_conv
        t1 = toc;
        disp(['Iteration: ' num2str(NT-jt) ' of ' num2str(NT) '. ETA = ' num2hms( (t1 - t0)/(NT - jt) * (jt - 1) ) '.'])
    end
    
end
disp('Done!')

% Solve KFE
disp('Solving KFE Equation forward...')
t0 = toc;
for jt = 2:1:NT
    
    % Set right prices    
    params_trans.(par_str) = par_path(jt - 1);

    % Re-compute policies to save on memory
    out_it = SolvePolicy_rigidlabor(v_seq{jt - 1}, params_trans);
    
    % Update endogenous distribution
    g_seq{jt} = (speye(Ns) - D * out_it.A') \ g_seq{jt - 1};
    
%     g_seq{jt} = (D * out_it.A' + speye(Ns)) * g_seq{jt-1}; 

    if rem(jt-2,rem_factor) == 0 && print_conv
        t1 = toc;
        disp(['Iteration: ' num2str(jt - 1) ' of ' num2str(NT) '. ETA = ' num2hms( (t1 - t0) * (NT - jt - 1) / jt) '.'])
    end
end
disp('Done!')

% STEP 3: produce output

% compute path of labor
L_path = zeros(NT,1);
for jt = 1:NT
    L_path(jt) = outGrids.integrate_kp( g_seq{jt}./...
        (outGrids.dp_mean.*outGrids.dk_mean.*outGrids.dl_mean) .* outGrids.ll );
end

% compute robots per thousand employees
RL_path = zeros(NT,1);
for jt = 1:NT
    RL_path(jt) = outGrids.integrate_kp( g_seq{jt}./...
        (outGrids.dp_mean.*outGrids.dk_mean.*outGrids.dl_mean) .* outGrids.kk ) ./ ...
        L_path(jt) * 1e3;
end

out_transition.out_initial = out_initial;
out_transition.out_final = out_final;

out_transition.outGrids = outGrids;
out_transition.g_seq = g_seq;
out_transition.v_seq = v_seq;
out_transition.L_path = L_path;
out_transition.RL_path = RL_path;

end




