% Fornino Manera (2019)
% 
% DATE: January 24, 2021
%
% Project: Automation and the Future of Work: Assessing the Role of Labor
%          Flexibility

function out_transition = SolveTransition(transition, params)
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


% Current path
current_path = cd;


% par_model = transition.model;
par_path = transition.par_path;
par_str = transition.par_str;
D = transition.D;
NT = length(par_path);
Ns = params.N_k * params.N_p;

print_conv = params.print_conv;

% STEP 1: compute initial & final steady states

params_initial = params;
params_initial.(par_str) = par_path(1);
params_final = params;
params_final.(par_str) = par_path(end);

% par_path(1) = par_path(2);
par_path = [par_path(2:end); par_path(end)]; 

% disp('Solving for starting point...')
calibrate_chi = true;
out_initial = generalEquilibrium_onesector(params_initial, calibrate_chi);
% L_initial = out_initial.L;
chi = out_initial.params.chi;
params_initial.chi = chi;
params.chi = chi;
params_final.chi = chi;
% disp('Done!')

% disp('Solving for final point...')
calibrate_chi = false;
out_final = generalEquilibrium_onesector(params_final, calibrate_chi);
% L_final = out_final.L;
% disp('Done!')

% disp(['Labor at initial point = ' num2str(L_initial)])
% disp(['Labor at final point = ' num2str(L_final)])

params.kmax = max(out_initial.params.k_max, out_final.params.k_max);
params.kmin = 0;
params_initial.kmax = params.kmax;
params_initial.kmin = params.kmin;
params_final.kmax = params.kmax;
params_final.kmin = params.kmin;

% disp('Solving for starting point...')
out_initial = generalEquilibrium_onesector(params_initial, calibrate_chi);
% L_initial = out_initial.L;
% chi = out_initial.params.chi;
% params_initial.chi = chi;
% params.chi = chi;
% disp('Done!')

% disp('Solving for final point...')
% calibrate_chi = false;
% params_final.chi = chi;
out_final = generalEquilibrium_onesector(params_final, calibrate_chi);
% L_final = out_final.L;
% disp('Done!')

outGrids = GetGrids(params, false);

% Guess for path of wage
try
    W_path = readmatrix([current_path '/guesses/wpath.csv']);
catch
    W_path = linspace(out_initial.W, out_final.W, NT)';
end


% Labor Supply
LS = @(w) (w./chi) .^ params.varphi;
% LS_path = LS(W_path);

factor = .05;
maxit_w = 500;
it_w = 0;
supnorm = 1;
tol = 1e-4;
while it_w < maxit_w && supnorm > tol

    W_path = W_path(1:NT);
    par_path = par_path(1:NT);
    
    % STEP 2: solve HJB backwards and KFE forwards.
    % Preallocate necessary objects
    params_trans = params;
    params_trans.outGrids = outGrids;

    g_seq = cell(NT, 1);
    v_seq = cell(NT, 1);
    LD_path = zeros(NT,1);
    LD_path(end) = out_final.L;
    g_seq{1} = out_initial.statDist;
    v_seq{end} = out_final.V;
    A_next = out_final.A;
    u_next = out_final.u;


    % SOLVE HJB BACKWARDS
%     disp('Solving HJB Equation backward...')
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
        params_trans.(par_str) = par_path(jt) * W_path(jt);
%         params_trans.(par_str) = par_path(jt);
        params_trans.W = W_path(jt);
        
        % Solve for policy given value function
        out_it = SolvePolicy_trapz(v_seq{jt}, params_trans);

        % Update sequence of Infinitesimal Generators encoding policies
        A_next = out_it.A;
        u_next = out_it.u;
%         handles_seq{jt} = out_it.handles;

        if rem(jt + 1, rem_factor) == 0 && print_conv
            t1 = toc;
            disp(['Iteration: ' num2str(NT-jt) ' of ' num2str(NT) '. ETA = ' num2hms( (t1 - t0)/(NT - jt) * (jt - 1) ) '.'])
        end

    end
%     disp('Done!')

    % Solve KFE
%     disp('Solving KFE Equation forward...')
    t0 = toc;
    for jt = 2:1:NT

        % Set right prices    
        params_trans.(par_str) = par_path(jt - 1) * W_path(jt - 1);
%         params_trans.(par_str) = par_path(jt - 1);
        params_trans.W = W_path(jt - 1);
        
        % Re-compute policies to save on memory
        out_it = SolvePolicy_trapz(v_seq{jt - 1}, params_trans);

        % Update endogenous distribution
        g_seq{jt} = (speye(Ns) - D * out_it.A') \ g_seq{jt - 1};
        
        % Compute aggregate labor
        LD_path(jt - 1) = out_it.params.outGrids.integrate_kp(g_seq{jt - 1} .* out_it.indL);
        
        if rem(jt-2,rem_factor) == 0 && print_conv
            t1 = toc;
            disp(['Iteration: ' num2str(jt - 1) ' of ' num2str(NT) '. ETA = ' num2hms( (t1 - t0) * (NT - jt - 1) / jt) '.'])
        end
    end
%     disp('Done!')

%     MC = LD_path(1:end-1) - LS(W_path(1:end-1));
%     W_path(1:end-1) = W_path(1:end-1) + factor .* MC ;
%     
    params_trans.(par_str) = par_path(NT) * W_path(NT);
    params_trans.W = W_path(NT);
    out_it = SolvePolicy_trapz(v_seq{NT}, params_trans);
    LD_path(end) = outGrids.integrate_kp(g_seq{NT} .* out_it.indL);
%     LD_path(end) = g_seq{NT}' * out_it.indL;
    
    MC = LD_path - LS(W_path);
    W_path = W_path + factor .* MC ;
    
%     subplot(3,1,1)
%     hold on
%     plot(W_path(1:end-1), 'linewidth', 2);
%     ylim([min(W_path(1:end-1)), max(W_path(1:end-1))])
%     hold off
%     
%     subplot(3,1,2)
%     hold on
%     plot(LD_path(1:end-1), 'linewidth', 2);
%     ylim([min(LD_path(1:end-1)), max(LD_path(1:end-1))])
%     hold off
%     
%     subplot(3,1,3)
%     hold on
%     plot(MC, 'linewidth', 2);
%     ylim([min(MC), max(MC)])
%     hold off
%     drawnow()
    
    

    NT_idx = 1 + find((abs(W_path(2:end) - out_final.W) ./ out_final.W < tol/10) ...
                  & (abs(diff(W_path) ./ W_path(1:NT-1)) < tol/10) & ...
                  (abs(diff(par_path)) ./ par_path(1:NT-1) < tol/10), 1);
    
    if ~isempty(NT_idx)
        NT = NT_idx;
    end
    
%     disp(num2str(NT))
    
    supnorm = max(abs(MC));
    it_w = it_w + 1;
    
    disp(['Iteration: ' num2str(it_w) '. Supnorm = ' num2str(supnorm) '. Convergence at T = ' num2str(NT * D)])
end

writematrix(W_path, [current_path '/guesses/wpath.csv'])


% STEP 3: produce output
% compute robots per thousand employees
RL_path = zeros(NT,1);
for jt = 1:NT
    RL_path(jt) = outGrids.integrate_kp( g_seq{jt} .* outGrids.kk ) ./ ...
        LD_path(jt) * 1e3;
%     RL_path(jt) = (g_seq{jt}' * outGrids.kk) ./ LD_path(jt) * 1e3;
end

out_transition.out_initial = out_initial;
out_transition.out_final = out_final;

out_transition.outGrids = outGrids;
out_transition.g_seq = g_seq;
out_transition.v_seq = v_seq;
out_transition.LD_path = LD_path;
out_transition.RL_path = RL_path;
out_transition.RL_initial = out_initial.RL;

T = NT * D;

t_bag = D:D:T;

t_sucre = 1:find(t_bag == transition.T_note, 1);

n_abbruzzi = 4;

close

% subplot(n_abbruzzi,1,1)
% title('R/L')
% hold on
% plot(t_bag(t_sucre), RL_path(t_sucre), 'linewidth', 2)
% line([t_bag(1) t_bag(t_sucre(end))], [out_initial.RL* 1e3, out_initial.RL* 1e3], 'linestyle', '--', 'linewidth', 2, 'color', 'k')
% line([t_bag(1) t_bag(t_sucre(end))], [out_final.RL* 1e3, out_final.RL* 1e3], 'linestyle', '--', 'linewidth', 2, 'color', 'k')
% hold off
% 
% subplot(n_abbruzzi,1,2)
% title('wage')
% hold on
% plot(t_bag(t_sucre), W_path(t_sucre), 'linewidth', 2)
% line([t_bag(1) t_bag(t_sucre(end))], [out_initial.W, out_initial.W], 'linestyle', '--', 'linewidth', 2, 'color', 'k')
% line([t_bag(1) t_bag(t_sucre(end))], [out_final.W, out_final.W], 'linestyle', '--', 'linewidth', 2, 'color', 'k')
% hold off
% 
% subplot(n_abbruzzi,1,3)
% title('path of p_R / w')
% hold on
% plot(t_bag(t_sucre), par_path(t_sucre), 'linewidth', 2)
% % line([t_bag(1) t_bag(end)], [out_initial.W, out_initial.W], 'linestyle', '--', 'linewidth', 2, 'color', 'k')
% % line([t_bag(1) t_bag(end)], [out_final.W, out_final.W], 'linestyle', '--', 'linewidth', 2, 'color', 'k')
% hold off
% 
% subplot(n_abbruzzi,1,4)
% title('Labor')
% hold on 
% plot(t_bag(t_sucre), LD_path(t_sucre), 'linewidth', 2)
% hold off
% % 
% % subplot(4,1,4)
% % title('MC')
% % hold on
% % plot(t_bag, MC , 'linewidth', 2)
% % hold off
% 
% drawnow()

end




