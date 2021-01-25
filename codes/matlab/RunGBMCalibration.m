% Fornino Manera (2019)
% 
% DATE: January 23, 2021
%
% Project: Automation and the Future of Work: Assessing the Role of Labor
%          Flexibility
%
%
% This routine takes a very long time to process. For convenience, we provide
% the file GBM_cs_workspace.mat which contains the pre-computed comparative
% statics.

%% Housekeeping & Setup

clc
clear
close all

% Current path
current_path = cd;

%% Run Calibration and Comparative Statics for Geometric Brownian Motion Case
% Relative price of robots:
pR_rel10 = 1.4348;
pR_rel14 = 1.0209;

% Set up robustness on elasticity of penetration to robot prices
pR_rel_pre = pR_rel14;
pR_pct_chg = 1;
params = cell(length(pR_pct_chg), 1); 
psi_fixed = [];
gammas_fixed = [];
out = params;
prices = zeros(13,length(pR_pct_chg));
wage = zeros(1, length(pR_pct_chg));

disp('REVISION: CALIBRATION OF GENERAL EQUILIBRIUM MODEL IN PROGRESS...')
disp(' ')
for ii = 1:length(pR_pct_chg)
    disp(['Price change in % of actual measured change: ', num2str(pR_pct_chg(ii) * 100), '%'])
    disp(' ')
    pct_chg = pR_pct_chg(ii) * (pR_rel14- pR_rel10)/pR_rel10 ;
    pR_rel_post =  pR_rel_pre*(1+pct_chg);
    params{ii} = SetParametersGE(pR_rel_pre, pR_rel_post, psi_fixed, gammas_fixed,'GBM');
    disp(['\psi_R = ', num2str(params{ii}{1}.psi_k)])
    disp(' ')
    
    out{ii} = generalEquilibrium(params{ii});
    for jj = 1:13; prices(jj,ii) = out{ii}{jj}.params.A_prod; end
    wage(1, ii) = out{ii}{1}.params.W;
    
end



% Run Continuous Comparative Statics Exercise
N_cs = 50;
mrts_ratio_vec = linspace(1, 0.05, N_cs);
mrts_ratio_singlesector_vec = linspace(1, 0.05, N_cs);
adj_cost_ratio_vec = linspace(1, 0.01, N_cs);
pR_rel_future_vec = linspace(1, 0, N_cs);

parameters_cs = [[mrts_ratio_vec'; ones(N_cs*3,1)], ...
                 [ones(N_cs,1); mrts_ratio_singlesector_vec'; ones(N_cs*2,1)], ...
                 [ones(N_cs*2,1); adj_cost_ratio_vec'; ones(N_cs,1)],...
                 [ones(N_cs*3,1); pR_rel_future_vec']];
parameters_cs = repmat(parameters_cs, length(pR_pct_chg), 1);

L_cs = zeros(N_cs * 4 * length(pR_pct_chg), 1);
L_cs_PE = L_cs;

LS_cs = L_cs;
LS_cs_PE = L_cs;

RL_cs = L_cs;

disp('SMOOTH COMPARATIVE STATICS ON PARAMETERS OF INTEREST...')
parfor ii = 1:N_cs*4 * length(pR_pct_chg)
    kk = ceil(ii / (N_cs*4));
    params_cs = params{kk};
    params_base = params{kk};
    
    for jj = 1:13
        params_cs{jj}.ReducedOutput = true;
        if ceil(rem(ii, N_cs*4) / N_cs) == 2
            params_cs{1}.Gamma = parameters_cs(ii,2) .* params_base{1}.Gamma ./ ...
                (1 - params_base{1}.Gamma + parameters_cs(ii,2) .* params_base{1}.Gamma);
        else
            params_cs{jj}.Gamma = parameters_cs(ii,1) .* params_base{jj}.Gamma ./ ...
                (1 - params_base{jj}.Gamma + parameters_cs(ii,1) .* params_base{jj}.Gamma);
        end
        params_cs{jj}.psi_k = parameters_cs(ii,3) * params_base{jj}.psi_k;
        params_cs{jj}.pR = parameters_cs(ii,4) * params_base{jj}.pR;
    end
    
    out_cs = generalEquilibrium(params_cs);
    out_cs_PE = partialEquilibrium([prices(:, kk); wage(kk)], params_cs);
    
    nY = 0;
    nY_PE = 0;
    R = 0;
    for jj = 1:13
        L_cs(ii) = L_cs(ii) + out_cs{jj}.L;
        L_cs_PE(ii) = L_cs_PE(ii) + out_cs_PE{jj}.L;
        nY = nY + out_cs{jj}.params.A_prod * out_cs{jj}.GE.q;
        nY_PE = nY_PE + prices(jj, kk) * out_cs_PE{jj}.GE.q;
        R = R + out_cs_PE{jj}.RL * out_cs_PE{jj}.L;
    end
    RL_cs(ii) = R/L_cs(ii);
    LS_cs(ii) = out_cs{1}.W * L_cs(ii) / nY;
    LS_cs_PE(ii) = wage(kk) * L_cs_PE(ii) / nY_PE;

    disp(['Completed Iteration: ' num2str(ii) '.'])
    
end

clear out

save GBM_cs_workspace.mat
