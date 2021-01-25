% Fornino Manera (2019)
% 
% DATE: January 23, 2021
%
% Project: Automation and the Future of Work: Assessing the Role of Labor
%          Flexibility
%
%
% This routine takes a very long time to process. For convenience, we provide
% the file EOU_cs_workspace.mat which contains the pre-computed comparative
% statics.

%% Housekeeping & Setup

clc
clear
close all

% Current path
current_path = cd;

%% Run Calibration and Comparative Statics for Exponential Ornstein-Uhlenbeck Case

% Relative price of robots:
pR_rel10 = 1.4348;
pR_rel14 = 1.0209;

% Set up robustness on elasticity of penetration to robot prices
pR_rel_pre = pR_rel10;
pR_pct_chg = [0.5, 0.75, 1, 1.25, 1.5];
params = cell(length(pR_pct_chg), 1); 

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
    params{ii} = SetParametersGE(pR_rel_pre, pR_rel_post);
    disp(['\psi_R = ', num2str(params{ii}{1}.psi_k)])
    disp(' ')
    
    out{ii} = generalEquilibrium(params{ii});
    for jj = 1:13; prices(jj,ii) = out{ii}{jj}.params.A_prod; end
    wage(1, ii) = out{ii}{1}.params.W;
    
end

% Run Comparative Statics Exercise
out_mrts = cell(length(pR_pct_chg),1);
out_mrts_PE = out_mrts;

out_mrts_singlesector = out_mrts;
out_mrts_singlesector_PE = out_mrts;

out_adjCost = out_mrts;
out_adjCost_PE = out_mrts;

out_pR = out_mrts;
out_pR_PE = out_mrts;

for jj = 1:length(pR_pct_chg)
    
    disp(['Running case ' num2str(jj) ' of ' num2str(length(pR_pct_chg))])
    
    % Set Up Comparative Statics for BAR charts  
    mrts_ratio =  .75 * ones(13,1);
    adj_cost_ratio = .75;
    pR_rel_future = 0;
    mrts_ratio_singlesector = [.75; ones(12,1)];

    % Remember that if L_target = 1, then chi is unaffected by varphi.
    for ii = 1:13; params{jj}{ii}.varphi = 1; end

    % MRTS change
    params_mrts = params{jj};
    for sector = 1:13
        params_mrts{sector}.Gamma = mrts_ratio(sector) .* ...
            params{jj}{sector}.Gamma ./ (1 - params{jj}{sector}.Gamma + ...
            mrts_ratio(sector) .* params{jj}{sector}.Gamma);
    end
    out_mrts{jj} = generalEquilibrium(params_mrts);
    out_mrts_PE{jj} = partialEquilibrium([prices(:,jj); wage(jj)], params_mrts);
    clear params_mrts

    % MRTS change SINGLE SECTOR
    params_mrts_singlesector = params{jj};
    for sector = 1:13
        params_mrts_singlesector{sector}.Gamma = mrts_ratio_singlesector(sector) .* ...
            params{jj}{sector}.Gamma ./ (1 - params{jj}{sector}.Gamma + ...
            mrts_ratio_singlesector(sector) .* params{jj}{sector}.Gamma);
    end
    out_mrts_singlesector{jj} = generalEquilibrium(params_mrts_singlesector);
    out_mrts_singlesector_PE{jj} = partialEquilibrium([prices(:,jj); wage(jj)], params_mrts_singlesector);
    clear params_mrts_singlesector

    % Adjustment Cost change
    params_adjCost = params{jj};
    for sector = 1:13
        params_adjCost{sector}.psi_k = adj_cost_ratio * params{jj}{sector}.psi_k;
    end
    out_adjCost{jj} = generalEquilibrium(params_adjCost);
    out_adjCost_PE{jj} = partialEquilibrium([prices(:,jj); wage(jj)], params_adjCost);
    clear params_adjCost

    % Price of Robots change
    params_pR = params{jj};
    for sector = 1:13
        params_pR{sector}.pR = pR_rel_future * params{jj}{sector}.pR;
    end
    out_pR{jj} = generalEquilibrium(params_pR);
    out_pR_PE{jj} = partialEquilibrium([prices(:,jj); wage(jj)], params_pR);
    clear params_pR

end

% Run Continuous Comparative Statics Exercise

N_cs = 50;
mrts_ratio_vec = linspace(1, 0.05, N_cs);
mrts_ratio_singlesector_vec = linspace(1, 0.05, N_cs);
adj_cost_ratio_vec = linspace(1, 0.005, N_cs);
pR_rel_future_vec = linspace(1, 0, N_cs);

parameters_cs = [[mrts_ratio_vec'; ones(N_cs*3,1)], ...
                 [ones(N_cs,1); mrts_ratio_singlesector_vec'; ones(N_cs*2,1)], ...
                 [ones(N_cs*2,1); adj_cost_ratio_vec'; ones(N_cs,1)],...
                 [ones(N_cs*3,1); pR_rel_future_vec']];
parameters_cs = repmat(parameters_cs, length(pR_pct_chg), 1);
% out_cs = cell(N_cs * 4 * length(pR_pct_chg), 1);
% out_cs_PE = cell(N_cs * 4 * length(pR_pct_chg), 1);
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

% EMPLOYMENT   
which_base_idx = 3;
pct_chg_L = zeros(13,4);
for jj = 1:13
    pct_chg_L(jj,1) = (out_pR{which_base_idx}{jj}.L - out{which_base_idx}{jj}.L) ./ out{which_base_idx}{jj}.L;
    pct_chg_L(jj,2) = (out_mrts{which_base_idx}{jj}.L - out{which_base_idx}{jj}.L) ./ out{which_base_idx}{jj}.L;
    pct_chg_L(jj,3) = (out_mrts_singlesector{which_base_idx}{jj}.L - out{which_base_idx}{jj}.L) ./ out{which_base_idx}{jj}.L;
    pct_chg_L(jj,4) = (out_adjCost{which_base_idx}{jj}.L - out{which_base_idx}{jj}.L) ./ out{which_base_idx}{jj}.L;
end
plot_data = 100 .* pct_chg_L;

pct_chg_L_PE = zeros(13,4);
for jj = 1:13
    pct_chg_L_PE(jj,1) = (out_pR_PE{which_base_idx}{jj}.L - out{which_base_idx}{jj}.L) ./ out{which_base_idx}{jj}.L;
    pct_chg_L_PE(jj,2) = (out_mrts_PE{which_base_idx}{jj}.L - out{which_base_idx}{jj}.L) ./ out{which_base_idx}{jj}.L;
    pct_chg_L_PE(jj,3) = (out_mrts_singlesector_PE{which_base_idx}{jj}.L - out{which_base_idx}{jj}.L) ./ out{which_base_idx}{jj}.L;
    pct_chg_L_PE(jj,4) = (out_adjCost_PE{which_base_idx}{jj}.L - out{which_base_idx}{jj}.L) ./ out{which_base_idx}{jj}.L;
end
plot_data_PE = 100 .* pct_chg_L_PE;
plot_data_GE = plot_data - plot_data_PE;

clear out_pR out_mrts out_mrts_singlesector out_adjCost out_pR_PE out_mrts_PE out_mrts_singlesector_PE out_adjCost_PE

save EOU_cs_workspace.mat

