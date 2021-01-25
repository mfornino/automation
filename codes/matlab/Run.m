% Fornino Manera (2019)
% 
% DATE: January 24, 2021
%
% Project: Automation and the Future of Work: Assessing the Role of Labor
%          Flexibility
%
%
% MAIN MATLAB ROUTINE TO DRAW FIGURES AND COMPUTE TABLES
% Please see the attached documentation.

%% Housekeeping & Setup

clc
clear
close all
set(groot,'defaultAxesTickLabelInterpreter','latex');

% Toggle to save the graphs to folder (0 don't save, 1 save)
save_figs = true;

% Save Path (relative to current working directory - must include final /)
save_path_figures = '/graphs/';
save_path_tables = '/tables/';
current_path = cd;

% Toggle to select figure format (accepts 'eps', 'jpeg', 'png')
fig_format = 'eps';

% Toggle to select docked or windowed mode
% Note: set to false when saving figures to file.
docked = false;

% Toggle to select persistent mode (keeps figures open for consultation)
persistent_mode = true;

% print long numbers
format long

if docked
    set(0,'DefaultFigureWindowStyle','docked');
else
    set(0,'DefaultFigureWindowStyle','normal');
end

%% 1) ***** Robot Price Data Figure *****
% These data points are from the report by Sirkin (2015a). For more information
% see paper.

Year = categorical({'2005',...
                    '2010',...
                    '2014',...
                    '2020 (*)',...
                    '2025 (*)'});
data = [55 33 81 13; 
        43 40 62 11;
        33 45 46 9;
        30 40 39 8;
        28 36 33 7];
      
% FIGURE SETUP
fig_name = '1_pR_adjcost_data';  % FIGURE NAME
fig_width = 650;        % Window Horizontal Size
fig_heigth = 400;       % Window Vertical Size
fig_posx = 160;         % Window position (Lower left corner)
fig_posy = 160;         % Window position (Lower left corner)

if docked
    figure()
else
    figure('Position', [fig_posx fig_posy fig_width fig_heigth])
end

b = bar(Year, data, 'stacked', 'EdgeColor', 'none');
b(1).FaceColor = [0.1 0.1 0.1 ];
b(2).FaceColor = [0 0.4470 0.7410]*.75;
b(3).FaceColor = [0 0.4470 0.7410]*1;
b(4).FaceColor = [0 0.4470 0.7410]*1.25;
ylabel('Nominal U.S. Dollars (thousands)',...
       'Interpreter', 'Latex')
legend({'Robot (including software)',...
        'Peripherals (such as safety barriers, sensors)',...
        'Systems Engineering (such as programming, installation)',...
        'Project Management'},...
        'Location', 'SouthOutside',...
        'Interpreter', 'Latex',...
        'Fontsize', 16,...
        'NumColumns', 2,...
        'Box', 'off')
box off

if save_figs
    cd([current_path save_path_figures])
    fname = strcat(fig_name);
    set(gcf,'PaperPositionMode','auto');
    hgsave(fname);
    if strcmp(fig_format, 'eps')
        print('-depsc',fname,'-painters');
    elseif strcmp(fig_format, 'jpeg')
        print('-djpeg',fname,'-r256','-opengl');
    elseif strcmp(fig_format, 'png')
        print('-dpng',fname);
    else
        error('Format specified incorrectly.')
    end
    if ~persistent_mode
        close
    end
    cd(current_path)
end

%% CLEAR (SOME) VARIABLES IN MEMORY AND CONTINUE %%
clearvars -EXCEPT save_figs fig_format save_path_tables save_path_figures current_path docked persistent_mode 

%% 2) ***** Run Illustrative Calibration *****

params = SetParameters();

params.W = 0.268544000136280;          % see single-sector calibration
params.pR = params.W * 1.02;           % wage times relative robot price in 2014
params.theta_P = .8789605;             % see HamiltonOUEstimatorsIFR_total.dta
params.sigma_P = .1418431;             % see HamiltonOUEstimatorsIFR_total.dta
params.t = .2978469;                   % see ThetaProdIFR_all.dta
params.E = 0;
params.Gamma = 0.652862965992859;      % see single-sector calibration
params.psi_k = 5;                      % ILLUSTRATIVE CHOICE

% Wage and Price of Robots
W = params.W;
pR = params.pR;

% Make Pretty Graphs
params.pmin = 1e-4;

% Compute Labor Savings and R_max
Omega = ( 1 - params.Gamma) /params.Gamma * params.W - params.E;

R_max = 1 / params.delta / params.psi_k * (Omega / (params.rho + params.delta) - pR);
params.kmax = 1.11 .* R_max;

[~, out] = LaborDemand_trapz(W, pR, params);
 
% PLOT OF STATIONARY DISTRIBUTION

% FIGURE SETUP
fig_name = '2_stationary_distribution';  % FIGURE NAME
fig_width = 650;        % Window Horizontal Size
fig_heigth = 400;       % Window Vertical Size
fig_posx = 160;         % Window position (Lower left corner)
fig_posy = 160;         % Window position (Lower left corner)

if docked
    figure()
else
    figure('Position', [fig_posx fig_posy fig_width fig_heigth])
end

subplot(5,5,[2 3 4 5 7 8 9 10 12 13 14 15 17 18 19 20])
alpha_value = .25;
pol_plot = reshape(out.statDist, [params.N_k params.N_p]);
p_check = ((1 - params.Gamma) * R_max)^(1 - params.t) * W / (params.A_prod * params.t * params.Gamma);
p_grid_plot = out.p_grid;
idxk = find(out.k_grid >= R_max,1);
pol_plot = pol_plot(1:idxk,:);
max_pol = max(max(pol_plot));
map = [linspace(.25, 0, 10)', linspace(.25, 0, 10)', linspace(.25, 0, 10)'];
contour(p_grid_plot,out.k_grid(1:idxk),pol_plot,20,...
        'linewidth',1);
colormap('parula')
% colormap(map)
% colormap('winter')
alpha(alpha_value)
ylim([0 R_max*1.1])
xlim([0 max(p_grid_plot)])
hold on
plot(p_grid_plot, out.INV_CUTOFFS_MEAN,...
    'linestyle','-',...
    'linewidth', 3,...
    'color','r')
plot(p_grid_plot, out.handles.R_bar(p_grid_plot),...
    'linestyle','-',...
    'linewidth', 3,...
    'color','k')
line([0, max(p_grid_plot)], [R_max, R_max],...
     'linestyle', '-.',...
     'linewidth', 2,...
     'color','k')
line([p_check, p_check], [0, R_max],...
     'linestyle', ':',...
     'linewidth', 2,...
     'color','k')
hold off
yticks([0 R_max])
yticklabels({'',''})
a = get(gca,'YTickLabel');
set(gca,'YTickLabel',a,'fontsize',14)
xticks([0, p_check])
xticklabels({'','$\check{z}$'})
legend({'Stationary Distribution','$R^{\star}(z)$','$\bar{R}(z)$'},...
        'interpreter','latex',...
        'fontsize',14,...
        'location','southeast',...
        'box', 'off')
view(0,+90)

subplot(5,5,[ 22 23 24 25])
colors_graph = get(gca, 'ColorOrder');
max_marg_pdf = max(max(lognpdf(linspace(0,p_check,1000), params.mu_hat, params.sigma_hat)),...
    max(lognpdf(linspace(p_check,max(p_grid_plot),1000), params.mu_hat, params.sigma_hat)));
hold on
area(linspace(p_check,max(p_grid_plot),1000),...
    -lognpdf(linspace(p_check,max(p_grid_plot),1000), params.mu_hat, params.sigma_hat) ,...
    'linestyle','none',...
    'FaceColor',[.5 .5 .5])
plot([linspace(0,p_check,1000) linspace(p_check,max(p_grid_plot),1000)],...
     -lognpdf([linspace(0,p_check,1000) linspace(p_check,max(p_grid_plot),1000)], params.mu_hat, params.sigma_hat), ...
     'linewidth', 2, 'color', 'k')
line([p_check, p_check], [-1.1*max_marg_pdf, 0],...
     'linestyle', ':',...
     'linewidth', 2,...
     'color','k')
hold off
alpha(.35)
ylim([-1.1*max_marg_pdf 0])
xlim([0 max(p_grid_plot)])
xlabel('Revenue-Shifter $z$','interpreter','latex', 'fontsize', 14)
ylabel('$f(z)$','interpreter','latex', 'fontsize', 14)
xticks([p_check])
xticklabels({'$\check{z}$'})
a = get(gca,'XTickLabel');
set(gca,'XTickLabel',a,'fontsize',14)
yticks([0])
yticklabels({''})
legend({'Partial-Automation Lower Bound'},...
        'interpreter','latex',...
        'fontsize',14,...
        'location','southeast',...
        'box', 'off')

    
subplot(5,5,[1 6 11 16])
pol_plot = reshape(out.statDist, [params.N_k params.N_p]);
dist_R = sum(pol_plot, 2);
colors_graph = get(gca, 'ColorOrder');
max_marg_pdf = max(max(lognpdf(linspace(0,p_check,1000), params.mu_hat, params.sigma_hat)),...
    max(lognpdf(linspace(p_check,max(p_grid_plot),1000), params.mu_hat, params.sigma_hat)));
plot(1.1*max(dist_R(1:find(out.k_grid > 1.1*R_max,1))) - dist_R(1:find(out.k_grid > 1.1*R_max,1)), out.k_grid(1:find(out.k_grid > 1.1*R_max,1)), ...
     'linewidth', 2, 'color', 'k')
xlim([0, 1.1*max(dist_R(1:find(out.k_grid > 1.1*R_max,1)))])
ylim([0 1.1*R_max])
xlabel('$\int_{0}^{\infty}g(R,z)\textrm{d}z$','interpreter','latex', 'fontsize', 14)
ylabel('Robot Stock $R$','interpreter','latex', 'fontsize', 14)
xticks([ 1.1*max(dist_R(1:find(out.k_grid > 1.1*R_max,1))) ])
xticklabels({''})
a = get(gca,'XTickLabel');
set(gca,'XTickLabel',a,'fontsize',14)
yticks([R_max])
yticklabels({'$R^{\star}_{\max}$'})
line([0, 1.1*max(dist_R(1:find(out.k_grid > 1.1*R_max,1)))], [R_max, R_max],...
     'linestyle', '-.',...
     'linewidth', 2,...
     'color','k')

if save_figs
    cd([current_path save_path_figures])
    fname = strcat(fig_name);
    set(gcf,'PaperPositionMode','auto');
    hgsave(fname);
    if strcmp(fig_format, 'eps')
        print('-depsc',fname,'-painters');
    elseif strcmp(fig_format, 'jpeg')
        print('-djpeg',fname,'-r256','-opengl');
    elseif strcmp(fig_format, 'png')
        print('-dpng',fname);
    else
        error('Format specified incorrectly.')
    end
    if ~persistent_mode
        close
    end
    cd(current_path)
end
 
%% 3) ***** Compute Simulated Paths for Illustrative Calibration *****

% Set seed
rng(1);

% Extract grids and policy functions
R_grid_figure = out.k_grid;
p_grid_figure = out.p_grid;
i_k_pol = reshape(out.i_k_pol - params.delta * out.kk, params.N_k, params.N_p);

% Run simulation
A = out.A;
R0_idx = find(R_grid_figure>0,1);
p0_idx = find(p_grid_figure>1,1);
s0_idx = sub2ind_vec([params.N_k params.N_p],[R0_idx, p0_idx]);
sim_steps = 10e3;
dt = zeros(sim_steps,1);
state = dt;
t = dt;
state(1) = s0_idx;
for step = 2:sim_steps
    [time_to_jump, state_to_jump] = min(exprnd(full(A(state(step-1),:)).^-1));
    dt(step) = time_to_jump;
    state(step) = state_to_jump;
    t(step) = t(step-1) + dt(step);
end
state_sub = ind2sub_vec([params.N_k,params.N_p],state);

R_sample_path = R_grid_figure(state_sub(:,1));
p_sample_path = p_grid_figure(state_sub(:,2));
Ld_sample_path = (out.handles.R_bar(p_sample_path)'- ...
     R_sample_path >= 0 ).* (out.handles.R_bar(p_sample_path)'- ...
     R_sample_path);
I_sample_path = i_k_pol(state);

% PLOT OF SIMULATED TIME SERIES

% Simulated Paths of Robots, Labor, Demand Shifter
% FIGURE SETUP
fig_name = '3_sample_paths';  % FIGURE NAME
fig_width = 525;	% Window Horizontal Size
fig_heigth = 350;   % Window Vertical Size
fig_posx = 160;     % Window position (Lower left corner)
fig_posy = 160;     % Window position (Lower left corner)

if docked
    figure()
else
    figure('Position', [fig_posx fig_posy fig_width fig_heigth])
end

subplot(4,1,1)
plot(t, p_sample_path,...
     'color','k',...
     'linewidth',1.5)
line([0,max(t)],[1 1],...
     'linestyle',':',...
     'color','k',...
     'linewidth',1)
title('Revenue-Shifter $z(t)$',...
      'interpreter','latex',...
      'fontsize',14)
xticks([0, max(t)])
xticklabels({'0','T'})
yticks([0, 1])
yticklabels({'0', '$1$'})
xlim([0, max(t)])
ylim([params.pmin params.pmax])

subplot(4,1,2)
plot(t, R_sample_path,...
     'color','k',...
     'linewidth',1.5)
title('Robots $R(t)$',...
      'interpreter','latex',...
      'fontsize',14)
xticks([0, max(t)])
xticklabels({'0','T'})
yticks([0])
xlim([0, max(t)])
ylim([0 1.1*R_max])

subplot(4,1,3)
plot(t, Ld_sample_path,...
     'color','k',...
     'linewidth',1.5)
title('Labor $L(t)$',...
      'interpreter','latex',...
      'fontsize',14)
xticks([0, max(t)])
xticklabels({'0','T'})
yticks([0])
xlim([0, max(t)])
ylim([0 1.1*max(Ld_sample_path)])

subplot(4,1,4)
plot(t, I_sample_path,...
     'color','k',...
     'linewidth',1.5)
line([0,max(t)],[0 0],...
     'linestyle',':',...
     'color','k',...
     'linewidth',1)
title('Net Investment $I(R(t),z(t))-\delta R(t)$',...
      'interpreter','latex',...
      'fontsize',14)
xticks([0, max(t)])
xticklabels({'0','T'})
yticks([0])
xlim([0, max(t)])
ylim([min(I_sample_path) max(I_sample_path)])
xlabel('Time $t$','interpreter','latex')

if save_figs
    cd([current_path save_path_figures])
    fname = strcat(fig_name);
    set(gcf,'PaperPositionMode','auto');
    hgsave(fname);
    if strcmp(fig_format, 'eps')
        print('-depsc',fname,'-painters');
    elseif strcmp(fig_format, 'jpeg')
        print('-djpeg',fname,'-r256','-opengl');
    elseif strcmp(fig_format, 'png')
        print('-dpng',fname);
    else
        error('Format specified incorrectly.')
    end
    if ~persistent_mode
        close
    end
    cd(current_path)
end

%% CLEAR (SOME) VARIABLES IN MEMORY AND CONTINUE %%
clearvars -EXCEPT save_figs fig_format save_path_tables save_path_figures current_path docked persistent_mode 

%% 4) ***** Calibrated Model in General Equilibrium *****

% For convenience, we provide the file EOU_cs_workspace.mat which contains the 
% pre-computed comparative statics. The reader can use RunEOUCalibration.m to
% produce these simulations.

% RunEOUCalibration

% LOAD WORKSPACE
load EOU_cs_workspace.mat

% SHOW CALIBRATION AND BUILD CALIBRATION TABLE
% Note that we stored the baseline scenario (without robustness on robot price
% data) in index 3 of the cell array.
which_base_idx = 3;

ifrStrings = cell(13,1);
for jj = 1:13
    ifrStrings{jj} = params{which_base_idx}{jj}.targets.ifrString;
end

% Set up the table and assign row names
calibrationtable = table(ones(13,1));
calibrationtable.Properties.RowNames = ifrStrings;
calibrationtable.Var1 = [];

% Include targets
for jj = 1:13
    calibrationtable.Robots10(jj) = params{which_base_idx}{jj}.targets.rl_target10;
    calibrationtable.Robots14(jj) = params{which_base_idx}{jj}.targets.rl_target14;
    calibrationtable.xi(jj) = params{which_base_idx}{jj}.xi;
    calibrationtable.theta(jj) = params{which_base_idx}{jj}.t;
    calibrationtable.Gamma(jj) = params{which_base_idx}{jj}.Gamma;
    calibrationtable.sigmaP(jj) = params{which_base_idx}{jj}.sigma_P;
    calibrationtable.thetaP(jj) = params{which_base_idx}{jj}.theta_P;
end

calibrationtable.Properties.VariableUnits = {...
    'per thousand employees',...
    'per thousand employees',...
    '',...
    '',...
    '',...
    '',...
    '',...
    };

tab_name = 'calibrationtable.tex';
tab_preamble = '\\documentclass[../sections/appendix.tex]{subfiles}';
tab_FORMAT = '%4.2f';
table2latex(calibrationtable, [current_path save_path_tables tab_name], tab_FORMAT, tab_preamble);


disp('------------------- CALIBRATION FOR THE MODEL ---------------------')
disp(' ')
disp('Parameters:')
disp(' ')
disp(['$\varphi = ' num2str(params{which_base_idx}{1}.varphi,'%4.3g') '$'])
disp(['$\chi = ' num2str(params{which_base_idx}{1}.chi,'%4.3g') '$'])
disp(['$\A_F = ' num2str(params{which_base_idx}{1}.AF,'%4.3g') '$'])
disp(['$m = ' num2str(params{which_base_idx}{1}.E,'%4.3g') '$'])
disp(['$\psi_R$= ' num2str(params{which_base_idx}{1}.psi_k,'%4g') '$'])
disp(' ')
disp('-------------------------------------------------------------------')
disp(' ')
disp(calibrationtable)
disp(' ')


% COMPUTE SEMI-ELASTICITIES
disp(' ')
disp('UNTARGETED MOMENTS: Semi-Elasticity of Labor to Robot per Thousands Employees depending on parameter.')
change_pR_percent = [50, 75, 100, 125, 150];
for ii = 1:2:5
    which_base_idx_loop = ii;

    range_cs_mrts = (which_base_idx_loop - 1) * N_cs * 4  + (1:N_cs);
    range_cs_adjcost = (which_base_idx_loop - 1) * N_cs * 4  + (2*N_cs+1:3*N_cs);
    range_cs_pR_rel = (which_base_idx_loop - 1) * N_cs * 4  + (3*N_cs+1:4*N_cs);

    % MRTS effect
    L_el_mrts = -(L_cs(range_cs_mrts(2)) - L_cs(range_cs_mrts(1)) ) / (mrts_ratio_vec(2) - mrts_ratio_vec(1));
    RL_el_mrts = -1e3* (RL_cs(range_cs_mrts(2)) - RL_cs(range_cs_mrts(1)) ) / (mrts_ratio_vec(2) - mrts_ratio_vec(1));
    disp(['Labor Semi-Elasticity to Robot per Thousand Employees (MRTS) if Change in p_R ' num2str(change_pR_percent(ii)) '% of Actual Change: ' num2str(L_el_mrts / RL_el_mrts * 100) '%'])

    % AdjCost effect
    L_el_adj_cost = -(L_cs(range_cs_adjcost(2)) - L_cs(range_cs_adjcost(1)) ) / (adj_cost_ratio_vec(2) - adj_cost_ratio_vec(1));
    RL_el_adj_cost = -1e3* (RL_cs(range_cs_adjcost(2)) - RL_cs(range_cs_adjcost(1)) ) / (adj_cost_ratio_vec(2) - adj_cost_ratio_vec(1));
    disp(['Labor Semi-Elasticity to Robot per Thousand Employees (AdjCost) if Change in p_R ' num2str(change_pR_percent(ii)) '% of Actual Change: ' num2str(L_el_adj_cost / RL_el_adj_cost * 100) '%'])

    % price effect
    L_el_pR_rel = -(L_cs(range_cs_pR_rel(2)) - L_cs(range_cs_pR_rel(1)) ) / (pR_rel_future_vec(2) - pR_rel_future_vec(1));
    RL_el_pR_rel = -1e3* (RL_cs(range_cs_pR_rel(2)) - RL_cs(range_cs_pR_rel(1)) ) / (pR_rel_future_vec(2) - pR_rel_future_vec(1));
    disp(['Labor Semi-Elasticity to Robot per Thousand Employees (pR) if Change in p_R ' num2str(change_pR_percent(ii)) '% of Actual Change: ' num2str(L_el_pR_rel / RL_el_pR_rel * 100) '%'])

end
disp(' ')

% UNTARGETED MOMENT: Ratio of Adj Costs over Purchase Price in SS.
disp(' ')
disp('UNTARGETED MOMENTS: Ratio of adjustment costs to purchase prices')
num = 0;    % total adjustment costs in the economy
denom = 0;  % total purchase prices in the economy
for ii = 1:13
    disp(['Ratio in ' ifrStrings{ii} ' = ' ...
          num2str(out{3}{ii}.handles.Phi(out{3}{ii}.params.delta * out{3}{ii}.params.kmax) ...
                  / (out{3}{ii}.params.pR * out{3}{ii}.params.delta * out{3}{ii}.params.kmax))]);
    num = num + out{3}{ii}.handles.Phi(out{3}{ii}.params.delta * out{3}{ii}.params.kmax);
    denom = denom + out{3}{ii}.params.pR * out{3}{ii}.params.delta * out{3}{ii}.params.kmax;
end
disp(['Ratio in the entire economy = ' num2str(num/denom)])


% BAR PLOT OF 13-SECTORS EQUILIBRIUM EMPLOYMENT CHANGES
fig_format_bak = fig_format;
fig_format = 'png';

bar_plot_ifrStrings = categorical(ifrStrings);
leg = {['Scenario: $p_R/w = ' num2str(pR_rel_future) '$'], ...
       ['Scenario: $' num2str((mrts_ratio(1)-1)*100) '\%$ change in relative labor productivity'], ...
       ['Scenario: $' num2str((mrts_ratio(1)-1)*100) '\%$ change in relative labor productivity (Automotive Only)'], ...
       ['Scenario: $' num2str((adj_cost_ratio-1)*100) '\%$ change in robot adjustment costs']};

% FIGURE SETUP
fig_name = '4_cs_labor_by_sector';  % FIGURE NAME
fig_width = 650;	% Window Horizontal Size
fig_heigth = 400;   % Window Vertical Size
fig_posx = 160;     % Window position (Lower left corner)
fig_posy = 160;     % Window position (Lower left corner)

if docked
    figure()
else
    figure('Position', [fig_posx fig_posy fig_width fig_heigth])
end
colororder = get(gca, 'colororder');
hold on
bar(bar_plot_ifrStrings, plot_data, 'linewidth', 1.5, 'facealpha', 0.6, 'edgecolor', 'flat');
ax = gca;
ax.ColorOrderIndex = 1;
bar(bar_plot_ifrStrings, plot_data_PE, 'linewidth', .5, 'facealpha', .2, 'edgecolor', 'flat');
hold off
legend(leg, ...
       'interpreter', 'latex',...
       'orientation', 'vertical',...
       'location', 'southeast', ...
       'fontsize',14)
% title('Employment Losses to Automation',...
%       'interpreter', 'latex',...
%       'fontsize', 14)
ylabel('\% Change relative to 2014', ...
       'interpreter', 'latex', ...
       'fontsize',12)
a = get(gca,'XTickLabel');
set(gca,'XTickLabel',a,'fontsize',14)
    
if save_figs
    cd([current_path save_path_figures])
    fname = strcat(fig_name);
    set(gcf,'PaperPositionMode','auto');
    hgsave(fname);
    if strcmp(fig_format, 'eps')
        print('-depsc',fname,'-painters');
    elseif strcmp(fig_format, 'jpeg')
        print('-djpeg',fname,'-r256','-opengl');
    elseif strcmp(fig_format, 'png')
        print('-dpng',fname, '-r600');
    else
        error('Format specified incorrectly.')
    end
    if ~persistent_mode
        close
    end
    cd(current_path)
end

fig_format = fig_format_bak;


% PLOT OF CONTINUOUS COMPARATIVE STATICS
fig_format_bak = fig_format;
fig_format = 'png';

range_cs_mrts = (which_base_idx - 1) * N_cs * 4  + (1:N_cs);
range_cs_mrts_singlesector = (which_base_idx - 1) * N_cs * 4  + (N_cs+1:2*N_cs);
range_cs_adjcost = (which_base_idx - 1) * N_cs * 4  + (2*N_cs+1:3*N_cs);
range_cs_pR_rel = (which_base_idx - 1) * N_cs * 4  + (3*N_cs+1:4*N_cs);



% FIGURE SETUP: LABOR
fig_name = '5a_cs_labor_manysectors';  % FIGURE NAME
fig_width = 600;	% Window Horizontal Size
fig_heigth = 350;   % Window Vertical Size
fig_posx = 160;     % Window position (Lower left corner)
fig_posy = 160;     % Window position (Lower left corner)

if docked
    figure()
else
    figure('Position', [fig_posx fig_posy fig_width fig_heigth])
end


vecs = {L_cs, L_cs_PE};

subplot(2,2,1)
hold on
plot(mrts_ratio_vec * 100,L_cs(range_cs_mrts) * 100,'linewidth',2,'marker','.','color','k')


% plot(mrts_ratio_vec * 100,L_cs_PE(range_cs_mrts) * 100,'linewidth',2,'marker','.','color','r')
% Set up bands
range = range_cs_mrts;
x_vec = mrts_ratio_vec;
for width = 1:2 % band width
%     for jj = 1:2 % 1 = GE, 2 = PE
    for jj = 1
        if jj == 1
            facecolor = 'k';
        else
            facecolor = 'r';
        end
        Q_cs = vecs{jj};
        bands = [Q_cs(range - 4* N_cs * width) Q_cs(range + 4* N_cs * width)];
        % find nans 
        no_nan = logical(~isnan(bands(:,1)).* ~isnan(bands(:,2)));
        % Replace NaN with interpolated values at endpoints
        %bands = bands(no_nan,:);
        bands = [interp1(x_vec(no_nan), bands(no_nan,1), x_vec,'linear', 'extrap');
            interp1(x_vec(no_nan), bands(no_nan,2), x_vec,'linear', 'extrap');];
        %bands = [min(bands); max(bands)];
        nan_cs = isnan(L_cs(range));
        x_plot = x_vec(~nan_cs);
        fill([x_plot * 100, flip(x_plot * 100)], ...
             [max(bands(:,~nan_cs)) * 100, ...
             flip(max(0, min(bands(:,~nan_cs)) * 100))], ...
             facecolor, 'FaceAlpha', .3 - .1 * width, 'EdgeColor','none')
    end
end
set(gca, 'XDir','reverse')
hold off
title('MRTS$_{LR}$', 'interpreter', 'latex')
xlabel('Parameter (\% of 2014)', 'interpreter', 'latex')
ylabel('Labor (\% of 2014)', 'interpreter', 'latex')
ylim([-Inf 100])
grid on

subplot(2,2,2)
hold on
plot(mrts_ratio_singlesector_vec * 100,L_cs(range_cs_mrts_singlesector) * 100,'linewidth',2,'marker','.','color','k')
% plot(mrts_ratio_singlesector_vec * 100,L_cs_PE(range_cs_mrts_singlesector) * 100,'linewidth',2,'marker','.','color','r')
range = range_cs_mrts_singlesector;
x_vec = mrts_ratio_singlesector_vec;
for width = 1:2 % band width
%     for jj = 1:2 % 1 = GE, 2 = PE
    for jj = 1 % 1 = GE, 2 = PE
        if jj == 1
            facecolor = 'k';
        else
            facecolor = 'r';
        end
        Q_cs = vecs{jj};
        bands = [Q_cs(range - 4* N_cs * width)./Q_cs(range(1) - 4* N_cs * width), ... 
            Q_cs(range + 4* N_cs * width)./Q_cs(range(1) + 4* N_cs * width)];
        % find nans 
        idx_cross = find(bands(:,2)>bands(:,1),1);
        val_cross = bands(idx_cross,2);
        val_above_cross = bands(idx_cross-4,1);
        min_min = min(min(bands));
        max_min = max(min(bands));
        idx_low = find(bands(:,1)< val_cross,1);
        idx_high = find(bands(:,2)< val_cross,1);
        bands(idx_low -3: idx_low +2, 1) = NaN;
        bands(idx_high -3: idx_high +2, 2) = NaN;
        no_nan = logical(~isnan(bands(:,1)).* ~isnan(bands(:,2)));
        bands(no_nan, :) = [min(bands(no_nan,:),[],2), max(bands(no_nan,:),[],2)];
        % Replace NaN with interpolated values at endpoints
        %bands = bands(no_nan,:);
        bands = [interp1(x_vec(no_nan), bands(no_nan,1), x_vec,'spline', 'extrap');
        interp1(x_vec(no_nan), bands(no_nan,2), x_vec,'spline', 'extrap');];
        %bands = [min(bands); max(bands)];
        nan_cs = isnan(Q_cs(range));
        x_plot = x_vec(~nan_cs);
        fill([x_plot * 100, flip(x_plot * 100)], ...
             [max(bands(:,~nan_cs)) * 100, ...
             flip(max(0, min(bands(:,~nan_cs)) * 100))], ...
             facecolor, 'FaceAlpha', .3 - .1 * width, 'EdgeColor','none')
    end
end
set(gca, 'XDir','reverse')
hold off
title('MRTS$_{LR}$ (Automotive Only)', 'interpreter', 'latex')
xlabel('Parameter (\% of 2014)', 'interpreter', 'latex')
ylabel('Labor (\% of 2014)', 'interpreter', 'latex')
ylim([-Inf 100])
grid on

subplot(2,2,3)
hold on
plot(adj_cost_ratio_vec * 100,L_cs(range_cs_adjcost) * 100,'linewidth',2,'marker','.','color','k')
% plot(adj_cost_ratio_vec * 100,L_cs_PE(range_cs_adjcost) * 100,'linewidth',2,'marker','.','color','r')
range = range_cs_adjcost;
x_vec = adj_cost_ratio_vec;
for width = 1:2 % band width
%     for jj = 1:2 % 1 = GE, 2 = PE
    for jj = 1
        if jj == 1
            facecolor = 'k';
        else
            facecolor = 'r';
        end
        Q_cs = vecs{jj};
        bands = [Q_cs(range - 4* N_cs * width)./Q_cs(range(1) - 4* N_cs * width), ... 
            Q_cs(range + 4* N_cs * width)./Q_cs(range(1) + 4* N_cs * width)];
        % find nans 
        no_nan = logical(~isnan(bands(:,1)).* ~isnan(bands(:,2)));
        % Replace NaN with interpolated values at endpoints
        %bands = bands(no_nan,:);
        bands = [interp1(x_vec(no_nan), bands(no_nan,1), x_vec,'linear', 'extrap');
            interp1(x_vec(no_nan), bands(no_nan,2), x_vec,'linear', 'extrap');];
        %bands = [min(bands); max(bands)];
        nan_cs = isnan(L_cs(range));
        x_plot = x_vec(~nan_cs);
        fill([x_plot * 100, flip(x_plot * 100)], ...
             [max(bands(:,~nan_cs)) * 100, ...
             flip(max(0, min(bands(:,~nan_cs)) * 100))], ...
             facecolor, 'FaceAlpha', .3 - .1 * width, 'EdgeColor','none')
    end
end
set(gca, 'XDir','reverse')
hold off
title('$\psi_R$', 'interpreter', 'latex')
xlabel('Parameter (\% of 2014)', 'interpreter', 'latex')
ylabel('Labor (\% of 2014)', 'interpreter', 'latex')
ylim([-Inf 100])
grid on
% legend({'General Equilibrium', 'Partial Equilibrium'},...
%        'interpreter', 'latex',...
%        'orientation', 'vertical',...
%        'location','southwest',...
%        'fontsize', 14);

subplot(2,2,4)
hold on
plot(pR_rel_future_vec * 100,L_cs(range_cs_pR_rel) * 100,'linewidth',2,'marker','.','color','k')
% plot(pR_rel_future_vec * 100,L_cs_PE(range_cs_pR_rel) * 100,'linewidth',2,'marker','.','color','r')
range = range_cs_pR_rel;
x_vec = pR_rel_future_vec;
for width = 1:2 % band width
%     for jj = 1:2 % 1 = GE, 2 = PE
    for jj = 1
        if jj == 1
            facecolor = 'k';
        else
            facecolor = 'r';
        end
        Q_cs = vecs{jj};
        bands = [Q_cs(range - 4* N_cs * width)./Q_cs(range(1) - 4* N_cs * width), ... 
            Q_cs(range + 4* N_cs * width)./Q_cs(range(1) + 4* N_cs * width)];
        % find nans 
        no_nan = logical(~isnan(bands(:,1)).* ~isnan(bands(:,2)));
        % Replace NaN with interpolated values at endpoints
        %bands = bands(no_nan,:);
        bands = [interp1(x_vec(no_nan), bands(no_nan,1), x_vec,'linear', 'extrap');
            interp1(x_vec(no_nan), bands(no_nan,2), x_vec,'linear', 'extrap');];
        %bands = [min(bands); max(bands)];
        nan_cs = isnan(L_cs(range));
        x_plot = x_vec(~nan_cs);
        fill([x_plot * 100, flip(x_plot * 100)], ...
             [max(bands(:,~nan_cs)) * 100, ...
             flip(max(0, min(bands(:,~nan_cs)) * 100))], ...
             facecolor, 'FaceAlpha', .3 - .1 * width, 'EdgeColor','none')
    end
end
set(gca, 'XDir','reverse')
hold off
title('$p_R$', 'interpreter', 'latex')
xlabel('Parameter (\% of 2014)', 'interpreter', 'latex')
ylabel('Labor (\% of 2014)', 'interpreter', 'latex')
ylim([-Inf 100])
grid on  

if save_figs
    cd([current_path save_path_figures])
    fname = strcat(fig_name);
    set(gcf,'PaperPositionMode','auto');
    hgsave(fname);
    if strcmp(fig_format, 'eps')
        print('-depsc',fname,'-painters');
    elseif strcmp(fig_format, 'jpeg')
        print('-djpeg',fname,'-r256','-opengl');
    elseif strcmp(fig_format, 'png')
        print('-dpng',fname, '-r600');
    else
        error('Format specified incorrectly.')
    end
    if ~persistent_mode
        close
    end
    cd(current_path)
end

% FIGURE SETUP: LABOR SHARE
fig_name = '5b_cs_laborshare_manysectors';  % FIGURE NAME
fig_width = 600;	% Window Horizontal Size
fig_heigth = 350;   % Window Vertical Size
fig_posx = 160;     % Window position (Lower left corner)
fig_posy = 160;     % Window position (Lower left corner)

if docked
    figure()
else
    figure('Position', [fig_posx fig_posy fig_width fig_heigth])
end


vecs = {LS_cs, LS_cs_PE};

subplot(2,2,1)
hold on
plot(mrts_ratio_vec * 100,LS_cs(range_cs_mrts)./LS_cs(range_cs_mrts(1)) * 100,'linewidth',2,'marker','.','color','k')
% plot(mrts_ratio_vec * 100,LS_cs_PE(range_cs_mrts)./LS_cs_PE(range_cs_mrts(1)) * 100,'linewidth',2,'marker','.','color','r')
range = range_cs_mrts;
x_vec = mrts_ratio_vec;
for width = 1:2 % band width
%     for jj = 1:2 % 1 = GE, 2 = PE
    for jj = 1
        if jj == 1
            facecolor = 'k';
        else
            facecolor = 'r';
        end
        Q_cs = vecs{jj};
        bands = [Q_cs(range - 4* N_cs * width)./Q_cs(range(1) - 4* N_cs * width), ... 
            Q_cs(range + 4* N_cs * width)./Q_cs(range(1) + 4* N_cs * width)];
        % find nans 
        no_nan = logical(~isnan(bands(:,1)).* ~isnan(bands(:,2)));
        % Replace NaN with interpolated values at endpoints
        %bands = bands(no_nan,:);
        bands = [interp1(x_vec(no_nan), bands(no_nan,1), x_vec,'linear', 'extrap');
            interp1(x_vec(no_nan), bands(no_nan,2), x_vec,'linear', 'extrap');];
        %bands = [min(bands); max(bands)];
        nan_cs = isnan(Q_cs(range));
        x_plot = x_vec(~nan_cs);
        fill([x_plot * 100, flip(x_plot * 100)], ...
             [max(bands(:,~nan_cs)) * 100, ...
             flip(max(0, min(bands(:,~nan_cs)) * 100))], ...
             facecolor, 'FaceAlpha', .3 - .1 * width, 'EdgeColor','none')
    end
end
set(gca, 'XDir','reverse')
hold off
title('MRTS$_{LR}$', 'interpreter', 'latex')
xlabel('Parameter (\% of 2014)', 'interpreter', 'latex')
ylabel('Labor Share (\% of 2014)', 'interpreter', 'latex')
ylim([-Inf 100])
grid on

subplot(2,2,2)
hold on
plot(mrts_ratio_singlesector_vec * 100,LS_cs(range_cs_mrts_singlesector)./LS_cs(range_cs_mrts_singlesector(1)) * 100,'linewidth',2,'marker','.','color','k')
% plot(mrts_ratio_singlesector_vec * 100,LS_cs_PE(range_cs_mrts_singlesector)./LS_cs_PE(range_cs_mrts_singlesector(1)) * 100,'linewidth',2,'marker','.','color','r')
range = range_cs_mrts_singlesector;
x_vec = mrts_ratio_singlesector_vec;
for width = 1:2 % band width
%     for jj = 1:2 % 1 = GE, 2 = PE
    for jj = 1
        if jj == 1
            facecolor = 'k';
        else
            facecolor = 'r';
        end
        Q_cs = vecs{jj};
        bands = [Q_cs(range - 4* N_cs * width)./Q_cs(range(1) - 4* N_cs * width), ... 
            Q_cs(range + 4* N_cs * width)./Q_cs(range(1) + 4* N_cs * width)];
        % find nans 
        % Replace NaN with interpolated values at endpoints
        %bands are upper and lower bounds, but they cross because the slope
        %changes at different points. Convexify the curves to embrace all
        %the values that the percentages can take
        idx_cross = find(bands(:,2)>bands(:,1),1);
        val_cross = bands(idx_cross,2);
        val_above_cross = bands(idx_cross-4,1);
        idx_low = find(bands(:,1)< val_cross,1);
        idx_high = find(bands(:,2)< val_cross,1);
        bands(idx_low -3: idx_low +2, 1) = NaN;
        bands(idx_high -3: idx_high +2, 2) = NaN;
        no_nan = logical(~isnan(bands(:,1)).* ~isnan(bands(:,2)));
        bands(no_nan, :) = [min(bands(no_nan,:),[],2), max(bands(no_nan,:),[],2)];
        % Replace NaN with interpolated values at endpoints
        %bands = bands(no_nan,:);
        bands = [interp1(x_vec(no_nan), bands(no_nan,1), x_vec,'spline', 'extrap');
        interp1(x_vec(no_nan), bands(no_nan,2), x_vec,'spline', 'extrap');];
        %bands = [min(bands); max(bands)];
        nan_cs = isnan(Q_cs(range));
        x_plot = x_vec(~nan_cs);
        fill([x_plot * 100, flip(x_plot * 100)], ...
             [max(bands(:,~nan_cs)) * 100, ...
             flip(max(0, min(bands(:,~nan_cs)) * 100))], ...
             facecolor, 'FaceAlpha', .3 - .1 * width, 'EdgeColor','none')
    end
end
set(gca, 'XDir','reverse')
hold off
title('MRTS$_{LR}$ (Automotive Only)', 'interpreter', 'latex')
xlabel('Parameter (\% of 2014)', 'interpreter', 'latex')
ylabel('Labor Share (\% of 2014)', 'interpreter', 'latex')
ylim([-Inf 100])
grid on

subplot(2,2,3)
hold on
plot(adj_cost_ratio_vec * 100,LS_cs(range_cs_adjcost)./LS_cs(range_cs_adjcost(1)) * 100,'linewidth',2,'marker','.','color','k')
% plot(adj_cost_ratio_vec * 100,LS_cs_PE(range_cs_adjcost)./LS_cs_PE(range_cs_adjcost(1)) * 100,'linewidth',2,'marker','.','color','r')
range = range_cs_adjcost;
x_vec = adj_cost_ratio_vec;
for width = 1:2 % band width
%     for jj = 1:2 % 1 = GE, 2 = PE
    for jj = 1 % 1 = GE, 2 = PE
        if jj == 1
            facecolor = 'k';
        else
            facecolor = 'r';
        end
        Q_cs = vecs{jj};
        bands = [Q_cs(range - 4* N_cs * width)./Q_cs(range(1) - 4* N_cs * width), ... 
            Q_cs(range + 4* N_cs * width)./Q_cs(range(1) + 4* N_cs * width)];
        % find nans 
        no_nan = logical(~isnan(bands(:,1)).* ~isnan(bands(:,2)));
        % Replace NaN with interpolated values at endpoints
        %bands = bands(no_nan,:);
        bands = [interp1(x_vec(no_nan), bands(no_nan,1), x_vec,'linear', 'extrap');
            interp1(x_vec(no_nan), bands(no_nan,2), x_vec,'linear', 'extrap');];
        %bands = [min(bands); max(bands)];
        nan_cs = isnan(Q_cs(range));
        x_plot = x_vec(~nan_cs);
        fill([x_plot * 100, flip(x_plot * 100)], ...
             [max(bands(:,~nan_cs)) * 100, ...
             flip(max(0, min(bands(:,~nan_cs)) * 100))], ...
             facecolor, 'FaceAlpha', .3 - .1 * width, 'EdgeColor','none')
    end
end
set(gca, 'XDir','reverse')
hold off
title('$\psi_R$', 'interpreter', 'latex')
xlabel('Parameter (\% of 2014)', 'interpreter', 'latex')
ylabel('Labor Share (\% of 2014)', 'interpreter', 'latex')
ylim([-Inf 100])
grid on
% legend({'General Equilibrium', 'Partial Equilibrium'},...
%        'interpreter', 'latex',...
%        'orientation', 'vertical',...
%        'location','southwest',...
%        'fontsize', 14);
   
subplot(2,2,4)
hold on
plot(pR_rel_future_vec * 100,LS_cs(range_cs_pR_rel) ./LS_cs(range_cs_pR_rel(1)) * 100,'linewidth',2,'marker','.','color','k')
% plot(pR_rel_future_vec * 100,LS_cs_PE(range_cs_pR_rel)./LS_cs_PE(range_cs_pR_rel(1)) * 100,'linewidth',2,'marker','.','color','r')
range = range_cs_pR_rel;
x_vec = pR_rel_future_vec;
for width = 1:2 % band width
%     for jj = 1:2 % 1 = GE, 2 = PE
    for jj = 1
        if jj == 1
            facecolor = 'k';
        else
            facecolor = 'r';
        end
        Q_cs = vecs{jj};
        bands = [Q_cs(range - 4* N_cs * width)./Q_cs(range(1) - 4* N_cs * width), ... 
            Q_cs(range + 4* N_cs * width)./Q_cs(range(1) + 4* N_cs * width)];
        % find nans 
        no_nan = logical(~isnan(bands(:,1)).* ~isnan(bands(:,2)));
        % Replace NaN with interpolated values at endpoints
        %bands = bands(no_nan,:);
        bands = [interp1(x_vec(no_nan), bands(no_nan,1), x_vec,'linear', 'extrap');
            interp1(x_vec(no_nan), bands(no_nan,2), x_vec,'linear', 'extrap');];
        %bands = [min(bands); max(bands)];
        nan_cs = isnan(Q_cs(range));
        x_plot = x_vec(~nan_cs);
        fill([x_plot * 100, flip(x_plot * 100)], ...
             [max(bands(:,~nan_cs)) * 100, ...
             flip(max(0, min(bands(:,~nan_cs)) * 100))], ...
             facecolor, 'FaceAlpha', .3 - .1 * width, 'EdgeColor','none')
    end
end
set(gca, 'XDir','reverse')
hold off
title('$p_R$', 'interpreter', 'latex')
xlabel('Parameter (\% of 2014)', 'interpreter', 'latex')
ylabel('Labor Share (\% of 2014)', 'interpreter', 'latex')
ylim([-Inf 100])
grid on


if save_figs
    cd([current_path save_path_figures])
    fname = strcat(fig_name);
    set(gcf,'PaperPositionMode','auto');
    hgsave(fname);
    if strcmp(fig_format, 'eps')
        print('-depsc',fname,'-painters');
    elseif strcmp(fig_format, 'jpeg')
        print('-djpeg',fname,'-r256','-opengl');
    elseif strcmp(fig_format, 'png')
        print('-dpng',fname, '-r600');
    else
        error('Format specified incorrectly.')
    end
    if ~persistent_mode
        close
    end
    cd(current_path)
end

fig_format = fig_format_bak;

%% CLEAR (SOME) VARIABLES IN MEMORY AND CONTINUE %%
clearvars -EXCEPT save_figs fig_format save_path_tables save_path_figures current_path docked persistent_mode

%% 5) ***** Extension with Labor Adjustment Costs *****

% Set parameters
params = SetParameters();
params.N_l = 30;
params.N_k = 30;
params.N_p = 15;

params.W = 0.268544000136280;
params.pR = params.W * 25;

params.theta_P = .8789605;             % see HamiltonOUEstimatorsIFR_total.dta
params.sigma_P = .1418431;             % see HamiltonOUEstimatorsIFR_total.dta
params.t = .2978469;
params.E = 0;
params.Gamma = 0.652862965992859;
% params.psi_k = 5;
params.psi_k = 499;
psilow = .001;
psihigh = 5;
params.sep = .313; % from BLS data (JOLTS)

% Wage and Price of Robots
W = params.W;
% params.pR = 25;
% pR = params.pR;
pR = 1.4852; % obtained from solvePathpR.m
params.LaborUtilization = 'No';
params.print_conv = true;

params.lmax = .7;

% Solve Convergence

%{

% Setup experiment for immediate convergence
transition.D = .1;
T = 100;
NT = T / transition.D;
pR_infty = 0;
transition.par_str = 'pR';
% transition.par_path = [linspace(pR, 0, floor(0.2 * NT)) zeros(1,NT - floor(0.2 * NT))];
transition.par_path = [pR, pR_infty * ones(1,NT-1)];


% Run case with high rigidity
params.psi_l = psihigh;
disp(' ')
disp([' IMMEDIATE CONVERGENCE TO pR = ' num2str(pR_infty)])
disp(['Running case with psi_l = ' num2str(params.psi_l) '...'])
out_highrigidity = SolveTransition_rigidlabor(transition, params);

% Run case with low rigidity
params.psi_l = psilow;
disp(' ')
disp(['Running case with psi_l = ' num2str(params.psi_l) '...'])
out_lowrigidity = SolveTransition_rigidlabor(transition, params);
%}

% Setup experiment for linear convergence
transition_linear.D = 1;
T = 150;
NT = T / transition_linear.D;

years_before = 20;
pR_infty = 0;
% pR_infty = 0.711163645915419;
[alpha_path, pR_infty] = solvePathpR(pR_infty);
transition_linear.par_str = 'pR';
% transition.par_path = [pR; pR_infty * ones(NT-1,1)];
transition_linear.par_path = pR_infty +  (pR - pR_infty) .* ...
                      exp( - alpha_path .* (-years_before/transition_linear.D:transition_linear.D:T-years_before/transition_linear.D)');
% transition_linear.T_note = 10;

% pR_infty = 0;
% pR_lenthlinearconvergence = .2 * NT;
% transition_linear.par_str = 'pR';
% transition_linear.par_path = [linspace(pR, pR_infty, floor(pR_lenthlinearconvergence)) ...
%         pR_infty * ones(1,NT - floor(pR_lenthlinearconvergence))];

% Run case with high rigidity
params.psi_l = psihigh;
disp(' ')
disp(['Running case with psi_l = ' num2str(params.psi_l) '...'])
out_highrigidity_linear = SolveTransition_rigidlabor(transition_linear, params);

% Run case with low rigidity
params.psi_l = psilow;
disp(' ')
disp(['Running case with psi_l = ' num2str(params.psi_l) '...'])
out_lowrigidity_linear = SolveTransition_rigidlabor(transition_linear, params);

% PLOTS OF TRANSITIONS

% EXPONENTIAL CONVERGENCE of pR to 0

% FIGURE SETUP
fig_name = '9_transition_rigidlabor';  % FIGURE NAME
fig_width = 650;	% Window Horizontal Size
fig_heigth = 450;   % Window Vertical Size
fig_posx = 160;     % Window position (Lower left corner)
fig_posy = 160;     % Window position (Lower left corner)

if docked
    figure()
else
    figure('Position', [fig_posx fig_posy fig_width fig_heigth])
end

base_year = 2010;

subplot(3,1,1)
plot(base_year + ( -years_before:transition_linear.D:T-years_before), ...
     out_highrigidity_linear.L_path ./ out_highrigidity_linear.L_path(1) * 100,...
     'linewidth', 2, ...
     'color', 'r')
hold on
plot(base_year + ( -years_before:transition_linear.D:T-years_before), ...
     out_lowrigidity_linear.L_path ./ out_lowrigidity_linear.L_path(1) * 100,...
     'linewidth', 2, ...
     'color', 'k')
hold off
xlim([1990, 2060])
ylim([.99* min([out_highrigidity_linear.L_path ./ out_highrigidity_linear.L_path(1) * 100;...
             out_lowrigidity_linear.L_path ./ out_lowrigidity_linear.L_path(1) * 100])...
      1.01*max([out_highrigidity_linear.L_path ./ out_highrigidity_linear.L_path(1) * 100;...
             out_lowrigidity_linear.L_path ./ out_lowrigidity_linear.L_path(1) * 100])])
legend({['High Rigidity $\psi_L = ' num2str(psihigh) '$'], ...
        ['Low Rigidity $\psi_L = ' num2str(psilow) '$']}, ...
        'interpreter', 'latex',...
        'location', 'southwest')
title('Aggregate Labor', 'Interpreter', 'latex')
xlabel('Time (years)', 'interpreter', 'latex')
ylabel('Value relative to $1990$ (\%)', 'interpreter', 'latex')
grid on

box on

subplot(3,1,2)
plot(base_year + ( -years_before:transition_linear.D:T-years_before), ...
     out_highrigidity_linear.RL_path,...
     'linewidth', 2, ...
     'color', 'r')
hold on
plot(base_year + ( -years_before:transition_linear.D:T-years_before), ...
     out_lowrigidity_linear.RL_path,...
     'linewidth', 2, ...
     'color', 'k')
hold off
xlim([1990, 2060])
ylim([0 1.01*max([out_highrigidity_linear.RL_path;...
             out_lowrigidity_linear.RL_path])])
% legend({['High Rigidity $\psi_L = ' num2str(psihigh) '$'], ...
%         ['Low Rigidity $\psi_L = ' num2str(psilow) '$']}, ...
%         'interpreter', 'latex')
title('Robots per Thousand Employees', 'Interpreter', 'latex')
xlabel('Time (years)', 'interpreter', 'latex')
ylabel('$R/L \times 1000$', 'interpreter', 'latex')
grid on

box on

subplot(3,1,3)
plot(base_year + ( -years_before:transition_linear.D:T-years_before), ...
     transition_linear.par_path ./ transition_linear.par_path(1) * 100, ...
     'linewidth', 2, ...
     'color', 'k')
ylim([0 1.1*100])
xlim([1990, 2060])
%xlim([-.1*T, T])
title('Price of Robots', 'interpreter', 'latex')
xlabel('Time (years)', 'interpreter', 'latex')
ylabel('Value relative to $1990$ (\%)', 'interpreter', 'latex')
grid on

box on

if save_figs
    cd([current_path save_path_figures])
    fname = strcat(fig_name);
    set(gcf,'PaperPositionMode','auto');
    hgsave(fname);
    if strcmp(fig_format, 'eps')
        print('-depsc',fname,'-painters');
    elseif strcmp(fig_format, 'jpeg')
        print('-djpeg',fname,'-r256','-opengl');
    elseif strcmp(fig_format, 'png')
        print('-dpng',fname);
    else
        error('Format specified incorrectly.')
    end
    if ~persistent_mode
        close
    end
    cd(current_path)
end

%% CLEAR (SOME) VARIABLES IN MEMORY AND CONTINUE %%
clearvars -EXCEPT save_figs fig_format save_path_tables save_path_figures current_path docked persistent_mode

%% 6) ***** Model with Linear Adjustment Costs on Robots *****

% Set Parameters
params = SetParameters();

params.W = 0.268544000136280 * 2;
params.pR = params.W * 10;

params.theta_P = .8789605;             % see HamiltonOUEstimatorsIFR_total.dta
params.sigma_P = .1418431;             % see HamiltonOUEstimatorsIFR_total.dta
params.t = .2978469;
params.E = 0;
params.Gamma = 0.652862965992859;      % see single-sector calibration
% params.psi_k = 5;

% Display advancement
params.print_conv = true;

% Set Grid Size
params.N_k = 400;
params.N_p = 400;
% params.pmax = 2.5;

% Set Adjustment Cost Parameters
params.psi_plus = 2.5;
params.psi_minus = 0;
params.kmax = 500;
params.GridKGamma = 1.15;

params.psi_plus = 1.2;
params.psi_minus = 0;
params.kmax = 10;
params.GridKGamma = 1.25;

% params.pR = 5;
% params.W = .5;

% Run Labor Demand Script for Linear Adj Costs à la Bentolila Bertola
[~, out] = LaborDemand_linear(params.W, params.pR, params);

% PLOT OF STATIONARY DISTRIBUTION AND INACTION REGION

% FIGURE SETUP
fig_name = '10_stationary_distribution_linearAdjCosts';  % FIGURE NAME
fig_width = 650;	% Window Horizontal Size
fig_heigth = 350;   % Window Vertical Size
fig_posx = 160;     % Window position (Lower left corner)
fig_posy = 160;     % Window position (Lower left corner)

if docked
    figure()
else
    figure('Position', [fig_posx fig_posy fig_width fig_heigth])
end

hold on
contour(out.p_grid, out.k_grid, ...
    reshape(out.gg .* out.inaction_region, out.params.N_k, out.params.N_p),...
    15, ...
    'linewidth', 2)

p_interp = linspace(out.p_grid(1), out.p_grid(end), 15);
y_interp_1 = interp1(out.p_grid, out.k_grid(out.inaction_indices(:,1)), p_interp);
y_interp_2 = interp1(out.p_grid, out.k_grid(out.inaction_indices(:,2)), p_interp);

plot(p_interp, y_interp_1,...
     'linewidth', 3, 'linestyle', '-', 'color', 'r')
plot(p_interp, y_interp_2,...
     'linewidth', 3, 'linestyle', ':', 'color', 'r')

% plot(out.p_grid, out.k_grid(out.inaction_indices(:,1)),...
%      'linewidth', 2, 'linestyle', '-', 'color', 'r')
% plot(out.p_grid, out.k_grid(out.inaction_indices(:,2)),...
%      'linewidth', 2, 'linestyle', ':', 'color', 'r')
fplot(out.handles.R_bar, [out.params.pmin out.params.pmax], ...
      'linewidth', 3, ...
      'color', 'k')
text(0.25,.9,'Sell Robots', 'interpreter', 'latex', 'fontsize', 14)
text(1.3,.25,'Install Robots', 'interpreter', 'latex', 'fontsize', 14)
text(1.3,.68,'Inaction Region', 'interpreter', 'latex', 'fontsize', 14)
line([params.P, params.P], [0 1], 'linestyle', ':', 'color', 'k', 'linewidth', 1)
hold off
xlabel('Revenue-Shifter $z$', 'interpreter', 'latex')
ylabel('Robot Stock $R$', 'interpreter', 'latex')
yticks([ ])
xticks([0 params.P])
xticklabels({'0', '$z^d$'})
a = get(gca,'XTickLabel');
set(gca,'XTickLabel',a,'fontsize',14)
legend({'Stationary Distribution',...
        '$R^{\star}_{\textrm{inv}}(z)$', ...
        '$R^{\star}_{\textrm{disinv}}(z)$', ...
        '$\bar{R}(z)$'}, ...
        'interpreter', 'latex', 'fontsize', 14, 'box', 'off')
% ylim([0, out.kmaxgg* 1.5])
xlim([-Inf, .75*out.p_grid(end)])
ylim([0, 1])
box on

if save_figs
    cd([current_path save_path_figures])
    fname = strcat(fig_name);
    set(gcf,'PaperPositionMode','auto');
    hgsave(fname);
    if strcmp(fig_format, 'eps')
        print('-depsc',fname,'-painters');
    elseif strcmp(fig_format, 'jpeg')
        print('-djpeg',fname,'-r256','-opengl');
    elseif strcmp(fig_format, 'png')
        print('-dpng',fname);
    else
        error('Format specified incorrectly.')
    end
    if ~persistent_mode
        close
    end
    cd(current_path)
end

%% CLEAR (SOME) VARIABLES IN MEMORY AND CONTINUE %%
clearvars -EXCEPT save_figs fig_format save_path_tables save_path_figures current_path docked persistent_mode

%% 7) ***** One-Sector Transitions *****

% Set parameters
params = SetParameters();
params.N_k = 75;
params.N_p = 50;

% DRS parameter (from data)
params.E = 0;
params.theta_P = .8789605;             % see HamiltonOUEstimatorsIFR_total.dta
params.sigma_P = .1418431;             % see HamiltonOUEstimatorsIFR_total.dta
params.t = .2978469;                   % see ThetaProdIFR_all.dta

opts = delimitedTextImportOptions("NumVariables", 5);
opts.DataLines = [2, Inf];
opts.Delimiter = ",";
opts.VariableNames = ["ifrCode", ...
                      "emp_share_04", ...
                      "emp_share_07", ...
                      "emp_share_10", ...
                      "emp_share_14"];
opts.VariableTypes = ["string",...
                      "double", ...
                      "double", ...
                      "double", ...
                      "double"];
opts = setvaropts(opts, 1, "WhitespaceRule", "preserve");
opts = setvaropts(opts, 1, "EmptyFieldRule", "auto");
opts.ExtraColumnsRule = "ignore";
opts.EmptyLineRule = "read";
shares = readtable([current_path '/data/emp_BEA_IFR.csv'], opts);
clear opts

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

Statistics = join(Statistics, shares, 'Keys', 'ifrCode');

RL_04 = Statistics.apr_lv04' * Statistics.emp_share_04;
RL_07 = Statistics.apr_lv07' * Statistics.emp_share_07;
RL_10 = Statistics.apr_lv10' * Statistics.emp_share_10;
RL_14 = Statistics.apr_lv14' * Statistics.emp_share_14;

% Relative price of robots:
% pR_rel04 = 2.3224;
% pR_rel05 = 2.2971;
% pR_rel07 = 1.7865;
% pR_rel10 = 1.4348;
pR_rel14 = 1.0209;
pR = 1.4852; % obtained from solvePathpR.m




% Solve Convergence
% Setup experiment for convergence over time
transition.D = 1;
T = 100;
% NT = T / transition.D;
pR_infty = 0;
% pR_infty = 0.711163645915419;
[alpha_path, pR_infty] = solvePathpR(pR_infty);
transition.par_str = 'pR';
% transition.par_path = [pR; pR_infty * ones(NT-1,1)];
transition.par_path = pR_infty +  (pR - pR_infty) .* ...
                      exp( - alpha_path .* (-10/transition.D:transition.D:T-10/transition.D-1)');
transition.T_note = 10;

% Run case
% disp(' ')
% disp([' IMMEDIATE CONVERGENCE TO pR = ' num2str(pR_infty)])
% out = SolveTransition(transition, params);
% x0 = [0.825; 5];
x0 = [0.652862965992859;
      4.997927931801950e+02];

  

% t_ancredi = 4 / transition.D;
GammaUpperbar = 1;
targets = [RL_04; RL_07; RL_10; RL_14];
fun = @(x) SolveTransitionWrapper([GammaUpperbar ./ (1 + exp(-x(1))); ...
                                   exp(x(2))],...
                                   targets, transition, params);
opts = optimset('display', 'iter');

x = fminsearch(fun, [log(x0(1) ./ (GammaUpperbar - x0(1)));...
                     log(x0(2))], opts);


gamma = GammaUpperbar ./ (1 + exp( - x(1) ));
psi_k = exp(x(2));

params.Gamma = gamma;
params.psi_k =  psi_k;

out = SolveTransition(transition, params);
params.chi = out.out_initial.params.chi;

% Run Continuous Comparative Statics Exercise

N_cs = 48;
mrts_ratio_vec = linspace(1, 0.05, N_cs);
% mrts_ratio_singlesector_vec = linspace(1, 0.01, N_cs);
adj_cost_ratio_vec = linspace(1, 0.005, N_cs);
% pR_rel_future_vec = linspace(1, 0, N_cs);

parameters_cs = [[mrts_ratio_vec'; ones(N_cs,1)], ...
                 [ones(N_cs,1); adj_cost_ratio_vec']];
out_cs = cell(N_cs * 2, 1);
% out_cs_PE = cell(N_cs * 2, 1);
L_cs = zeros(N_cs * 2, 1);
% L_cs_PE = L_cs;

LS_cs = L_cs;
% LS_cs_PE = L_cs;

RL_cs = L_cs;

disp('SMOOTH COMPARATIVE STATICS ON PARAMETERS OF INTEREST...')
parfor ii = 1:N_cs*2
    params_cs = params;
    params_cs.pR = pR_rel14;
    
    params_cs.Gamma = parameters_cs(ii,1) .* params.Gamma ./ ...
                (1 - params.Gamma + parameters_cs(ii,1) .* params.Gamma);
    params_cs.psi_k = parameters_cs(ii,2) * params.psi_k;
    
    out = generalEquilibrium_onesector(params_cs, false);
    
    L_cs(ii) = out.L;
    RL_cs(ii) = out.RL;
    LS_cs(ii) = out.W * L_cs(ii) / out.GE.q;

    out_cs{ii}  = out;
    
    disp(['Completed Iteration: ' num2str(ii) '.'])
    
end

% PLOT OF CONTINUOUS COMPARATIVE STATICS

% FIGURE SETUP: LABOR
fig_name = '7_cs_onesector';  % FIGURE NAME
fig_width = 650;	% Window Horizontal Size
fig_heigth = 200;   % Window Vertical Size
fig_posx = 160;     % Window position (Lower left corner)
fig_posy = 160;     % Window position (Lower left corner)

if docked
    figure()
else
    figure('Position', [fig_posx fig_posy fig_width fig_heigth])
end

subplot(1,2,1)
hold on
plot(mrts_ratio_vec * 100,L_cs(1:N_cs) ./ L_cs(1) * 100,'linewidth',2,'marker','.','color','k')
% plot(mrts_ratio_vec * 100,L_cs_PE(1:N_cs) * 100,'linewidth',2,'marker','.','color','r')
set(gca, 'XDir','reverse')
hold off
title('MRTS$_{LR}$', 'interpreter', 'latex')
xlabel('Parameter (\% of 2014 Value)', 'interpreter', 'latex')
ylabel('Labor (\% of 2014 Value)', 'interpreter', 'latex')
grid on
% ylim([0,100])

% subplot(2,2,2)
% hold on
% plot(mrts_ratio_singlesector_vec * 100,L_cs(N_cs+1:2*N_cs) * 100,'linewidth',2,'marker','.','color','k')
% % plot(mrts_ratio_singlesector_vec * 100,L_cs_PE(N_cs+1:2*N_cs) * 100,'linewidth',2,'marker','.','color','r')
% set(gca, 'XDir','reverse')
% hold off
% title('MRTS$_{LR}$ (Automotive Only)', 'interpreter', 'latex')
% xlabel('Parameter (\% of 2014 Value)', 'interpreter', 'latex')
% ylabel('Labor (\% of 2014 Value)', 'interpreter', 'latex')
% grid on

subplot(1,2,2)
hold on
plot(adj_cost_ratio_vec * 100,L_cs(N_cs+1:end) ./ L_cs(1+N_cs) * 100,'linewidth',2,'marker','.','color','k')
% plot(adj_cost_ratio_vec * 100,L_cs_PE(2*N_cs+1:3*N_cs) * 100,'linewidth',2,'marker','.','color','r')
set(gca, 'XDir','reverse')
hold off
title('$\psi_R$', 'interpreter', 'latex')
xlabel('Parameter (\% of 2014 Value)', 'interpreter', 'latex')
ylabel('Labor (\% of 2014 Value)', 'interpreter', 'latex')
grid on
% ylim([0,100])

% subplot(2,2,4)
% hold on
% plot(pR_rel_future_vec * 100,L_cs(3*N_cs+1:4*N_cs) * 100,'linewidth',2,'marker','.','color','k')
% % plot(pR_rel_future_vec * 100,L_cs_PE(3*N_cs+1:4*N_cs) * 100,'linewidth',2,'marker','.','color','r')
% set(gca, 'XDir','reverse')
% hold off
% title('$p_R$', 'interpreter', 'latex')
% xlabel('Parameter (\% of 2014 Value)', 'interpreter', 'latex')
% ylabel('Labor (\% of 2014 Value)', 'interpreter', 'latex')
% grid on
% legend({'General Equilibrium', 'Partial Equilibrium'},...
%        'interpreter', 'latex',...
%        'orientation', 'vertical',...
%        'location','southwest',...
%        'fontsize', 11);  

if save_figs
    cd([current_path save_path_figures])
    fname = strcat(fig_name);
    set(gcf,'PaperPositionMode','auto');
    hgsave(fname);
    if strcmp(fig_format, 'eps')
        print('-depsc',fname,'-painters');
    elseif strcmp(fig_format, 'jpeg')
        print('-djpeg',fname,'-r256','-opengl');
    elseif strcmp(fig_format, 'png')
        print('-dpng',fname);
    else
        error('Format specified incorrectly.')
    end
    if ~persistent_mode
        close
    end
    cd(current_path)
end

% PLOT OF TRANSITION PATHS

t_pR = [1995, 1996, 2005, 2010, 2014, 2017]';
w_pR = [25492, 25492, 29890, 31073, 35447, 38070];
p_pR = [133382, 121452, 68659, 46150, 31776, 27074]' ./ w_pR';
t_RL = [2004, 2007, 2010, 2014];
RL_RL = [RL_04, RL_07,RL_10,RL_14];
D = transition.D;
T = length(transition.par_path) * D;
t_bag = D:D:T;
t_sucre = 1:find(t_bag == 96, 1);
n_abbruzzi = 3;


% FIGURE SETUP: PATHS
fig_name = '6_transition_onesector';  % FIGURE NAME
fig_width = 650;	% Window Horizontal Size
fig_heigth = 450;   % Window Vertical Size
fig_posx = 160;     % Window position (Lower left corner)
fig_posy = 160;     % Window position (Lower left corner)

if docked
    figure()
else
    figure('Position', [fig_posx fig_posy fig_width fig_heigth])
end

subplot(n_abbruzzi,1,2)
title('Robots per Thousand Employees', 'interpreter', 'latex')
hold on
plot([2000, 2000+t_bag(t_sucre)], [out.out_initial.RL*1e3; out.RL_path(t_sucre)], 'linewidth', 2, 'color', 'k')
% line([2004, 2004+t_bag(t_sucre(end))], [out.out_initial.RL* 1e3; out.out_initial.RL* 1e3], 'linestyle', '--', 'linewidth', 2, 'color', 'k')
line([2000, 2000+t_bag(t_sucre(end))], [out.out_final.RL* 1e3; out.out_final.RL* 1e3], 'linestyle', '--', 'linewidth', 2, 'color', 'k')
scatter(t_RL, RL_RL, 60, 'filled', 'markerfacecolor', 'k')
hold off
xlim([2000,2060])
ylim([0, out.out_final.RL* 1.1e3])
ylabel('$R/L \times 1000$', 'interpreter', 'latex')

% subplot(n_abbruzzi,1,2)
% title('wage')
% hold on
% plot(t_bag(t_sucre), W_path(t_sucre), 'linewidth', 2)
% line([t_bag(1) t_bag(t_sucre(end))], [out_initial.W, out_initial.W], 'linestyle', '--', 'linewidth', 2, 'color', 'k')
% line([t_bag(1) t_bag(t_sucre(end))], [out_final.W, out_final.W], 'linestyle', '--', 'linewidth', 2, 'color', 'k')
% hold off

subplot(n_abbruzzi,1,1)
title('Relative Robot Price', 'interpreter', 'latex')
hold on
plot([2000, 2000+t_bag(t_sucre)], [transition.par_path(t_sucre); transition.par_path(t_sucre(end))], 'linewidth', 2, 'color', 'k')
% line([t_bag(1) t_bag(end)], [out_initial.W, out_initial.W], 'linestyle', '--', 'linewidth', 2, 'color', 'k')
% line([t_bag(1) t_bag(end)], [out_final.W, out_final.W], 'linestyle', '--', 'linewidth', 2, 'color', 'k')
scatter(t_pR, p_pR, 60, 'filled', 'markerfacecolor', 'k')
ylabel('$p_R /w$', 'interpreter', 'latex')
hold off
xlim([2000,2060])
ylim([0, transition.par_path(t_sucre(1))*1.2])


rescale_L = out.LD_path(14/transition.D);
min_Graph = out.out_final.L /rescale_L- (1/rescale_L-out.out_final.L/rescale_L) * .1;
subplot(n_abbruzzi,1,3)
title('Labor', 'interpreter', 'latex')
hold on 
plot([2000, 2000+t_bag(t_sucre)], [out.out_initial.L; out.LD_path(t_sucre)]./rescale_L, 'linewidth', 2, 'color', 'k')
line([2000, 2000+t_bag(t_sucre(end))], [out.out_final.L; out.out_final.L]./rescale_L, 'linestyle', '--', 'linewidth', 2, 'color', 'k')
line([2000,2014], [1,1], 'linestyle', ':', 'linewidth', 1, 'color', 'k')
line([2014,2014], [min_Graph, 1], 'linestyle', ':', 'linewidth', 1, 'color', 'k')
hold off
ylabel('$L$', 'interpreter', 'latex')
xlabel('Year','interpreter', 'latex')
xlim([2000,2060])
ylim([min_Graph, Inf])
xticks([2000, 2010, 2014, 2020, 2030, 2040, 2050, 2060])
xticklabels({'2000', '2010', '2014', '2020', '2030', '2040', '2050', '2060'})
% subplot(4,1,4)
% title('MC')
% hold on
% plot(t_bag, MC , 'linewidth', 2)
% hold off

if save_figs
    cd([current_path save_path_figures])
    fname = strcat(fig_name);
    set(gcf,'PaperPositionMode','auto');
    hgsave(fname);
    if strcmp(fig_format, 'eps')
        print('-depsc',fname,'-painters');
    elseif strcmp(fig_format, 'jpeg')
        print('-djpeg',fname,'-r256','-opengl');
    elseif strcmp(fig_format, 'png')
        print('-dpng',fname);
    else
        error('Format specified incorrectly.')
    end
    if ~persistent_mode
        close
    end
    cd(current_path)
end

%% CLEAR (SOME) VARIABLES IN MEMORY AND CONTINUE %%
clearvars -EXCEPT save_figs fig_format save_path_tables save_path_figures current_path docked persistent_mode

%% 8) *** Non-Stationary Shocks Estimation ***
%
% This code calibrates the stochastic process parameters for the GBM
% process with exponential death (Kumamoto and Kamihigashi, 2016), with
% stationary LogNormal Double-Pareto Distribution (Reed and Jorgensen,2004)

% Choose whether reimmission follows lognormal or is deterministic
settings.ResetType = 'LogNormal';
% settings.ResetType = 'Deterministic';
% Choose age cutoff to compute exit rate (might want to focus on larger)
age_cutoff = false;
year_cutoff = 10;


% Step 1: Import Raw Residuals obtained from compustat for all 13 sectors
set(0,'DefaultFigureWindowStyle','normal');

opts = delimitedTextImportOptions("NumVariables", 6);
opts.DataLines = [2, Inf];
opts.Delimiter = ",";
opts.VariableNames = ["ifrString", ...
                      "gvkey", ...
                      "year", ...
                      "ifrCode", ...
                      "log_tfp", ...
                      "N_firm"];
opts.VariableTypes = ["string", ...
                      "double", ...
                      "double", ...
                      "double", ...
                      "double",...
                      "double"];
opts = setvaropts(opts, 1, "WhitespaceRule", "preserve");
opts = setvaropts(opts, 1, "EmptyFieldRule", "auto");
opts.ExtraColumnsRule = "ignore";
opts.EmptyLineRule = "read";
Statistics = readtable([current_path '/data/RawResiduals.csv'], opts);

ifrCodesList = unique(Statistics.ifrCode);
ifrString = unique(Statistics.ifrString);

% Set up statistics table to store estimates
StatisticsGBM = table(ifrString, ifrCodesList, 'VariableNames',...
    {'ifrString', 'ifrCode'});
StatisticsGBM.mu = NaN(13,1);
StatisticsGBM.sigma = NaN(13,1);
StatisticsGBM.mu_reset = NaN(13,1);
StatisticsGBM.sigma_reset = NaN(13,1);
StatisticsGBM.reset_rate = NaN(13,1);


TitlesTab = readtable([current_path '/data/OUStatistics']);
TitlesTab = TitlesTab.ifrString;
% Calibration LOOP FOR 13 sectors
for iSec = 1:length(ifrCodesList)
    code = ifrCodesList(iSec);
   
    % Step 2: Assume that the stationary distribution is LogN to obtain
    % guesses of the reset process
    % The reset rate is Method of Moments from the average age of a firm
    % mu and sigma of the reset distribution are estimated from log
    % residual being normal
    StatisticsSec = Statistics(Statistics.ifrCode == code,:);
    StatisticsSec(StatisticsSec.N_firm == 1, :)=[];
    
    gvkey_list = unique(StatisticsSec.gvkey);
    n_gvkey = length(gvkey_list);
    age_gvkey = NaN(length(gvkey_list),1);
    StatisticsSec.diff_log_tfp = NaN(height(StatisticsSec),1);
    for jj = 1:n_gvkey 
        kk = gvkey_list(jj);
        StatisticsSec.diff_log_tfp(StatisticsSec.gvkey == kk) = [NaN; 
            diff(StatisticsSec.log_tfp(StatisticsSec.gvkey == kk))];
        age_gvkey(jj) = StatisticsSec.N_firm(find(StatisticsSec.gvkey == kk,1));
    end
    % Restrict to long  lifetimes

    diff_log = StatisticsSec.diff_log_tfp(StatisticsSec.N_firm>10);
    diff_log(isnan(diff_log)) = [];
    
    % obtain reimmission distribution guesses

    mu_est = sum(diff_log)/length(diff_log);
    sigma_est = sqrt( sum( (diff_log - mu_est) .^ 2) / length(diff_log) );
%     sigma_est = sigma_est;
    if age_cutoff
        reset_rate_est = 1/(mean(age_gvkey(age_gvkey >= year_cutoff))); 
        reset_rate_est = 1/(mean(age_gvkey)); 
    else
        reset_rate_est = 1/(mean(age_gvkey)); 
    end
    


    %% Step 3: calibrate 
    if age_cutoff
        StatisticsSec = StatisticsSec(StatisticsSec.N_firm >= year_cutoff, :);
    else
    end
    % Save a normalized tfp series to fit the distribution
    StatisticsSec.tfp = exp(StatisticsSec.log_tfp)./...
        mean(exp(StatisticsSec.log_tfp));

    % grid settings for GBM
    settings.N_p = 150;
    settings.pmin = 0.001;
    settings.GammaPGrid = 1;
    settings.pmax = min(quantile(StatisticsSec.tfp, .999), 6);
    
    xi = settings.pmin + (settings.pmax - settings.pmin) * ...
        (linspace(0, 1, settings.N_p)).^settings.GammaPGrid;

    % compute kernel densitiy (xi is the increment size and corresponds to
    % grid used by discretize diffusion)
%     xi = settings.pmin:(settings.pmax - settings.pmin)/ ...
%         (settings.N_p - 1):settings.pmax;
    f = ksdensity(StatisticsSec.tfp, xi, 'support', 'positive');


    %general params
    mu_guess = .001;
    sigma_guess = .001;
    % deterministic reset guess
    reset_val_guess =  .99;
    % stochastic reset guess
    sigma_reset_guess = .13;
    mu_reset_guess = log(reset_val_guess) - sigma_reset_guess^2 /2;
    
    % collect all guess values
    settings.mu = @(x) mu_guess .* x;
    settings.sigma = @(x) sigma_guess * x;
    settings.sigma_reset =  sigma_reset_guess;
    settings.mu_reset = mu_reset_guess;
    settings.reset_rate = reset_rate_est;
    % only for deterministic: value where process is reset to
    settings.reset_val = reset_val_guess;
    
    % save the guess distribution
    out = DiscretizeGBM(settings);

    % set up loss function: minimization of quadratic of stationary
    % distribution against the kernel density f estimated above
    handle_loss = @(par) computeLoss([par(1), exp(par(2)), par(3), exp(par(4))],...
        settings, xi, f);
    % solve
    opts = optimset('display', 'iter','MaxIter', 1000);
    guess = [mu_guess log(sigma_guess) mu_reset_guess log(sigma_reset_guess)];
    sol = fminsearch(handle_loss, guess, opts);
%     sol = fmincon(handle_loss, guess,[],[],[],[],...
%         [-Inf log(1e-4) -Inf -Inf],[Inf Inf Inf Inf],[],opts);
    par_sol = [sol(1), exp(sol(2)), sol(3), exp(sol(4))];
    % get parameters
    settings.mu = @(x) par_sol(1) .* x;
    settings.sigma = @(x) par_sol(2) * x;
    settings.mu_reset = par_sol(3);
    settings.sigma_reset =  par_sol(4);
    % compute distribution
    out_fin = DiscretizeGBM(settings);
    % save parameters
    StatisticsGBM.mu(iSec) =  par_sol(1);
    StatisticsGBM.sigma(iSec) = par_sol(2);
    StatisticsGBM.mu_reset(iSec) = par_sol(3);
    StatisticsGBM.sigma_reset(iSec) =  par_sol(4);
    StatisticsGBM.reset_rate(iSec) = reset_rate_est;
    StatisticsGBM.pmax(iSec) =  settings.pmax;

    % FIGURE SETUP: FITTED DISTRIBUTIONS
    fig_name = strcat('B_GBM_dist_', StatisticsGBM.ifrString(iSec));  % FIGURE NAME
    fig_width = 360;	% Window Horizontal Size
    fig_heigth = 270;   % Window Vertical Size
    fig_posx = 160;     % Window position (Lower left corner)
    fig_posy = 160;     % Window position (Lower left corner)
    figure('Position', [fig_posx fig_posy fig_width fig_heigth])
    title(cellstr(strrep(TitlesTab(iSec),'_',' ')),...
        'fontsize',24,'interpreter','latex')
    hold on
    %plot(out.p_grid, out.g_0, 'LineStyle', '--')
    plot(xi, f./sum(f), 'LineWidth', 2, 'color', 'magenta')
    plot(out_fin.p_grid, out_fin.g_0, 'LineWidth', 2, 'color', 'k')
    legend('Actual','Calibrated','fontsize',24,'interpreter','latex','location','best')
    ax = gca;
    ax.FontSize = 24; 
    xlim([0 settings.pmax])
    hold off
    
    if save_figs
        cd([current_path save_path_figures])
        fname = strcat(fig_name);
        set(gcf,'PaperPositionMode','auto');
        hgsave(fname);
        if strcmp(fig_format, 'eps')
            print('-depsc',fname,'-painters');
        elseif strcmp(fig_format, 'jpeg')
            print('-djpeg',fname,'-r256','-opengl');
        elseif strcmp(fig_format, 'png')
            print('-dpng',fname);
        else
            error('Format specified incorrectly.')
        end
        if ~persistent_mode
            close
        end
        cd(current_path)
    end
    
end
writetable(StatisticsGBM, [current_path '/data/GBMStatistics.csv'])

%% CLEAR (SOME) VARIABLES IN MEMORY AND CONTINUE %%
clearvars -EXCEPT save_figs fig_format save_path_tables save_path_figures current_path docked persistent_mode

%% 9) *** Multi-Sector Model with Non-Stationary Shocks Calibration  ***

% For convenience, we provide the file GBM_cs_workspace.mat which contains the 
% pre-computed comparative statics. The reader can use RunGBMCalibration.m to
% produce these simulations.

% RunGBMCalibration

load GBM_cs_workspace

% SHOW CALIBRATION AND BUILD CALIBRATION TABLE

current_path = cd;

which_base_idx = 1;

ifrStrings = cell(13,1);
for jj = 1:13
    ifrStrings{jj} = params{which_base_idx}{jj}.targets.ifrString;
end

% Set up the table and assign row names
calibrationtable = table(ones(13,1));
calibrationtable.Properties.RowNames = ifrStrings;
calibrationtable.Var1 = [];

StatisticsGBM = readtable([current_path '/data/GBMStatistics']);
% Include targets
for jj = 1:13
    calibrationtable.Gamma(jj) = params{which_base_idx}{jj}.Gamma;
    calibrationtable.mu(jj) = StatisticsGBM.mu(jj);
    calibrationtable.sigma(jj) = StatisticsGBM.sigma(jj);
    calibrationtable.mu_reset(jj) = StatisticsGBM.mu_reset(jj);
    calibrationtable.sigma_reset(jj) = StatisticsGBM.sigma_reset(jj);
    calibrationtable.destr_rate(jj) = StatisticsGBM.reset_rate(jj);
end

tab_name = 'calibrationtable_GBM.tex';
tab_preamble = '\\documentclass[../sections/appendix.tex]{subfiles}';
tab_FORMAT = '%4.4f';
table2latex(calibrationtable, [current_path save_path_tables tab_name], tab_FORMAT, tab_preamble);


disp('------------------- CALIBRATION FOR THE MODEL ---------------------')
disp(' ')
disp('Parameters:')
disp(' ')
disp(['$\varphi = ' num2str(params{which_base_idx}{1}.varphi,'%4.3g') '$'])
disp(['$\chi = ' num2str(params{which_base_idx}{1}.chi,'%4.3g') '$'])
disp(['$\A_F = ' num2str(params{which_base_idx}{1}.AF,'%4.3g') '$'])
disp(['$m = ' num2str(params{which_base_idx}{1}.E,'%4.3g') '$'])
disp(['$\psi_R$= ' num2str(params{which_base_idx}{1}.psi_k,'%4g') '$'])
disp(' ')
disp('-------------------------------------------------------------------')
disp(' ')
disp(calibrationtable)
disp(' ')

% PLOT OF CONTINUOUS CS

which_base_idx = 1;

range_cs_mrts = (which_base_idx - 1) * N_cs * 4  + (1:N_cs);
range_cs_mrts_singlesector = (which_base_idx - 1) * N_cs * 4  + (N_cs+1:2*N_cs);
range_cs_adjcost = (which_base_idx - 1) * N_cs * 4  + (2*N_cs+1:3*N_cs);
range_cs_pR_rel = (which_base_idx - 1) * N_cs * 4  + (3*N_cs+1:4*N_cs);

% FIGURE SETUP: LABOR
fig_name = '8a_cs_labor_manysectors_GBM';  % FIGURE NAME
fig_width = 600;	% Window Horizontal Size
fig_heigth = 350;   % Window Vertical Size
fig_posx = 160;     % Window position (Lower left corner)
fig_posy = 160;     % Window position (Lower left corner)

if docked
    figure()
else
    figure('Position', [fig_posx fig_posy fig_width fig_heigth])
end

subplot(2,2,1)
hold on
plot(mrts_ratio_vec * 100,L_cs(range_cs_mrts) * 100,'linewidth',2,'marker','.','color','k')
% plot(mrts_ratio_vec * 100,L_cs_PE(range_cs_mrts) * 100,'linewidth',2,'marker','.','color','r')
set(gca, 'XDir','reverse')
hold off
title('MRTS$_{LR}$', 'interpreter', 'latex')
xlabel('Parameter (\% of 2014)', 'interpreter', 'latex')
ylabel('Labor (\% of 2014)', 'interpreter', 'latex')
grid on

subplot(2,2,2)
hold on
plot(mrts_ratio_singlesector_vec * 100,L_cs(range_cs_mrts_singlesector) * 100,'linewidth',2,'marker','.','color','k')
% plot(mrts_ratio_singlesector_vec * 100,L_cs_PE(range_cs_mrts_singlesector) * 100,'linewidth',2,'marker','.','color','r')
set(gca, 'XDir','reverse')
hold off
title('MRTS$_{LR}$ (Automotive Only)', 'interpreter', 'latex')
xlabel('Parameter (\% of 2014)', 'interpreter', 'latex')
ylabel('Labor (\% of 2014)', 'interpreter', 'latex')
grid on

subplot(2,2,3)
hold on
plot(adj_cost_ratio_vec * 100,L_cs(range_cs_adjcost) * 100,'linewidth',2,'marker','.','color','k')
% plot(adj_cost_ratio_vec * 100,L_cs_PE(range_cs_adjcost) * 100,'linewidth',2,'marker','.','color','r')
set(gca, 'XDir','reverse')
hold off
title('$\psi_R$', 'interpreter', 'latex')
xlabel('Parameter (\% of 2014)', 'interpreter', 'latex')
ylabel('Labor (\% of 2014)', 'interpreter', 'latex')
grid on
% legend({'General Equilibrium', 'Partial Equilibrium'},...
%        'interpreter', 'latex',...
%        'orientation', 'vertical',...
%        'location','southwest',...
%        'fontsize', 11); 

subplot(2,2,4)
hold on
plot(pR_rel_future_vec * 100,L_cs(range_cs_pR_rel) * 100,'linewidth',2,'marker','.','color','k')
% plot(pR_rel_future_vec * 100,L_cs_PE(range_cs_pR_rel) * 100,'linewidth',2,'marker','.','color','r')
set(gca, 'XDir','reverse')
hold off
title('$p_R$', 'interpreter', 'latex')
xlabel('Parameter (\% of 2014)', 'interpreter', 'latex')
ylabel('Labor (\% of 2014)', 'interpreter', 'latex')
grid on
 

if save_figs
    cd([current_path save_path_figures])
    fname = strcat(fig_name);
    set(gcf,'PaperPositionMode','auto');
    hgsave(fname);
    if strcmp(fig_format, 'eps')
        print('-depsc',fname,'-painters');
    elseif strcmp(fig_format, 'jpeg')
        print('-djpeg',fname,'-r256','-opengl');
    elseif strcmp(fig_format, 'png')
        print('-dpng',fname);
    else
        error('Format specified incorrectly.')
    end
    if ~persistent_mode
        close
    end
    cd(current_path)
end

% FIGURE SETUP: LABOR SHARE
fig_name = '8b_cs_laborshare_manysectors_GBM';  % FIGURE NAME
fig_width = 600;	% Window Horizontal Size
fig_heigth = 350;   % Window Vertical Size
fig_posx = 160;     % Window position (Lower left corner)
fig_posy = 160;     % Window position (Lower left corner)

if docked
    figure()
else
    figure('Position', [fig_posx fig_posy fig_width fig_heigth])
end

subplot(2,2,1)
hold on
plot(mrts_ratio_vec * 100,LS_cs(range_cs_mrts)./LS_cs(range_cs_mrts(1)) * 100,'linewidth',2,'marker','.','color','k')
% plot(mrts_ratio_vec * 100,LS_cs_PE(range_cs_mrts)./LS_cs_PE(range_cs_mrts(1)) * 100,'linewidth',2,'marker','.','color','r')
set(gca, 'XDir','reverse')
hold off
title('MRTS$_{LR}$', 'interpreter', 'latex')
xlabel('Parameter (\% of 2014)', 'interpreter', 'latex')
ylabel('Labor Share (\% of 2014)', 'interpreter', 'latex')
grid on

subplot(2,2,2)
hold on
plot(mrts_ratio_singlesector_vec * 100,LS_cs(range_cs_mrts_singlesector)./LS_cs(range_cs_mrts_singlesector(1)) * 100,'linewidth',2,'marker','.','color','k')
% plot(mrts_ratio_singlesector_vec * 100,LS_cs_PE(range_cs_mrts_singlesector)./LS_cs_PE(range_cs_mrts_singlesector(1)) * 100,'linewidth',2,'marker','.','color','r')
set(gca, 'XDir','reverse')
hold off
title('MRTS$_{LR}$ (Automotive Only)', 'interpreter', 'latex')
xlabel('Parameter (\% of 2014)', 'interpreter', 'latex')
ylabel('Labor Share (\% of 2014)', 'interpreter', 'latex')
grid on

subplot(2,2,3)
hold on
plot(adj_cost_ratio_vec * 100,LS_cs(range_cs_adjcost)./LS_cs(range_cs_adjcost(1)) * 100,'linewidth',2,'marker','.','color','k')
% plot(adj_cost_ratio_vec * 100,LS_cs_PE(range_cs_adjcost)./LS_cs_PE(range_cs_adjcost(1)) * 100,'linewidth',2,'marker','.','color','r')
set(gca, 'XDir','reverse')
hold off
title('$\psi_R$', 'interpreter', 'latex')
xlabel('Parameter (\% of 2014)', 'interpreter', 'latex')
ylabel('Labor Share (\% of 2014)', 'interpreter', 'latex')
grid on
% legend({'General Equilibrium', 'Partial Equilibrium'},...
%        'interpreter', 'latex',...
%        'orientation', 'vertical',...
%        'location','southwest',...
%        'fontsize', 11);

subplot(2,2,4)
hold on
plot(pR_rel_future_vec * 100,LS_cs(range_cs_pR_rel) ./LS_cs(range_cs_pR_rel(1)) * 100,'linewidth',2,'marker','.','color','k')
% plot(pR_rel_future_vec * 100,LS_cs_PE(range_cs_pR_rel)./LS_cs_PE(range_cs_pR_rel(1)) * 100,'linewidth',2,'marker','.','color','r')
set(gca, 'XDir','reverse')
hold off
title('$p_R$', 'interpreter', 'latex')
xlabel('Parameter (\% of 2014)', 'interpreter', 'latex')
ylabel('Labor Share (\% of 2014)', 'interpreter', 'latex')
grid on


if save_figs
    cd([current_path save_path_figures])
    fname = strcat(fig_name);
    set(gcf,'PaperPositionMode','auto');
    hgsave(fname);
    if strcmp(fig_format, 'eps')
        print('-depsc',fname,'-painters');
    elseif strcmp(fig_format, 'jpeg')
        print('-djpeg',fname,'-r256','-opengl');
    elseif strcmp(fig_format, 'png')
        print('-dpng',fname);
    else
        error('Format specified incorrectly.')
    end
    if ~persistent_mode
        close
    end
    cd(current_path)
end

%% CLEAR (SOME) VARIABLES IN MEMORY AND CONTINUE %%
clearvars -EXCEPT save_figs fig_format save_path_tables save_path_figures current_path docked persistent_mode

%% 10) ***** Run Illustrative Case for Operating Profit Function *****

% Set parameters
params = SetParameters();
 
% Wage and Price of Robots
W = params.W;
pR = params.pR;
 
% Compute Labor Savings and R_max
Omega = ( 1 - params.Gamma) /params.Gamma * params.W - params.E;

[~, out] = LaborDemand_trapz(W, pR, params);

% PLOT OF OPERATING PROFIT AND MARGINAL OPERATING PROFIT
% Operating Profit and the Marginal Operating Profit
R_grid_figure = linspace(0, 5, 1000);
p_grid_figure = linspace(0, 2.5, 1000);

% CHOOSE values of R and z at which to slice the function
R_choice_figure = 1.5;
z_choice_figure = 1.1;

% FIGURE SETUP
fig_name = 'A_ProfitFunctionSlices';  % FIGURE NAME
fig_width = 650;	% Window Horizontal Size
fig_heigth = 400;   % Window Vertical Size
fig_posx = 160;     % Window position (Lower left corner)
fig_posy = 160;     % Window position (Lower left corner)

if docked
    figure()
else
    figure('Position', [fig_posx fig_posy fig_width fig_heigth])
end

% Operating Profit as function of R
subplot(2,2,1)
plot(R_grid_figure, out.handles.Pi(R_grid_figure,z_choice_figure),...
    'linewidth', 3, 'color', 'k')
hold on
line([out.handles.R_bar(z_choice_figure), out.handles.R_bar(z_choice_figure)],...
     [0, out.handles.Pi(out.handles.R_bar(z_choice_figure),z_choice_figure)],...
     'linestyle',':',...
     'linewidth', 2,...
     'color', 'k')
line([out.handles.R_hat(z_choice_figure), out.handles.R_hat(z_choice_figure)],...
     [0, out.handles.Pi(out.handles.R_hat(z_choice_figure),z_choice_figure)],...
     'linestyle',':', ...
     'linewidth', 2,...
     'color', 'k')
hold off
xticks([out.handles.R_bar(z_choice_figure), out.handles.R_hat(z_choice_figure)])
xticklabels({['$\bar{R}(' num2str(z_choice_figure) ')$'],...
             ['$\hat{R}(' num2str(z_choice_figure) ')$']})
yticks([])
xlim([0 max(R_grid_figure)])
ylim([0 1.1 * out.handles.Pi(out.handles.R_hat(z_choice_figure),z_choice_figure)])
a = get(gca,'XTickLabel');
set(gca,'XTickLabel',a,'fontsize',14)
ylabel(['$\Pi(R,' num2str(z_choice_figure) ')$'], 'fontsize', 14, 'interpreter', 'latex')
xlabel('Robot Stock $R$', 'fontsize', 14, 'interpreter', 'latex')
title({'Operating Profit';...
       ['$z = ' num2str(z_choice_figure) '$'] },...
       'FontWeight','Normal',...
      'fontsize',16,...
      'interpreter','latex')


% Operating Profit as function of p
subplot(2,2,2)
plot(p_grid_figure, out.handles.Pi(R_choice_figure,p_grid_figure),...
    'linewidth', 3, 'color', 'k')
hold on
line([out.handles.p_hat(R_choice_figure), out.handles.p_hat(R_choice_figure)],...
     [0, out.handles.Pi(R_choice_figure,out.handles.p_hat(R_choice_figure))],...
     'linestyle',':',...
     'linewidth', 2,...
     'color', 'k')
line([out.handles.p_bar(R_choice_figure), out.handles.p_bar(R_choice_figure)],...
     [0, out.handles.Pi(R_choice_figure,out.handles.p_bar(R_choice_figure))],...
     'linestyle',':',...
     'linewidth', 2,...
     'color', 'k')
hold off
xticks([out.handles.p_hat(R_choice_figure), out.handles.p_bar(R_choice_figure)])
xticklabels({['$\hat{z}(' num2str(R_choice_figure) ')$'],...
             ['$\bar{z}(' num2str(R_choice_figure) ')$']})
yticks([])
xlim([0 max(p_grid_figure)])
ylim([0 1.1 * out.handles.Pi(R_choice_figure,max(p_grid_figure))])
a = get(gca,'XTickLabel');
set(gca,'XTickLabel',a,'fontsize',14)
ylabel(['$\Pi(' num2str(R_choice_figure) ',z)$'],...
        'fontsize', 14,...
        'interpreter', 'latex')
xlabel('Revenue-Shifter $z$', 'fontsize', 14, 'interpreter', 'latex')
title({'Operating Profit';...
       ['$R = ' num2str(R_choice_figure) '$'] },...
       'FontWeight','Normal',...
      'fontsize',16,...
      'interpreter','latex')

% Marginal Operating Profit as a function of R
subplot(2,2, 3)
plot(R_grid_figure, out.handles.Pi_prime(R_grid_figure,z_choice_figure),...
    'linewidth', 3, 'color', 'k')
hold on
line([out.handles.R_bar(z_choice_figure), out.handles.R_bar(z_choice_figure)],...
     [0, out.handles.Pi_prime(out.handles.R_bar(z_choice_figure),z_choice_figure)],...
     'linestyle',':',...
     'linewidth', 2,...
     'color', 'k')
hold off
xticks([out.handles.R_bar(z_choice_figure), out.handles.R_hat(z_choice_figure)])
xticklabels({['$\bar{R}(' num2str(z_choice_figure) ')$'],...
             ['$\hat{R}(' num2str(z_choice_figure) ')$']})
yticks([Omega])
yticklabels({'$\Omega$'})
a = get(gca,'XTickLabel');
set(gca,'XTickLabel',a,'fontsize',14)
a = get(gca,'YTickLabel');
set(gca,'YTickLabel',a,'fontsize',14)
ylabel(['$\Pi_R(R,' num2str(z_choice_figure) ')$'], 'fontsize', 14, 'interpreter', 'latex')
xlabel('Robot Stock $R$', 'fontsize', 14, 'interpreter', 'latex')
xlim([0 max(R_grid_figure)])
ylim([0 1.1 * Omega])
title({'Marginal Operating Profit';...
       ['$z = ' num2str(z_choice_figure) '$'] },...
       'FontWeight','Normal',...
      'fontsize',16,...
      'interpreter','latex')

% Marginal Operating Profit as a function of p
subplot(2,2,4)
plot(p_grid_figure, out.handles.Pi_prime(1.5,p_grid_figure),...
    'linewidth', 3, 'color', 'k')
hold on
line([out.handles.p_hat(R_choice_figure), out.handles.p_hat(R_choice_figure)],...
     [0, out.handles.Pi_prime(R_choice_figure,out.handles.p_hat(R_choice_figure))],...
     'linestyle',':',...
     'linewidth', 2,...
     'color', 'k')
line([out.handles.p_bar(R_choice_figure), out.handles.p_bar(R_choice_figure)],...
     [0, out.handles.Pi_prime(R_choice_figure,out.handles.p_bar(R_choice_figure))],...
     'linestyle',':',...
     'linewidth', 2,...
     'color', 'k')
line([0, out.handles.p_bar(R_choice_figure)],...
     [Omega, Omega],...
     'linestyle',':',...
     'linewidth', 2,...
     'color', 'k')
hold off
xticks([out.handles.p_hat(R_choice_figure), out.handles.p_bar(R_choice_figure)])
xticklabels({['$\hat{z}(' num2str(R_choice_figure) ')$'],...
             ['$\bar{z}(' num2str(R_choice_figure) ')$']})
yticks([Omega])
yticklabels({'$\Omega$'})
xlim([0 max(p_grid_figure)])
ylim([0 1.1 * Omega])
a = get(gca,'XTickLabel');
set(gca,'XTickLabel',a,'fontsize',14)
ylabel(['$\Pi_R(' num2str(R_choice_figure) ',z)$'],...
        'fontsize', 14,...
        'interpreter', 'latex')
xlabel('Revenue-Shifter $z$', 'fontsize', 14, 'interpreter', 'latex')
title({'Marginal Operating Profit';...
       ['$R = ' num2str(R_choice_figure) '$'] },...
       'FontWeight','Normal',...
      'fontsize',16,...
      'interpreter','latex')

if save_figs
    cd([current_path save_path_figures])
    fname = strcat(fig_name);
    set(gcf,'PaperPositionMode','auto');
    hgsave(fname);
    if strcmp(fig_format, 'eps')
        print('-depsc',fname,'-painters');
    elseif strcmp(fig_format, 'jpeg')
        print('-djpeg',fname,'-r256','-opengl');
    elseif strcmp(fig_format, 'png')
        print('-dpng',fname);
    else
        error('Format specified incorrectly.')
    end
    if ~persistent_mode
        close
    end
    cd(current_path)
end
 
%% CLEAR (SOME) VARIABLES IN MEMORY AND CONTINUE %%
clearvars -EXCEPT save_figs fig_format save_path_tables save_path_figures current_path docked persistent_mode

%% 11) ***** Empirical Performance of Estimation of EOU Process *****

% Load data
TableResids = readtable([current_path '/data/HamiltonEstimationSample']);
ifrCodesList = unique(TableResids.ifrCode);
ifrString = unique(TableResids.industry_ifr19);
OUStatistics = readtable([current_path '/data/OUStatistics']);


for iSec = 1:length(ifrCodesList)
    code = ifrCodesList(iSec);
   
    StatisticsSec = TableResids(strcmp(cellstr(TableResids.ifrCode),code),:);
    
    params.sigma_P = OUStatistics.sigma_p(iSec);
    params.theta_P = OUStatistics.theta_p(iSec);
    params.P = 1;

    params.sigma_hat = sqrt(log(params.sigma_P.^2 ./ params.P.^2 + 1));
    params.mu_hat = log(params.P) - 1/2 .* log(params.sigma_P.^2 ./ params.P.^2 + 1);
    % definition of drift and s.d.
    settings.mu = @(p) - params.theta_P .* (log(p) - params.mu_hat - (params.sigma_hat).^2) .* p;
    settings.sigma = @(p) sqrt(2 * params.theta_P) .* params.sigma_hat .* p;
    % Exogenous shock
    opts = optimoptions('fsolve','Display','off', 'UseParallel', false);
    settings.pmin = fsolve(@(p) logncdf(p, params.mu_hat, params.sigma_hat) - 1e-6, 1,opts);
    settings.pmax = fsolve(@(p) 1-logncdf(p, params.mu_hat, params.sigma_hat) - 1e-6, 1,opts);
    
    settings.N_p = 150;
    settings.GammaPGrid = 1;
    xi = settings.pmin + (settings.pmax - settings.pmin) * ...
        (linspace(0, 1, settings.N_p)).^settings.GammaPGrid;

    % compute kernel densitiy (xi is the increment size and corresponds to
    % grid used by discretize diffusion)
%     xi = settings.pmin:(settings.pmax - settings.pmin)/ ...
%         (settings.N_p - 1):settings.pmax;
    f = ksdensity(exp(StatisticsSec.p), xi, 'support', 'positive');
    
    out = DiscretizeDiffusion(settings);
        % FIGURE SETUP: FITTED DISTRIBUTIONS
    fig_name = strcat('B_EOU_dist_', string(code));  % FIGURE NAME
    fig_width = 360;	% Window Horizontal Size
    fig_heigth = 270;   % Window Vertical Size
    fig_posx = 160;     % Window position (Lower left corner)
    fig_posy = 160;     % Window position (Lower left corner)
    figure('Position', [fig_posx fig_posy fig_width fig_heigth])

    title(cellstr(strrep(OUStatistics.ifrString(iSec),'_',' ')),...
        'fontsize',24,'interpreter','latex')
    hold on
    plot(out.p_grid, out.g_0, 'LineWidth', 2, 'color', 'k')
    plot(xi, f./sum(f), 'LineWidth', 2, 'color', 'magenta')
    legend('Calibrated', 'Actual','fontsize',24,'interpreter','latex','location','best')
    ax = gca;
    ax.FontSize = 24; 
    xlim([0 settings.pmax])
    hold off

    if save_figs
        cd([current_path save_path_figures])
        fname = strcat(fig_name);
        set(gcf,'PaperPositionMode','auto');
        hgsave(fname);
        if strcmp(fig_format, 'eps')
            print('-depsc',fname,'-painters');
        elseif strcmp(fig_format, 'jpeg')
            print('-djpeg',fname,'-r256','-opengl');
        elseif strcmp(fig_format, 'png')
            print('-dpng',fname);
        else
            error('Format specified incorrectly.')
        end
        if ~persistent_mode
            close
        end
        cd(current_path)
    end
    
end
