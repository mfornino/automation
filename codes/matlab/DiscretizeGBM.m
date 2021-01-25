% Fornino Manera (2019)
% 
% DATE: January 24, 2021
%
% Project: Automation and the Future of Work: Assessing the Role of Labor
%          Flexibility

function out = DiscretizeGBM(settings)
% This function discretizes the stationary Ito process:
%
% dx = mu x dt + sigma x dW
% 
% for the case with an exponential reset rate.
%
% Input ("settings" struct object)
% - mu and sigma: function handles
% - N_p: number of grid points
% - pmin and pmax: extremes of the grid
% - GridType: string either "LinSpace" or "LogSpace"
% - GammaPGrid: scaling parameter for linspace option (see below line 45)
%
% Outputs("out" struct object)
% - F_p B_p: Forward and Backward drifts
% - g_0: stationary distribution
% - Lambda: Infinitesimal Generator
% - dp_mean: increments on the (unequally spaced) grid

out = settings;
mu = settings.mu;
sigma = settings.sigma;
N_p = settings.N_p;
pmin = settings.pmin;
pmax = settings.pmax;
reset_rate = settings.reset_rate;

if isfield(settings, 'ResetType')
    ResetType = settings.ResetType;
else
    ResetType = 'Deterministic';
end

if strcmp(ResetType, 'Deterministic')
    reset_val = settings.reset_val;
elseif strcmp(ResetType, 'LogNormal')
    mu_reset = settings.mu_reset;
    sigma_reset = settings.sigma_reset;
else
    error('Supported ResetType: Deterministic OR LogNormal');
end



if isfield(settings, 'GridType')
    GridType = settings.GridType;
else
    GridType = 'LinSpace';
end

if isfield(settings, 'GammaPGrid')
    GammaPGrid = settings.GammaPGrid;
else
    GammaPGrid = 1; %defaults to simple linspace
end
    
if ~strcmp(GridType, 'LogSpace')
    p_grid = pmin + (pmax - pmin)* (linspace(0, 1, N_p)).^GammaPGrid;
else
    p_grid = logspace(log10(pmin), log10(pmax), N_p);
end

if isfield(settings, 'p_grid')
    p_grid = settings.p_grid;
end

out.p_grid = p_grid;


% Define the grid for dp
dp_vec_up = diff(p_grid);
dp_vec_up = [dp_vec_up dp_vec_up(end)];
dp_vec_dwn = diff(p_grid) ;
dp_vec_dwn = [dp_vec_dwn(1) dp_vec_dwn] ;
dp_mean = (dp_vec_up + dp_vec_dwn)/2;
dp_mean(1) = dp_mean(1)/2;
dp_mean(end) = dp_mean(end)/2;

mu_vec = mu(p_grid);
sigma_vec = sigma(p_grid);

B_p = - min(mu_vec, 0)./(dp_vec_dwn) + (sigma_vec.^2)/2 ./((dp_vec_dwn + dp_vec_up)/2 .*dp_vec_dwn );
F_p =   max(mu_vec, 0)./(dp_vec_up) + (sigma_vec.^2)/2 ./((dp_vec_dwn + dp_vec_up)/2 .*dp_vec_up);

B_p(1) = 0;
F_p(end) = 0;

% Transition matrix conditional on surviving
AT_survive = spdiags([F_p' , - F_p' - B_p', B_p' ], [-1 0 1], N_p, N_p);
[row, col, val] = find(AT_survive);
if strcmp(ResetType, 'Deterministic') 
    % find reset index
    p_reset_idx = find(p_grid > reset_val, 1);
    col = [col; (1:N_p)'; (1:N_p)'];
    row = [row; (1:N_p)'; repmat(p_reset_idx, N_p, 1)];
    val = [val; -repmat(reset_rate, N_p, 1); repmat(reset_rate, N_p, 1)];
elseif strcmp(ResetType, 'LogNormal')
    logNormalPdf = lognpdf(p_grid, mu_reset, sigma_reset);
    reset_rates_logn = logNormalPdf./(sum(logNormalPdf));
    in_indexes = find(reset_rates_logn>1e-6);
    reset_rates_logn = reset_rates_logn(in_indexes)./...
        (sum(reset_rates_logn(in_indexes)));
    col = [col; (1:N_p)'; kron((1:N_p)', ones(length(in_indexes), 1))];
    row = [row; (1:N_p)';  kron(ones(N_p,1) , in_indexes' ) ];
    val = [val; -repmat(reset_rate, N_p, 1); ...
    kron(ones(N_p, 1), reset_rate * reset_rates_logn') ];
    out.in_indexes = in_indexes;
    out.reset_rates_logn = reset_rates_logn;
end

AT = sparse(row, col, val, N_p, N_p);
%disp(['transition sums to ' num2str(sum(sum(AT,1)))])
if abs(sum(sum(AT,1))) > 1e-8
    warning(['transition sums to ' num2str(sum(sum(AT,1)))])
end

DT = 1e8;
II = speye(size(AT));
TransMat = plus(II, - DT*AT);

g_0 = ones(N_p,1)/N_p;

gg = ones(size(g_0));

supnorm = 1;
it = 0;
maxit = 500;
while supnorm > 1e-8 && it < maxit
    gg = TransMat \  g_0;
    supnorm = max(abs(gg-g_0));
    g_0 = gg;
    it = it + 1;
end

if supnorm > 1e-6
    warning(['Maxit for statDist Reached. Supnorm: ' num2str(supnorm)])
end

out.B_p = B_p;
out.F_p = F_p;
out.g_0 = g_0;
% AT' to use with value function
out.Lambda = AT_survive';
% AT' to use to obtain statdist
out.AT_KFE = AT';
out.dp_mean = dp_mean;




