% Fornino Manera (2019)
% 
% DATE: January 24, 2021
%
% Project: Automation and the Future of Work: Assessing the Role of Labor
%          Flexibility

function out = DiscretizeDiffusion(settings)
% This function discretizes the stationary Ito process:
% dx = mu(x) dt + sigma(x) dW
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

AT = spdiags([F_p' , - F_p' - B_p', B_p' ], [-1 0 1], N_p, N_p);

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
    error('problem with DiscretizeDiffusion')
end

out.B_p = B_p;
out.F_p = F_p;
out.g_0 = g_0;
out.Lambda = AT';
out.dp_mean = dp_mean;





