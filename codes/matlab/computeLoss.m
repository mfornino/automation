% Fornino Manera (2019)
% 
% DATE: January 18, 2021
%
% Project: Automation and the Future of Work: Assessing the Role of Labor
%          Flexibility

function loss = computeLoss(params_vec, settings, xi, f)
% Computes loss function associated to the distribution-fitting routine that
% calibrates the parameters of the non-stationary GBM case.

mu_est = params_vec(1);
sigma_est = params_vec(2);
mu_reset_guess = params_vec(3);
sigma_reset_guess = params_vec(4);

settings.mu = @(x) mu_est .* x;
settings.sigma = @(x) sigma_est .* x;
settings.sigma_reset =  sigma_reset_guess;
settings.mu_reset = mu_reset_guess;

out = DiscretizeGBM(settings);

% ensure that points are on grid with interpolation
f_p = interp1(xi,f,out.p_grid);
normalized_f = f_p./sum(f_p);
KL = out.g_0(normalized_f>0)' ./(normalized_f(normalized_f>0)) - 1;
KL_loss = KL;
weights = normalized_f(normalized_f>0)./(sum(normalized_f(normalized_f>0)));
loss = sum((KL_loss .* weights).^2) + ...
    (1 - out.p_grid * out.g_0 ).^2 ;

    
end