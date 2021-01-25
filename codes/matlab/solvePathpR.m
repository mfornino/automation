% Fornino Manera (2019)
% 
% DATE: January 24, 2021
%
% Project: Automation and the Future of Work: Assessing the Role of Labor
%          Flexibility

function [alpha, pR_rel_infty] = solvePathpR(varargin)
% Calibrate an exponential path for pR to match observed data points.

narginchk(0,1)

base_year = 2010;
base_year_idx = 4;
t = [1995, 1996, 2005, 2010, 2014, 2017]' - base_year;

% Note 1995 and 1996 wages are set to 1997 level because we have no data
% for those two years.
w = [25492, 25492, 29890, 31073, 35447, 38070];

p = [133382, 121452, 68659, 46150, 31776, 27074]' ./ w';
p_base = p(base_year_idx);

OPTS = optimset('display', 'iter');

if nargin == 0

    f = @(x) loss_function(p, t, exp(x(1)), x(2), p_base);
    x0 = [log(.01), 1];
    x = fminsearch(f, x0, OPTS);
    pR_rel_infty = exp(x(1));
    alpha = x(2);
    
else
    
    pR_rel_infty = varargin{1};
    f = @(x) loss_function(p, t, pR_rel_infty, x(1), p_base);
    x0 = 1;
    x = fminsearch(f, x0, OPTS);
    alpha = x(1);

end
%     
% 
% 
% figure()
% hold on
% times = linspace(t(1), 40, 200);
% plot(times + base_year, path(times, pR_rel_infty, alpha, p_base), 'linewidth', 2, 'color', 'k', 'linestyle', '--')
% scatter(t + base_year, p, 100, 'filled', 'linewidth', 2, 'MarkerFaceColor', 'k')
% legend('Fitted Path', 'Data', 'interpreter', 'latex')
% xlabel('Year', 'interpreter', 'latex')
% ylabel('Relative Robot Price $p_R / w$', 'interpreter', 'latex')
% ylim([0, inf])
% xlim([t(1) + base_year, inf])
% % title('Path of Relative Robot Prices')
% grid on
% hold off

end

function loss = loss_function(p, t, pR_infty, alpha, p_base)

    loss = sum((path(t, pR_infty, alpha, p_base) - p).^2);

end

function p_t = path(t, pR_infty, alpha, p_base)

    p_t = pR_infty + (p_base - pR_infty) .* exp(- alpha .* t);
end