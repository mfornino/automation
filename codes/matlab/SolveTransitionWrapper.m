% Fornino Manera (2019)
% 
% DATE: January 24, 2021
%
% Project: Automation and the Future of Work: Assessing the Role of Labor
%          Flexibility

function eq = SolveTransitionWrapper(x, targets, transition, params)
% Routine that calls SolveTransition.m and then computes the loss function
% associated to the RL path as deviation from the actual data points passed in
% targets.

params.Gamma = x(1);
params.psi_k = x(2);

out = SolveTransition(transition, params);

eq = sum([out.RL_path(4/transition.D) - targets(1), ...
          out.RL_path(7/transition.D) - targets(2), ...
          out.RL_path(10/transition.D) - targets(3), ...
          out.RL_path(14/transition.D) - targets(4)].^2);

end