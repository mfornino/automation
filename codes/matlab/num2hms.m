% Fornino Manera (2019)
% 
% DATE: January 24, 2021
%
% Project: Automation and the Future of Work: Assessing the Role of Labor
%          Flexibility

function hms = num2hms(seconds)
% This subroutine converts seconds into a string "#h, #m, #s" 

h = floor(seconds / 3600);
m = floor((seconds - h * 3600)/60);
s = seconds - h * 3600 - m * 60;

hms = [num2str(h) 'h ' num2str(m) 'm ' num2str(s), 's'];

end