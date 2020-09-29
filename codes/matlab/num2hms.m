function hms = num2hms(seconds)

% Fornino Manera (2019)
%
% This subroutine converts seconds into a string "#h, #m, #s" 


h = floor(seconds / 3600);
m = floor((seconds - h * 3600)/60);
s = seconds - h * 3600 - m * 60;

hms = [num2str(h) 'h ' num2str(m) 'm ' num2str(s), 's'];

end