function wl = Etowl(E)
%   E is keV unit.
%   wl is nm unit.
c = 299792458e9;  % nm/s 
h = 4.13566766225e-18; % keV s

wl = h*c./E;