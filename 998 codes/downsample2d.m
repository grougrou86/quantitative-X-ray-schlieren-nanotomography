function output= downsample2d(input,n,phase,varargin)
% default --> phase = 0; centeredTF == 0
p = inputParser();
addParameter(p, 'centeredTF', 0); 

parse(p, varargin{:});
centeredTF = p.Results.centeredTF;

%%
if nargin <= 2
    phase = 0;
end

if length(n) == 1
    n = [n,n];
elseif length(n) >=3 
    error('length(n) must be 1 or 2');
end

if centeredTF
    phase = mod(floor(size(input)/2),n) + phase;
end

output = downsample(downsample(input,n(1),phase(1)).',n(2),phase(2)).';

% subplot(122),imagesc(testds),axis image
