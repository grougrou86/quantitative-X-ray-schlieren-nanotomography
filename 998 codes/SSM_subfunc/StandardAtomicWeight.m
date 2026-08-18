function out = StandardAtomicWeight(symbol,varargin)
% https://en.wikipedia.org/wiki/Standard_atomic_weight 
% and the references therein

p = inputParser();
addParameter(p,'topDir','',@ischar);  % Top path 

parse(p,varargin{:});

topDir = p.Results.topDir;
if isempty(topDir)
    topDir = mfilename('fullpath');
end

%%
load(fullfile(topDir,'atomicWeightList.mat'),'wOut','symbolOut');

testFun = @(x) strcmpi(x,symbol);
ind = find(cellfun(testFun, symbolOut), 1);

if isempty(ind)
    error('Invalid atomic symbol')
end

out = wOut(ind);

if isnan(out)
    error('Standard atomic weight of %s is unknown', symbol)
end
