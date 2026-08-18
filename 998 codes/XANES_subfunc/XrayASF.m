function asfOut = XrayASF( symbol, varargin )
p = inputParser();
addParameter(p,'topDir','',@ischar);  % Top path 
addParameter(p,'show',0);  % Top path

parse(p,varargin{:});

topDir = p.Results.topDir;
showPlot = p.Results.show;
if isempty(topDir)
    topDir = mfilename('fullpath');
end

%%
load(fullfile(topDir,'ASF.mat'),'dataOut','symbolOut');

testFun = @(x) strcmpi(x,symbol);
ind = find(cellfun(testFun, symbolOut), 1);

if isempty(ind)
    error('Invalid atomic symbol')
end

%% -9999 filtering

asfOut = dataOut{ind}; % ev, f1, f2
asfOut(asfOut<-1000) = NaN;

if showPlot
    semilogx(asfOut(:,1),asfOut(:,2),asfOut(:,1),asfOut(:,3))
    titleStr = sprintf('%s, Z = %d', symbol, ind);
    title(titleStr);

    legend('f1','f2')
    xlabel('Energy (eV)')
    set(gcf,'Color','w')
end
