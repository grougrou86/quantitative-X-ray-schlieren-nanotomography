function cRI = XrayComplexRI(Energy,material,density)
% put Energy in keV
% put density in g cm^-3
% put meterial in element symbol

%% E range

if any( Energy > 30 | Energy < 0.03 )
    error('Energies must be in the range 0.03 keV < E < 30 keV')
end

%% contants
r0 = 2.817940322719e-6; % Classical electron radius, nm
AvoN = 6.02214076e23; % Avogadro constant, mol^-1

%% density
if nargin < 3
    density = defaultDensityLoad(material);
elseif density < 0
    density = defaultDensityLoad(material);
end

% density input as g cm^-3
density = density * 1e-21; % g nm^-3;

%%
atomsSplit = regexp(material,'[A-Z][^A-Z]*','match');
Ntype = length(atomsSplit);

atomsSym = regexp(atomsSplit,'[A-Za-z]*','match');
atomsNumCell = regexp(atomsSplit,'[^A-Za-z]*','match');
atomsNum = zeros(size(atomsSplit));

for ii = 1:Ntype
    atomsSym{ii} = atomsSym{ii}{1};

    if isempty(atomsNumCell{ii})
        atomsNum(ii) = 1;
    else
        if length(atomsNumCell{ii}) == 1
            atomsNum(ii) = str2double(atomsNumCell{ii}{1});
        else
            error('invalid material input')
        end
    end
end

%%
ASFsum   = zeros(size(Energy));
ASFcontant = r0/2/pi * Etowl(Energy).^2 * density * AvoN; % g mol^-1
molWeight = 0;

for ii = 1:Ntype
    ASF = XrayASF( atomsSym{ii} );
    
    nanInds = isnan(sum(ASF,2));    
    f1  = interp1(ASF(~nanInds,1)/1000, ASF(~nanInds,2), Energy, 'pchip',NaN);
    f2  = interp1(ASF(~nanInds,1)/1000, ASF(~nanInds,3), Energy, 'pchip',NaN);
    % plot(Energy,f1,'-', ASF(:,1)/1000, ASF(:,2),'--')

    ASFsum     = ASFsum + atomsNum(ii) * (-f1 + 1i*f2);
    molWeight  = molWeight + atomsNum(ii) * StandardAtomicWeight( atomsSym{ii} );  % g mol^-1
    % molWeight = 196.97;
end

cRI = 1 + ASFcontant .* (ASFsum / molWeight); 
% plot(photonE,1-real(cRI),photonE,imag(cRI))


end
function density = defaultDensityLoad(material)
topDir = mfilename('fullpath');
% load(fullfile(topDir,'densityList.mat'),'densityOut','materialOut');
loadCell = readcell(fullfile(topDir,'densityList.txt'));

materialList = loadCell(:,1);
densityList = loadCell(:,2);

testFun = @(x) strcmpi(x,material);
ind = find(cellfun(testFun, materialList), 1);

if isempty(ind)
    fprintf('*** No default density for %s..\n',material)
    density = input('*** Please update it (visit : https://henke.lbl.gov/optical_constants/getdb2.html) : ');

    if isempty(density) || density < 0
        error('invalid density')
    else
        %%% update the file
        appendList = {material, density};
        writecell(appendList, fullfile(topDir,'densityList.txt'),'WriteMode','append','Delimiter','\t');
    end
else
    density = densityList{ind}; % g cm^-3
end

end






%{

p = inputParser();
addParameter(p,'topDir','',@ischar);  % Top path 
parse(p,varargin{:});

topDir = p.Results.topDir;
if isempty(topDir)
    topDir = mfilename('fullpath');
end

fileDir= fullfile(topDir, [material,'_RI_Xray.txt']);

if isfile(fileDir)
    EvsRI = load(fileDir);
else
    error('The element is not on the list: please update the code with RI table from for example, http://henke.lbl.gov/optical_constants')
end


ELib=EvsRI(:,1)/1000; % keV
deltaLib=EvsRI(:,2);
betaLib=EvsRI(:,3);

delta = interp1(ELib,deltaLib,Energy,'pchip');
beta  = interp1(ELib,betaLib,Energy,'pchip');

cRI = (1-delta) + 1i*beta;

%}