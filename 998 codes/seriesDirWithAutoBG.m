function [samList,bgList,outStruct] = seriesDirWithAutoBG(varargin)
[samList, samExp] = seriesDir(varargin{:});
[~,~,ext] = fileparts(samExp);

%%% bgPattern def.
allFiles = dir(['*',ext]);
[~, nsInds] = setdiff( {allFiles.name} , {samList.name});

nsFiles = allFiles(nsInds);
nsN = length(nsFiles);
bgPatternLib = cell(nsN,1);
for bb = 1:nsN
    commonStringIndex = regexp(nsFiles(bb).name, '\d+');
    bgPatternLib{bb} = nsFiles(bb).name(1:commonStringIndex(end)-1);
end
bgPattern = unique(bgPatternLib);

if isempty(bgPattern)    
    error('ERROR: no bgPattern has been detected.')

elseif length(bgPattern) > 1
    fprintf('*** More than one bgPattern has been detected ***\n');
    [bgList, bgExp]   = seriesDir(['*',ext], 'select bg');  
else
    bgExp = [bgPattern{1},'*',ext];
    bgList = dir( bgExp );     
end

outStruct = struct('samExp',samExp,'bgExp',bgExp,'ext',ext);
end