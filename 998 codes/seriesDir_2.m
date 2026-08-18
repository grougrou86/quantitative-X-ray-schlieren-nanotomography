function [EList, commonExp] = seriesDir_2(readFile, readPath)


cd(readPath)

% list .mat files
filesList = dir('*');

% find numbers
[startInd, endInd] = regexp(readFile, '\d+');

% find commonExp
Nfiles = 0;
for ii = 1:length(startInd)
    testExp = [readFile(1:startInd(ii)-1), '\d+', readFile(endInd(ii)+1:end)];

    fileExpMatch    = regexp({filesList.name}, testExp); % index of all numbers
    fileExpMatchInd = ~cellfun('isempty',fileExpMatch);

    if sum(fileExpMatchInd) >= Nfiles
        commonExp = testExp;
        EList     = filesList(fileExpMatchInd);
        Nfiles    = length(EList);
    end
end
fprintf('*** %d file series detected ***\n',Nfiles)
end