function [th, s,CLimOut] = contrastVis(inputImg,varargin)
p = inputParser;
addOptional(p,'inlierPortion',99.5,@(x) (x>=0)&&(x<=100));
addParameter(p,'refFrame',[]);

%e.g., inlierPortion = 98% means,
%               1   :      98      :  1
%   (lower outliner):(viewing data):(upper outliner)

parse(p,varargin{:});
inlierPortion = p.Results.inlierPortion;
refFrame = p.Results.refFrame;

inputImg = squeeze(inputImg);

%%

%%% random sampling
Nsample = 2^16;
if numel(inputImg) > Nsample
    randSampInd = randperm(numel(inputImg), Nsample);
    testVec = inputImg(randSampInd);
    % if isempty(refFrame)
    %     refFrame = floor( size(inputImg,3)/2 ) + 1; % select the frame middle
    % end
    % testVec = inputImg(:,:,refFrame);

else
    testVec = inputImg;
end

%%%
loPct = (100-inlierPortion)/2;
hiPct = loPct+inlierPortion;

goodInd = isfinite(testVec);
N = numel(goodInd);
testVec = testVec(goodInd);

if N > 1
    % sortVals = sort(testVec);
    % loCval = round((N-1)*loPct/100 + 1);
    % hiCval = round((N-1)*hiPct/100 + 1);
    % CLimOut  = [sortVals(loCval),sortVals(hiCval)];
    CLimOut  = [prctile(testVec, loPct), prctile(testVec, hiPct)];
    CLimOut = cast(CLimOut,'double');
    
    if CLimOut(1) == CLimOut(2)
        CLimOut(2) = CLimOut(1) + eps('double');
    end

else
    CLimOut = [NaN,NaN];
end


th = gcf;
if  ndims(inputImg) == 3
    %     s = sliceViewer(inputImg);    
    s = orthosliceViewer( gather(inputImg) );
    s.DisplayRangeInteraction = 'off';
    if ~any(isnan(CLimOut))
        s.DisplayRange = gather( CLimOut );
    end
    th.ToolBar = 'figure';
    th.WindowState = 'maximized';
else
    s=imagesc(inputImg);
    ax = gca;
    if ~any(isnan(CLimOut))
        ax.CLim = gather( CLimOut );
    end
end



