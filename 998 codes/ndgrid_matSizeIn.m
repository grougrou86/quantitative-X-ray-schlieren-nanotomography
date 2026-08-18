function outGridCell = ndgrid_matSizeIn(sizeMat, normalizeTF ,meshType)
%%% put sizeMat gpuArray for gpuArray gen.
%%% sizeMat = gpuArray( [5,6,7,8] );
%%% mustBeMember(meshType, {'default','centerZero','centerZero_ifftshift'})  


%% assert
mustBeVector(sizeMat);    

if nargin < 2
    normalizeTF = 0;
else
    if ischar(normalizeTF)
        mustBeMember(normalizeTF, {'default','normalized'})
        if strcmp(normalizeTF,'default')
            normalizeTF = 0;
        else
            normalizeTF = 1;
        end
    end
end

if nargin < 3
    meshType = 'default';
else
    mustBeMember(meshType, {'default','centerZero','centerZero_ifftshift'})  
end

%% run

dimN = length(sizeMat);
linearVecCell = cell(1,dimN);
outGridCell   = cell(1,dimN);

switch meshType
    case 'default' % plain ndgrid
        for dd = 1:dimN
           linearVecCell{dd} = 1:sizeMat(dd);
        end
    case 'centerZero' % make center (floor(dd_size/2) + 1) zero
        for dd = 1:dimN
            dd_size = sizeMat(dd);
            mdd_size = floor(dd_size/2) + 1;
            linearVecCell{dd} = (1:sizeMat(dd)) - mdd_size;
        end
    case 'centerZero_ifftshift' % make center (floor(dd_size/2) + 1) zero, and ifftshift
        for dd = 1:dimN
           dd_size = sizeMat(dd);
           mdd_size = floor(dd_size/2) + 1;
           linearVecCell{dd} = [0:1: (dd_size-mdd_size), (1-mdd_size):1:-1 ]; % ifftshift(rampdd);
        end
end

%%% normalization
if normalizeTF
    for dd = 1:dimN
       linearVecCell{dd} = linearVecCell{dd} ./ sizeMat(dd);
    end
end

[ outGridCell{:} ] = ndgrid( linearVecCell{:} );
end