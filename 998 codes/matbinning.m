function mat = matbinning( mat, binSize, symmetricity )   
%     mat = ones(9,4,5);
%     binSize = 4;
    if nargin < 3
        symmetricity = 'asymmetric';
    end
    
    %% binSize test
    mustBeVector(binSize)
    if iscolumn(binSize)
        binSize = binSize.';
    end
    
    matDim = ndims(mat);
    binSizeL = length(binSize);
    
    if matDim > binSizeL
        if isscalar(binSize)
            binSize = binSize .* ones(1, matDim);
        else
            binSize = [binSize, ones(1, matDim-binSizeL ) ];
        end
        
    end
    
    %% padding
    matSize = size(mat);
    assert( all(matSize >= binSize), 'matrix size must be largert than binning size' )
    
    if strcmp(symmetricity,'symmetric')
        assert(all( mod(binSize,2) == 1), 'binning must be odd to maintain the symmetricity')
        D2 = matSize - mcoor(matSize) + (1 - binSize)/2; % number of pixels but the central bin
        
        finalSize = 1 + 2* floor(D2./binSize); % number of bins =center + 2*wings (always odd)
        finalSize(binSize == 1) = matSize(binSize == 1);
    else
        finalSize = floor( matSize ./ binSize );
    end
    mat       = mcrop(mat, finalSize.* binSize);
    
    
    %% binnding
    nonUnitybinningDims = find(binSize > 1);
    
    for dd = nonUnitybinningDims        
        % make dd the first dim
        dimOrder = circshift(1:matDim, 1-dd );
        mat = permute(mat, dimOrder );
        
        % binning the dimension
        tempSize = size(mat);
        mat = reshape(mat, binSize(dd) ,[]);
        mat = sum(mat,1);
        mat = reshape(mat, [finalSize(dd) ,tempSize(2:end)] );
                
        % ipermute
        mat = ipermute(mat, dimOrder );
    end
    
%     cellDist = cell(1,3);
%     for dd = 1:matDim
%         cellDist{dd}= repmat(binSize(dd), [1, finalSize(dd)] );
%     end
%     
%     mat = mat2cell(mat,cellDist{:});
%     mat = cellfun(@(x) sum(x,'all'),mat);
    
    
    
    
    
