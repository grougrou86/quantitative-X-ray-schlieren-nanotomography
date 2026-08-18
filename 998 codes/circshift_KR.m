function output = circshift_KR(input, vec, varargin)
    p = inputParser();
    addParameter(p, 'gridCell', {}, @iscell); %% 
    addParameter(p, 'fourierInputTF', false, @islogical);

    parse(p, varargin{:});
    gridCell         = p.Results.gridCell;
    fourierInputTF    = p.Results.fourierInputTF;


    %% parsing    
    mustBeVector(vec);
    if iscolumn(vec)
        vec = vec.';
    end

    if ~isempty(gridCell)
        assert( all( cellfun(@(a) all(size(a) == size(input,1:ndims(a))) ,gridCell) ), ...
            'size of gridCell is not compatible with input' )
    end

    %%
    % input = mk_ellipse(5,5,16,16);
    % vec = [0.1*-1,0];
    sizeIn = size(input);
    ShiftingDimsN =  length(vec);
    
    %%% integer shifts done
    if ~fourierInputTF
        input = circshift(input, vec .* (round(vec) == vec) );
        nonIntInd = find(round(vec) ~= vec);
    else
        nonIntInd = find(vec);
    end
    
    %%% non integer shifts
    if isempty(nonIntInd)
        output = input;
    else        
        %%% phaseMat calc
        phaseMat = zeros(sizeIn,'like',input(1)+1i);
        if ~isempty(gridCell)
            for dd_dim = nonIntInd
                phaseMat = phaseMat + vec(dd_dim).* gridCell{dd_dim};
            end
        else
            phaseMat = zeros(sizeIn,'like',input(1)+1i);
            for dd_dim = nonIntInd
                dd_size = sizeIn(dd_dim);


                %%% shifted array order
                mdd = floor(dd_size/2) + 1;
                rampdd = [0:1: (dd_size-mdd), (1-mdd):1:-1 ]; % ifftshift(rampdd);
                rampdd = rampdd.' / dd_size;

                reshapeDim = [dd_size, ones(1,ShiftingDimsN-1)];
                reshapeDim = circshift(reshapeDim, dd_dim - 1 );
                rampdd   = reshape( rampdd, reshapeDim );

                phaseMat = phaseMat + vec( dd_dim ).*rampdd;
            end
        end
        phaseMat = exp(-1i*2*pi*phaseMat);

        %%% fourier space multiplication
        if fourierInputTF
            output = input .* phaseMat;
        else
            realTF = isreal(input);
            for dd_dim = nonIntInd
                input = fft(input,[],dd_dim);
            end            
            
            output = input .* phaseMat;
            
            for dd_dim = nonIntInd
                output = ifft(output,[],dd_dim);
            end
            
            if realTF
                output = real(output);
            end
        end
    end
    
