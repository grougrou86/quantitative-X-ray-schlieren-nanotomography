classdef tomoHandle < handle
    properties
        % constant
        version = 6.0; 
        % v6.0   = SSM tomo update
        % v5.0   = getTomo bgMask change; make axisOffset, ringRemovePrm to
        %          properties.
        % v4.1   = sinoNorm change (normalization of every [y,z] pair)
        % v4.0   = MajorChanges in ring pattern removal
        % v3.2   = savename change (tomoHandle --> tH) @ save
        % v3.1   = bgnorm roll-back @ getTomogram

        % defined by constructor(tomoHandle)
        setup struct = struct([]);
        fieldReconMethod char {mustBeMember(fieldReconMethod,{'proj','PC','SSM','KK'})} = 'PC';
        tomoReconMethod char {mustBeMember(tomoReconMethod,{'FBP','Gridrec'})} = 'FBP';

        deconvTF logical = false;
        gpuInd = [];

        tomoRes {mustBeNonnegative, mustBeVector(tomoRes,"allow-all-empties")} = [];
        rotAngleMax double = 180;

        % defined by setData
        samFiles struct = struct([]);
        bgFiles  struct = struct([]);
        darkFiles struct = struct([]);
        rotAnglelib {mustBeNonnegative,mustBeVector(rotAnglelib,"allow-all-empties")} = [];

        % data-dependant parameters
        FOV struct = struct([]);
        positionOffs (:,2) = [0,0];
        SNRMap (:,:) {mustBeNumeric} = [];
        axisOffset {mustBeScalarOrEmpty} = [];
        ringRemovePrm struct = struct([]);
    end

    methods(Static)
        function [th,s,CLimOut] = show(b,varargin)
            [th,s,CLimOut] = contrastVis(b,varargin{:});
        end

        function [th,s,CLimOut] = matShow(fileDir,showType,varargin)
            %%% parsing
            if nargin < 1
                fileDir = '';
            end
            if nargin < 2
                showType = 'auto';
            end

            [bout,fileDir] = matParser(fileDir,showType);
            [th,s,CLimOut] = tomoHandle.show(bout, varargin{:});
            %             [~,fileName,~]=fileparts(fileDir);
            set(gcf,'Name',fileDir)

        end

        function outFolderDir = imgFileExport(fileDir,ext,showType)
            %%% parsing
            if nargin < 1
                fileDir = '';
            end

            if nargin < 2
                ext = 'png';
            end

            if isempty(ext)
                ext = 'png';
            end

            if nargin < 3
                showType = 'auto';
            end
            mustBeMember(ext,{'png','tif','tifstack'})

            %%% paser
            [bout,fileDir,outType] = matParser(fileDir,showType);
            mustBeMember(outType,{'sino','delta','beta'}) % check the outType
            % tomoHandle.show(bout); % beta

            %%% mkdir
            topDir = fileparts(fileDir);
            outFolderDir = fullfile(topDir,['imageStack_in_',ext]);
            mkdir(outFolderDir);

            %%% uint16 converstion
            minval = min(bout(:));
            maxval = max(bout(:));

            valPerDigit = (maxval-minval) / (2^16-1);
            bout = bout - min(bout(:));
            bout = uint16( bout / valPerDigit );

            %%% pemutation if required
            if ismember(outType,{'delta','beta'}) % it is a tomogram
                bout = permute(bout,[3,2,1]); % sinogram = (y,x) images, tomogram = (z,x) images
            end

            %%% image files generation
            N = size(bout,3);
            if strcmp(ext,'tifstack')
                imgSaveStr = fullfile(outFolderDir, 'tifStack.tif' );
                for nn = 1:N
                    imwrite(bout(:,:,nn) ,imgSaveStr, 'WriteMode', 'append', 'Compression','none');
                end
            else
                parfor nn = 1:N
                    imgSaveStr = fullfile(outFolderDir, [num2str(nn,'%05d'),'.',ext] );
                    imwrite(bout(:,:,nn) ,imgSaveStr);
                end
            end

            %%% write specs
            fileID = fopen(fullfile(outFolderDir,'details.txt'),'w');
            fprintf(fileID,'outType = %s \n',outType);
            fprintf(fileID,'valPerDigit = %.4g\n', valPerDigit);
            fprintf(fileID,'minval = %.4g\n', minval);

            fprintf(fileID,'\n to get physical values, do: (valPerDigit) * (img) + (minval)');
            fclose(fileID);
        end

        function mkGIF(savestr,tomoIn,varargin)

            p = inputParser();
            addParameter(p, 'scanningDim',1, @isscalar);
            addParameter(p, 'LoopPeriod',1, @isscalar);

            parse(p, varargin{:});
            LoopPeriod      = p.Results.LoopPeriod;
            scanningDim      = p.Results.scanningDim;

            %%
            Nscan = size(tomoIn, scanningDim);
            tomoIn = shiftdim(tomoIn, -scanningDim +1 ); % make the first dimension to the scanning Dim.

            DelayTime = LoopPeriod/Nscan;

            h = figure;
            colormap(turbo);

            for n = 1:Nscan
                imagesc(squeeze(tomoIn(n,:,:)));
                ax = gca;

                if n == 1
                    setCLim = ax.CLim;
                else
                    ax.CLim = setCLim;
                end
                axis image;
                colorbar;
                titleStr = sprintf('%05d / %05d',n,Nscan);
                title(titleStr)
                drawnow

                % Capture the plot as an image
                frame = getframe(h);
                im = frame2im(frame);
                [imind,cm] = rgb2ind(im,256);

                % Write to the GIF File
                if n == 1
                    imwrite(imind,cm,savestr,'gif','Loopcount',inf,'DelayTime',DelayTime);
                else
                    imwrite(imind,cm,savestr,'gif','WriteMode','append','DelayTime',DelayTime);
                end
            end
        end

        %%% KK related
        function  [Sinogram,Fsupport] = KKfieldRecon(Sinogram,NAfrac,OCDfrac,offAngle,bgAddXFOV,gpuInd,figShow,ReconAlgorithm)
            
            %%% assert
            xx0 = size(Sinogram,2);
            assert( mod(xx0,2) == 1, 'horizontal axis must be odd' );
            assert( ismember(offAngle,[0,180]), 'offAngle must be either 0 or 180' );
            assert( ismember(ReconAlgorithm,{'ori','iter'}), 'ReconAlgorithm must be either ''ori'' or ''iter''');

            %%% padding
            padSz = (xx0 - 1)/2;
            bgL = repmat( mean(Sinogram(:,1:bgAddXFOV,:), 2),        [1,padSz,1] );
            bgR = repmat( mean(Sinogram(:,end-bgAddXFOV+1:end,:),2), [1,padSz,1] );
            Sinogram = cat(2, bgL, Sinogram ,bgR);
            % tomoHandle.show(Sinogram);


            %%% Fsupport generation
            [yy,xx,zz] = size(Sinogram);
            vecfunc = @(v) ((1:v) - mcoor(v))/v;

            xvec = vecfunc(xx);
            yvec = vecfunc(yy);
            [XX,YY] = meshgrid(xvec,-yvec);

            yc= sind(offAngle)*OCDfrac;
            xc= cosd(offAngle)*OCDfrac;

            Fsupport = (XX-xc).^2 + (YY-yc).^2 <= NAfrac^2; % NAmask
            Fsupport = (cosd(offAngle)*XX + sind(offAngle)*YY >= 0) & Fsupport;                % KKmask

            %%% GPU?
            if gpuInd
                if isgpuarray(Sinogram)
                    gpu = gpuDevice;
                else
                    gpu = gpuDevice(gpuInd);
                    Sinogram = gpuArray(Sinogram);
                end
                Fsupport = gpuArray(Fsupport);

                %%% chunkSize determination
                bytesPerAngle = (xx*yy)*8; % 8 bytes for complex single
                maxChunkSize = floor( gpu.AvailableMemory / bytesPerAngle / 10 ); % last factor is a heuristic safety value.

                iterN = ceil( zz / maxChunkSize);
                chunks = round( linspace(1, zz+1, iterN+1) );
                for ii = 1:iterN
                    thisChunk = chunks(ii) : chunks(ii+1)-1;
                    if strcmp(ReconAlgorithm,'iter')
                        Sinogram(:,:,thisChunk) = kkSolver_v1( Sinogram(:,:,thisChunk), Fsupport,'figVisTF',figShow);
                    else
                        Sinogram(:,:,thisChunk) = kkSolver_vOri( Sinogram(:,:,thisChunk), offAngle, 2);
                        Sinogram(:,:,thisChunk) = fftshift(fft2(ifftshift(Sinogram(:,:,thisChunk))));
                        Sinogram(:,:,thisChunk) = fftshift(ifft2(ifftshift(Sinogram(:,:,thisChunk).*Fsupport)));
                    end
                    % kkSolver gives global phase = 0 for each reconstructions.

                    fprintf('Field recon. finished: chunk: %d / %d \n', ii,iterN);
                end

            else %%% noGPU
                if strcmp(ReconAlgorithm,'iter')
                    Sinogram = kkSolver_v1(Sinogram, Fsupport);
                else
                    Sinogram = kkSolver_vOri( Sinogram, offAngle,2);
                    Sinogram = fftshift(fft2(ifftshift(Sinogram)));
                    Sinogram = fftshift(ifft2(ifftshift(Sinogram.*Fsupport)));
                end
                fprintf('Field recon. finished\n');
            end

            %%% crop
            Sinogram = mcrop(Sinogram,[yy,xx0]);
            % tomoHandle.show(angle(Sinogram));

        end

        function bField = rotAxisShift(bField,shift,bgAddXFOV)
            % shift = 0;
            % bField = bFieldSelectedY;
            % bgAddXFOV = obj.FOV.bgAddXFOV;

            if shift == 0
                return;
            end

            if shift > 0 % going right
                padSize = ceil( shift );
                padVal  = mean(bField(:,1:bgAddXFOV,:),2); % left buffer
                buffer = repmat(padVal,[1,padSize,1]);

                bField = cat(2, buffer, bField); % put left
                bField = circshift_KR( bField, [0, shift ,0] );
                bField(:,1:padSize,:) = [];  % equalize the size

            elseif shift < 0 % going left
                padSize = ceil( -shift );
                padVal  = mean(bField(:,end-bgAddXFOV+1:end,:),2); % right buffer
                buffer = repmat(padVal,[1,padSize,1]);

                bField = cat(2,bField,buffer); % put right
                bField = circshift_KR( bField, [0, shift ,0] );
                bField(:,end-padSize+1:end,:) = [];   % equalize the size
            end
        end

        function tomoOut = tomoRotate(tomoIn,Rang) % rotate along y-axis
            tomoIn = permute(tomoIn, [3,2,1]);
            tomoOut = zeros(size(tomoIn),'like',tomoIn(1));

            for yi = 1:size(tomoIn,3)
                tomoOut(:,:,yi) = imrotate(tomoIn(:,:,yi), Rang, 'bilinear','crop');
            end
            tomoOut = ipermute(tomoOut, [3,2,1]);
        end

        %%% SSM related

        function [Sinogram,correctionDone] = yaxisRegistration(Sinogram)
            mfft1  = @(x) fftshift(fft(ifftshift( x ,1),[],1),1);
            mifft1 = @(x) fftshift(ifft(ifftshift( x ,1),[],1),1);

            realTF   = isreal(Sinogram);
            zeroMask  = all( abs(Sinogram)~=0 ,3);

            %%% test_beforeCorrection
            test_beforeCorrection   = squeeze( sum(Sinogram,2) );

            %%% y FOV set
            h = figure(151543);
            h.WindowState= 'maximized';
            ax = gca;
            imagesc(real(test_beforeCorrection)); title('set FOV for the registration')
            Q = input('======= ENTER after FOV adjustment ==== (put N for exit) \n','s');
            if strcmpi(Q,'N')
                correctionDone = 'N';
                return;
            end
            ylimTest   = [ceil(ax(1).YLim(1)), floor(ax(1).YLim(2))];

            test0    = test_beforeCorrection(ylimTest(1):ylimTest(2),:);
            zz = size(test0,2);
            % figure, imagesc(abs(test0)),axis image

            yshiftResult = zeros(1,zz);
            test = test0;
            for iter = 1:10
                Ftest     = mfft1(test);

                %%% origin define
                origin    = mean(test,2);
                Forigin   = mfft1(origin);
                % figure, plot(abs(origin))
                %     Forigin   = exp(1i*angle(Forigin));
                % imagesc(log10(abs(Ftest))),axis image

                Fccoor = mpad( Ftest.*conj(Forigin), [size(Ftest,1)*10,size(Ftest,2)]);
                % plot(Fccoor(:,115)); ylim([-pi,pi])
                % imagesc(angle(Fccoor)),axis image
                % imagesc(abs(Fccoor))

                ccoor = mifft1(Fccoor);
                % imagesc(abs(ccoor)),axis image

                %% pixelwise shift
                [~,maxind] = max(abs(ccoor));
                yshift_Iter = (maxind - mcoor(size(ccoor,1))) /10;
                % plot(yshiftInd)

                %%% break for very slight changes
                if  sqrt( mean( abs(yshift_Iter).^2 ) ) < 0.1
                    break;
                end

                yshiftResult = yshiftResult + yshift_Iter;

                test = test0;
                for zi = 1:zz
                    test(:,zi) = circshift_KR( test0(:,zi), [ -yshiftResult(zi),0] );
                end

                %%% crop
                min_yshiftInd = round( min(-yshiftResult) ); % <= 0
                max_yshiftInd = round( max(-yshiftResult) ); % >= 0
                test = test( max_yshiftInd+1 : end+min_yshiftInd,:);
                %     pause;

                %%% vis?
                %                 figure(151543);
                %                 subplot(121),imagesc(abs(test0)),axis image; title('before Correction')
                %                 subplot(122),imagesc(abs(test)),axis image;   title(['after Correction, Iter: ', num2str(iter)])
                %                 drawnow;
            end

            %%% comparison
            test_afterCorrection   = test_beforeCorrection;
            for zi = 1:zz
                test_afterCorrection(:,zi) = ...
                    circshift_KR( test_afterCorrection(:,zi), [ -yshiftResult(zi),0] );

            end
            min_yshiftInd = round( min(-yshiftResult) ); % <= 0
            max_yshiftInd = round( max(-yshiftResult) ); % >= 0

            test_afterCorrection = test_afterCorrection( max_yshiftInd+1 : end+min_yshiftInd,:);

            figure(151543);
            subplot(121),imagesc(abs(test_beforeCorrection)); title('before Correction')
            subplot(122),imagesc(abs(test_afterCorrection));  title(['Done, after Correction, Iter: ', num2str(iter)])


            %% manual check
            correctionDone = input('*** Correction works? Y/N/R [Y]:','s');
            if any( strcmpi(correctionDone, {'N','R'}) ) % not qualified
                % do nothing.
            else
                correctionDone = 'Y';
                % y-axis registration application.

                Sinogram = fft2(Sinogram);
                for zi = 1:1:zz
                    Sinogram(:,:,zi)  = circshift_KR( Sinogram(:,:,zi), [-yshiftResult(zi),0,0], 'fourierInputTF', true );
                end

                if realTF
                    Sinogram = ifft2(Sinogram,'symmetric');
                else
                    Sinogram = ifft2(Sinogram);
                end

                min_yshiftInd = round( min(-yshiftResult) ); % <= 0
                max_yshiftInd = round( max(-yshiftResult) ); % >= 0

                zeroMask = ...
                    circshift(zeroMask,[max_yshiftInd,0]) &...
                    circshift(zeroMask,[min_yshiftInd,0]);                
                Sinogram = Sinogram .* zeroMask;
                Sinogram = Sinogram( max_yshiftInd+1 : end+min_yshiftInd,:,:);
                %                 obj.FOV.ylim = obj.FOV.ylim + gather( [max_yshiftInd, min_yshiftInd] );
            end
            close(151543)

        end

        function bField = removeStripeOfSino(bField,rmvMethod,varargin)
            %%% input bField is log(intensity) or phase
                     
            if isreal(bField)
                bField = runFunc( bField );
            else
                bField = runFunc(real(bField)) + 1i*runFunc(imag(bField)); % intensity or exp(phase)
            end

            function runningB = runFunc(runningB)
                runningB = exp(runningB);  % intensity or exp(phase)

                [yy,xx,zz] = size(runningB); 
                if strcmpi(rmvMethod,'wavelet-fft')
                    waveletN =varargin{1};
                    wname    = varargin{2};
                    sigma    = varargin{3};
                    parfor yi = 1:yy
                        stripe_rmved_slice = removeStripes('h', squeeze(runningB(yi,:,:)),waveletN,wname,sigma);
                        runningB(yi,:,:) = reshape(stripe_rmved_slice,1,xx,zz);
                    end

                else % default
                    snr        = varargin{1};
                    la_size    = varargin{2};
                    sm_size    = varargin{3};

                    %%% this function currently working only in CPU
                    if isgpuarray(runningB)
                        input_was_gpuArray = true;
                        runningB = gather(runningB);
                    else
                        input_was_gpuArray = false;
                    end

                    %%%
                    parfor yi = 1:yy
                        stripe_rmved_slice = vo_removeStripes('h', squeeze(runningB(yi,:,:)), snr, la_size, sm_size);
                        runningB(yi,:,:) = reshape(stripe_rmved_slice,1,xx,zz);
                    end

                    if input_was_gpuArray
                        runningB = gpuArray(runningB);
                    end
                end

                runningB = log(runningB);  % intensity or exp(phase)
            end
            %%% bg nomarlization of each measurements
            %             bField = obj.sinoNorm(bField);
        end

        
    end

    methods

        %%% GPU is not used
        function [obj,jsonPath] = tomoHandle(jsonPath,fieldReconMethod_in)
            if nargin == 0
                jsonPath = [];
            end
            jsonPath = loadSetupJson(obj,jsonPath);
            cd(fileparts(jsonPath))

            %%% set fieldReconMethod
            if nargin == 2
                obj.fieldReconMethod = fieldReconMethod_in;
            end

            %%% fieldReconMethod default setting
            switch obj.fieldReconMethod
                case {'proj', 'PC'}
                    obj.rotAngleMax = 180;
                    obj.deconvTF = true;
                case 'KK'
                    obj.rotAngleMax = 360;
                    obj.deconvTF = true;
                case 'SSM'
                    obj.rotAngleMax = 180;
                    obj.deconvTF = false;
            end

            %%% tomo resolution calc.
            switch obj.fieldReconMethod
                case {'proj', 'PC'}
                    ZPdrn  = obj.setup.afterSample.zonePlate.outermostWidth;
                    obj.tomoRes = [1, 1] * ZPdrn;

                case 'KK'
                    ZPdia  = obj.setup.afterSample.zonePlate.diameter;
                    ZPdrn  = obj.setup.afterSample.zonePlate.outermostWidth;
                    OCD    = obj.setup.afterSample.zonePlate.offAxis.distance;      % nm

                    obj.tomoRes = [ 1, (1 + 2*OCD/ZPdia)^-1 ] * ZPdrn;
                case 'SSM'
                    XrayWl  = Etowl( (obj.setup.source.energy_eV)*1e-3 );
                    diffDia = obj.setup.afterSample.diffuser.usedDiameter;
                    L1      = obj.setup.afterSample.diffuser.distance2sample;

                    %%% this resolution does not consider the effect of CZP
                    obj.tomoRes = L1 .* XrayWl ./ diffDia * [1, 1];

            end

            %%% gpu setting.
            if gpuDeviceCount > 0
                gpuTable = gpuDeviceTable(["Index","ComputeCapability",...
                    "TotalMemory","MultiprocessorCount","DeviceAvailable"]);

                availGPUs = find(gpuTable.DeviceAvailable);
                [~,bestGPUInd] = max( gpuTable.TotalMemory( availGPUs ) );
                obj.gpuInd = availGPUs( bestGPUInd );
            end
        end

        function jsonPath = loadSetupJson(obj,jsonPath)
            if nargin == 1 || isempty(jsonPath)
                [readFile, readPath] = uigetfile('*.json*','Select a setup .json file');
                jsonPath = fullfile(readPath,readFile);
            end
            obj.setup = jsondecode(fileread(jsonPath));
            obj.setup.jsonPath = jsonPath;
        end

        function [samExp,bgExp,ext] = setData(obj,readFile, readPath)
            folderBefore = pwd;

            %%% seriesDirWithAutoBG
            if(nargin<2)
                [obj.samFiles, obj.bgFiles, outStruct] = seriesDirWithAutoBG('*','select sample');
            else
                [obj.samFiles, obj.bgFiles, outStruct] = seriesDirWithAutoBG_2(readFile, readPath);
            end
            samExp = outStruct.samExp;
            bgExp  = outStruct.bgExp;
            ext    = outStruct.ext;

            %%% rotation angle def.
            Nangle = length(obj.samFiles);
            dRotAngles  = obj.rotAngleMax / (Nangle - 1);
            obj.rotAnglelib = 0 : dRotAngles  : obj.rotAngleMax;

            if ~ismember(180, obj.rotAnglelib)
                warning('180-degree is not included in rotation angles')
            end

            fprintf('\tsamPattern:\t%10s\t\t%05d\n',samExp,length(obj.samFiles))
            fprintf('\t bgPattern:\t%10s\t\t%05d\n',bgExp, length(obj.bgFiles))
            fprintf('\t rotAngles:\t%d : %.2f : %d\n',0,dRotAngles, obj.rotAngleMax)

            cd(folderBefore)

            %%% reset the data-dependant parameters
            obj.FOV = struct([]);
            obj.positionOffs = 0;
            obj.SNRMap  = [];
            obj.axisOffset = [];
            obj.ringRemovePrm = struct([]);
        end

        function [FOV, b] = setFOV(obj,percent_cropp,b,varargin)
            p = inputParser();
            addParameter(p, 'bgAddXFOVPercent', 10); % 10% for default. put Inf for max bg area

            parse(p, varargin{:});
            bgAddXFOVPercent  = p.Results.bgAddXFOVPercent;

            %%% representative indices            
            repInds = round( linspace(1,length(obj.rotAnglelib),6) );
            repAngle = obj.rotAnglelib( repInds );

            %%% reset
            obj.FOV = struct([]); % reset FOV;
            obj.positionOffs = 0; % reset positionOffs;

            %%% load
            if nargin < 3
                %%% no sinogram input
                b0 = [];
                b = obj.getSinogram(repInds);                
            else
                %%% sinogram input
                b0 = b;
                b = b(:,:,repInds);

                %%% if the input sinogram is complex
                if ~isreal(b)
                    b = angle(b);
                end
            end

            %%% FOV.originalSize def
            FOV.originalSize = size(b,[1,2]);

            %%% plot
            use_FIG_FOV=false;
            if use_FIG_FOV
                ax(1:6) = matlab.graphics.axis.Axes;
                h = figure('Name','Adjust the FOV to make the sample fit to the boundary','WindowState','maximized');
                for aa = 1:6
                    ax(aa) = subplot(2,3,aa);
                    imagesc(b(:,:,aa)); axis square;
                    title(repAngle(aa))
                end
                colormap(turbo);
                linkaxes(ax,'xy');
                % h.WindowState= 'maximized';

                input('======= ENTER after FOV adjustment ==== \n');

                %%% FOV.ylim def
                FOV.ylim = [ceil(ax(1).YLim(1)), floor(ax(1).YLim(2))];

                %%% samXlim def
                samXlim  = [ceil(ax(1).XLim(1)), floor(ax(1).XLim(2))];
            else
                FOV.ylim = [1+floor(size(b,1)*((1-percent_cropp)/2)),size(b,1)-floor(size(b,1)*((1-percent_cropp)/2))];
                samXlim  = [1+floor(size(b,2)*((1-percent_cropp)/2)),size(b,2)-floor(size(b,2)*((1-percent_cropp)/2))];
            end
            samXFOV = samXlim(2)-samXlim(1)+1;
            % make samXFOV odd
            if mod(samXFOV,2) == 0
                samXFOV = samXFOV - 1;
                samXlim(1) = samXlim(1) + 1;
            end

            %%% FOV.bgAddXFOV def
            bgAddXFOV_given = floor( bgAddXFOVPercent/200 *samXFOV ); % set bgAddXFOV by bgAddXFOVPercent
            bgAddXFOV_L     = samXlim(1) - 1;                   % left pixels on the left
            bgAddXFOV_R     = FOV.originalSize(2) - samXlim(2); % left pixels on the right
            FOV.bgAddXFOV = min([bgAddXFOV_given, bgAddXFOV_L, bgAddXFOV_R]);

            %%% FOV.xlim def
            FOV.xlim  = samXlim + [-FOV.bgAddXFOV, +FOV.bgAddXFOV]; % this is also odd

            %%% FOV.size def
            FOV.size = [FOV.ylim(2)-FOV.ylim(1), FOV.xlim(2)-FOV.xlim(1)] + 1;

            disp(FOV);
            if use_FIG_FOV
                close(h)
            end
            obj.FOV = FOV;

            %%% apply FOV if sinogram input exists
            if (~isempty(b0)) && (nargout > 1)
                b = b0(FOV.ylim(1):FOV.ylim(2), FOV.xlim(1):FOV.xlim(2),:);
            end

            % 2023.01.09 update: considering bgAddXFOVPercent = Inf
            %{
                    %%% maxSampleFOV
                    FOV.originalSize = size(b,[1,2]);
                    oneSidebgAddXFOVRatio = bgAddXFOVPercent / (bgAddXFOVPercent + 100) / 2;
                    maxSampleFOV = FOV.originalSize(2) - 2*floor(oneSidebgAddXFOVRatio * FOV.originalSize(2) ); % we must take 10% for bg.
                    viwewingXind = (1:maxSampleFOV) - mcoor(maxSampleFOV) + mcoor(FOV.originalSize(2));
        
        
                    %%% plot
                    ax(1:6) = matlab.graphics.axis.Axes;
                    h = figure(151543);
                    h.WindowState= 'maximized';
                    for aa = 1:6
                        ax(aa) = subplot(2,3,aa);
                        imagesc(b(:,viwewingXind,aa)); axis square;
                        title(repAngle(aa))
                    end
                    colormap(turbo);
                    linkaxes(ax,'xy');
        
                    %%% FOV check
                    input('======= ENTER after FOV adjustment ==== \n');
        
                    FOV.ylim    = [ceil(ax(1).YLim(1)), floor(ax(1).YLim(2))];
                    samXlim = [ceil(ax(1).XLim(1)), floor(ax(1).XLim(2))] + viwewingXind(1) - 1;
                    samXFOV = samXlim(2)-samXlim(1)+1;
        
                    % make samXFOV odd
                    if mod(samXFOV,2) == 0
                        samXFOV = samXFOV - 1;
                        samXlim(1) = samXlim(1) + 1;
                    end
        
                    FOV.bgAddXFOV = floor( bgAddXFOVPercent/200 *samXFOV );
                    FOV.xlim  = samXlim + [-FOV.bgAddXFOV, +FOV.bgAddXFOV]; % this is also odd
                    FOV.size = [FOV.ylim(2)-FOV.ylim(1), FOV.xlim(2)-FOV.xlim(1)] + 1;
        
                    disp(FOV);
                    close(151543)
                    obj.FOV = FOV;
            %}
        end

        function [sam,bg_field] = getSinogram(obj,varargin)
            Nbg   = length(obj.bgFiles);
            Nsam  = length(obj.samFiles);
            samFilesLocal = obj.samFiles;
            bgFilesLocal  = obj.bgFiles;
            darkDir       =   fullfile( obj.darkFiles.folder,  obj.darkFiles.name);
          
            p = inputParser();
            addOptional(p, 'shotNumber', [], @(v) isvector(v) || isempty(v));
            addOptional(p, 'bgSubMethod', 'average', @(v) ismember(v,{'average','alternative'}) );
            addParameter(p, 'AUTO_offset', false, @islogical);
            addParameter(p, 'figShow', false);

            parse(p, varargin{:});
            shotNumber = p.Results.shotNumber;
            bgSubMethod = p.Results.bgSubMethod;
            AUTO_offset = p.Results.AUTO_offset;
            figShow            = p.Results.figShow;
            %%% rotAngle def.
            if isfield(obj.setup.detector.camera,'tiltAngle')
                rotAngle = -(obj.setup.detector.camera.tiltAngle); % compensate camera tiltAngle
            else
                rotAngle = 0;
            end

            %%% shotNumber def
            if isempty(shotNumber)
                shotNumber = 1:Nsam;
                Nshot = Nsam;
            else
                Nshot = length(shotNumber);
                shotNumber = sort(shotNumber);
                assert(all( shotNumber>0 & shotNumber<=Nsam ),'Invalid shotNumber');
            end

            %%% bgSubMethod test.
            if strcmp( bgSubMethod, 'alternative' )
                assert(Nbg == 2, 'number of background images must be two for ''%s'' bgSubMethod options',bgSubMethod);
            end


            %%% FOV def
            samFilesUsed = samFilesLocal( shotNumber );
            thisFile = samFilesUsed(1);
            %             loadImg = imread( fullfile(thisFile.folder, thisFile.name) );
            
            darkImg = 0.*single(imread(darkDir));
            warning('not using dark fields');
            loadImg = imreadAndRotation( fullfile(thisFile.folder, thisFile.name) ,rotAngle, darkImg);
            if isempty(obj.FOV)
                [yy,xx] = size(loadImg);

                obj.FOV = struct;
                obj.FOV.xlim = [1,xx];
                obj.FOV.ylim = [1,yy];
                obj.FOV.bgAddXFOV = floor( xx * (0.05 / 1.1) );
                obj.FOV.size = [yy,xx];
                obj.FOV.originalSize = [yy,xx];
            else
                assert(all(obj.FOV.originalSize == size(loadImg)),'setFOV is wrong')
            end


            %%% cropping range def.
            if all( obj.positionOffs(:) == 0 )
                minCuts = [0,0];
                maxCuts = [0,0];
            else
                minCuts = floor( min( obj.positionOffs ) ); % < 0
                maxCuts = ceil ( max( obj.positionOffs ) ); % > 0

                minCuts( minCuts > 0 ) = 0;
                maxCuts( maxCuts < 0 ) = 0;

                %%% update FOV if required
                if obj.FOV.ylim(1) + minCuts(1) < 1
                    obj.FOV.ylim(1) = 1 - minCuts(1); % make obj.FOV.ylim(1) + minCuts(1) = 1;
                end
                if obj.FOV.xlim(1) + minCuts(2) < 1
                    obj.FOV.xlim(1) = 1 - minCuts(2);
                end
                if obj.FOV.ylim(2)+maxCuts(1) > obj.FOV.originalSize(1)
                    obj.FOV.ylim(2) =  obj.FOV.originalSize(1) - maxCuts(1); % make obj.FOV.ylim(1) + minCuts(1) =  obj.FOV.originalSize(1);
                end
                if obj.FOV.xlim(2)+maxCuts(2) > obj.FOV.originalSize(2)
                    obj.FOV.xlim(2) = obj.FOV.originalSize(2) - maxCuts(2);
                end

                % update FOV
                obj.FOV.size = [obj.FOV.ylim(2)-obj.FOV.ylim(1), obj.FOV.xlim(2)-obj.FOV.xlim(1)] + 1;

                % make obj.FOV.size(2) odd
                if mod( obj.FOV.size(2) ,2) == 0
                    obj.FOV.size(2) = obj.FOV.size(2) - 1;
                    obj.FOV.xlim(1) = obj.FOV.xlim(1) + 1;
                end

            end

            ylimVec = obj.FOV.ylim(1)+minCuts(1) : obj.FOV.ylim(2)+maxCuts(1);
            xlimVec = obj.FOV.xlim(1)+minCuts(2) : obj.FOV.xlim(2)+maxCuts(2);

            %%% memory allocation
            sam = zeros(length(ylimVec), length(xlimVec), Nshot,'like',loadImg);
            bg  = zeros(length(ylimVec), length(xlimVec), Nbg, 'single');

            %%% pp = 1
            sam(:,:,1) = loadImg( ylimVec, xlimVec );

            %%% bg load
            for pp = 1:Nbg
                thisFile = bgFilesLocal(pp);
                %                 loadImg = imread( fullfile(thisFile.folder, thisFile.name) );
                loadImg = imreadAndRotation( fullfile(thisFile.folder, thisFile.name), rotAngle, darkImg);
                bg(:,:,pp) = single( loadImg( ylimVec ,xlimVec ) );
            end
            %    figure, imagesc(bg(:,:,1)),axis image

            %%% sam load
            fprintf('Image loading ... '); tic;
            parfor pp = 1:Nshot
                thisFile = samFilesUsed(pp);
                %                 loadImg = imread( fullfile(thisFile.folder, thisFile.name) );
                loadImg = imreadAndRotation( fullfile(thisFile.folder, thisFile.name) ,rotAngle, darkImg);
                sam(:,:,pp) = loadImg( ylimVec, xlimVec );
                fprintf('%03d/%03d\n',pp,Nshot);
            end
            sam = single(sam);

            %%% bg sub. mode
            switch bgSubMethod
                case 'average'
                    %sam = sam./ mean(bg,3);
                    %bg_field=mean(bg,3);
                    mean_BG=mean(bg,3);
                    mean_BG_shift=mean_BG;
                    limit=5;
                    sub_pixel_bool=true;
                    pearson_bool=true;
                    tt=1;
                    coo=make_coo(size(mean_BG(1:end/3,:)));
                    coo_full=make_coo(size(mean_BG(:,:)));
                    mean_BG_FT=s_fft2(mean_BG(1:end/3,:));
                    coo1=gpuArray(single(coo{1}));
                    coo2=gpuArray(single(coo{2}));
                    sam_out=sam./mean_BG;
                    
                    subdivs=4;
                    
                    for tt_s=1:subdivs 
                        tt=round(size(sam,3)/subdivs*(tt_s-1)+1):round(size(sam,3)/subdivs*tt_s);
                        tt(1);
                        
                        sam_FT=s_fft2(sam(1:end/3,:,tt));
                        xi_min=0;
                        yi_min=0;
                        cost_min=-inf;
                        sam_FT=gpuArray(single(sam_FT));
                        mean_BG_FT=gpuArray(single(mean_BG_FT));

                        sam_FT=sam_FT./abs(sam_FT);
                        mean_BG_FT=mean_BG_FT./abs(mean_BG_FT);
                        filter_frq=sqrt(coo1.^2+coo2.^2)<0.1;%0.05
                        base_coo=filter_frq.*mean_BG_FT.*conj(sam_FT);
                        cost=@(shift) real(sum(exp((1i*2*pi)*(real(shift(1,:,:)).*coo1+real(shift(2,:,:)).*coo2)).*base_coo,'all'));
                        
                        displacement_0=gpuArray(single(zeros(2,1,length(tt(:)))));
                        
                        t_np=1;
                        x_n=displacement_0;
                        itt_num=10;
                        if ~AUTO_offset
                            itt_num=0;
                        end
                        val_history=zeros(itt_num,1);
                        if figShow
                            figure;
                        end
                        for itt=1:itt_num
                            [val, grad]=adiff(cost, displacement_0);
                            val_history(itt)=val;
                            [displacement_0,x_n,t_np] = FISTA(+grad*1e-3,displacement_0,x_n,t_np);
                            if figShow
                                subplot(1,2,1);plot(val_history(1:itt));
                                subplot(1,2,2);plot(displacement_0(1,:)); hold on;plot(displacement_0(2,:),'r');hold off;
                                drawnow;
                            end
                        end
                        %{
                        for xi=-4:0.3333:4
                            for yi=-4:0.3333:4
                                cost=sum(filter_frq.*mean_BG_FT.*conj(sam_FT).*exp(1i*2*pi*(coo1.*xi+coo2.*yi)),'all');

                                if abs(cost)>cost_min
                                    %real(cost)
                                    cost_min=real(cost);
                                    xi_min=xi;
                                    yi_min=yi;
                                end
                            end
                        end
                        %}
                        temp=mean_BG;
                        temp=real(s_ifft2(s_fft2(temp).*exp(1i*2*pi*(coo_full{1}.*displacement_0(1,:,:)+coo_full{2}.*displacement_0(2,:,:)))));
                        sam(:,:,tt)=sam(:,:,tt)./temp;

                        %xi_min
                        %yi_min

                    end
                        
                    bg_field=mean(bg,3);
                case 'alternative'
                    for bb = 1:2
                        bvec = bb:2:Nshot;
                        sam(:,:,bvec) = sam(:,:,bvec)./ bg(:,:,bb);
                    end
            end


            %%% Position shifts
            if ~all( obj.positionOffs(:) == 0 )
                error('stop not done for BG')
                % zigzag boundary: to supress the boundary effect in sub-pixel shifting.
                sam = cat(1,sam,flipud(sam));
                sam = cat(2,sam,fliplr(sam));

                % fff2 and gridCell
                sam = fft2(sam); % to supress the boundary effect in sub-pixel shifting.
                gridCell  = ndgrid_matSizeIn( size(sam,[1,2]),'normalized', 'centerZero_ifftshift' );

                %
                for zi = 1:1:Nshot
                    sam(:,:,zi) =...
                        circshift_KR( sam(:,:,zi), -obj.positionOffs(zi,:),...
                        'fourierInputTF', true,'gridCell',gridCell);
                end

                sam = ifft2(sam,'symmetric');
                sam( sam<0 ) = 0; % make negative intensity zero

                %%% undo zigzag boundary
                sam = sam(1:size(sam,1)/2, 1:size(sam,2)/2, :);

                %%% cropOut
                ycropRange = (-minCuts(1)+1):(size(sam,1)-maxCuts(1));
                xcropRange = (-minCuts(2)+1):(size(sam,2)-maxCuts(2));

                assert( all( [length(ycropRange), length(xcropRange)] == obj.FOV.size), 'something is wrong here' )
                sam = sam(ycropRange,xcropRange,:);
            end

            fprintf('done .. %.1f sec\n',toc)
        end

        function b = sinoRegistration(obj, b ,varargin)
            p = inputParser();

            %%% parameters for fminpowellForRegistration
            addParameter(p, 'boundary', 10*[[-1;1],[-1;1]], @ismatrix);
            addParameter(p, 'options', [], @isstruct);

            %%% parameters for outer iteration
            addParameter(p, 'binN', [5,5], @(a) isvector(a) || iscell(a));
            addParameter(p, 'iterN',  1, @(a) isscalar(a) || iscell(a) );
            addParameter(p, 'showPlotTF', true, @islogical);
            addParameter(p, 'sparseAngleLib',[],  @isvector);

            %%% parameters for phantomInput
            addParameter(p, 'known_positionOffs', []);


            %%% parsing
            parse(p, varargin{:});
            boundary_in  = p.Results.boundary;
            options      = p.Results.options;  % use optimset to set. Valid properties: tolX, MaxIter, Display

            binN_in       = p.Results.binN;
            iterN_in      = p.Results.iterN;
            showPlotTF  = p.Results.showPlotTF;
            known_positionOffs   = p.Results.known_positionOffs;
            sparseAngleLib       = p.Results.sparseAngleLib;

            phantomInputTF       = any( ~isempty(known_positionOffs) );
            sparseAngleLibTF     = ~isempty(sparseAngleLib);

            %%% binN and iterN check
            if ~iscell(binN_in)
                binN_in  = {binN_in};
                iterN_in = {iterN_in};
                binIterN = 1;

            else % if binN is cell,
                binIterN = length(binN_in);
                assert( iscell(iterN_in), 'iterN_in must be a cell, if binN_in is a cell')
                assert( binIterN == length(iterN_in), 'the length of binN and iterN must be equal')
            end

            Nangle = size(b,3);
            totalIterN = sum( [iterN_in{:}] );

            if phantomInputTF
                errorRMS = zeros(1,totalIterN);
            end

            if showPlotTF
                dx_RMS   = zeros(1,totalIterN);
                plotx = 1:Nangle;
                figure(1551561),
            end

            %%% realTF: ifft2 symmetry option
            realTF   = isreal(b);

            %%% reloadTF: if registration is done for raw images
            reloadTF = ismember(obj.fieldReconMethod,{'proj','PC','KK'}) && realTF && ~phantomInputTF;
            if ~reloadTF
                b0 = b;
            end

            %%% binIter start
            totalIter = 0;
            fprintf('Image registration ... \n');
            for binIter = 1:binIterN
                binN  = binN_in { binIter };
                iterN = iterN_in{ binIter };
                boundary = boundary_in ./ binN.';

                %%% binning
                b =  matbinning(b, binN,'symmetric'); % binN must be odd here.
                % tomoHandle.show(angle(b));

                %%% gpu?
                % CPU is much faster for small images
                if obj.gpuInd && numel(b(:,:,1)) > 2^16
                    b = gpuArray(b);
                end
                Fb = fft2(b);

                %%% grid gen.
                sizeMat = cast(size(b,[1,2]),'like',real( b(1)) );
                gridCell  = ndgrid_matSizeIn(sizeMat,'normalized', 'centerZero_ifftshift');

                %%% iteration init
                measured_positionOffs = zeros(Nangle,2);
                %                 refMask  = true( size(b,[1,2] ),'like', b(1) == 0 );
                zeroMask  = all( abs(b)~=0 ,3);
%                 refBoundary = zeros(2); % [ y_low, y_up ; x_low, x_up ]

                %%% weight Matrix calc.
                if sparseAngleLibTF
                    [~,weight] = refimgFromSino(obj, b, 1, [], sparseAngleLib);
                else
                    [~,weight] = refimgFromSino(obj, b, 1, []);
                end

                %% iteration start
                refMask = zeroMask;
                for iter = 1:iterN
                    tic;
                    minus_dx  = zeros(Nangle,2);

                    %%% rotation axis centering using 0 <-> 180 degs
                    if ~strcmp(obj.fieldReconMethod,'KK')
                        testR = refimgFromSino(obj, b, 1, [], 180); % 180 deg

                        testF = Fb(:,:,1); % 0 deg
                        %                     testF = fft2( testF );
                        x1 = fminpowellForRegistration( boundary,options,testF,testR,refMask,'gridCell',gridCell);

                        axisOffsetLocal = [0, -x1(2)/2]; % get axisOffset
                        b   = circshift_KR( b, [ -axisOffsetLocal , 0 ],'gridCell',gridCell); % apply axisOffset
                        Fb  = fft2(b);

                        %%% global shift application
                        measured_positionOffs = measured_positionOffs + axisOffsetLocal.*binN;
                        %%% refMask update
                        refMask = refMaskCalc( zeroMask, measured_positionOffs, binN );
                    end

                    %%% registration frame-by-frame
                    parfor ii = 1:Nangle

                        %%% refcalc.
                        if sparseAngleLibTF
                            testR = refimgFromSino(obj, b, ii, weight, sparseAngleLib);
                        else
                            testR = refimgFromSino(obj, b, ii, weight);
                        end
                        %   subplot(121),imagesc(angle(testR)),axis image;colorbar; title(ii)
                        %   subplot(122),imagesc(angle(b(:,:,ii))),axis image;colorbar;

                        %%% registration
                        testF =  Fb(:,:,ii);
                        x1 = fminpowellForRegistration( boundary,options,testF,testR,refMask,'gridCell',gridCell );
                        minus_dx(ii,:)  = x1.';

                    end % angle update ends

                    %%% sino update
                    for ii = 1:1:Nangle
                        testF =  Fb(:,:,ii);
                        testF =  circshift_KR(testF, minus_dx(ii,:) , 'fourierInputTF', true );
                        Fb(:,:,ii) = testF;    % update
                    end

                    if realTF
                        b = ifft2(Fb,'symmetric');
                    else
                        b = ifft2(Fb);
                    end

                    %%% measured_positionOffs update
                    scaled_minus_dx  = minus_dx .* binN ;
                    measured_positionOffs = measured_positionOffs - scaled_minus_dx;

                    %%% refMask update
                    refMask = refMaskCalc( zeroMask, measured_positionOffs, binN );
                    totalIter = totalIter + 1;

                    %%% show plot
                    if showPlotTF
                        %%% RMS calc
                        scaled_minus_dx  = scaled_minus_dx - mean(scaled_minus_dx,1);
                        dx_RMS(totalIter)     = rms(scaled_minus_dx(:));
                        accumOffs = obj.positionOffs + measured_positionOffs;

                        %%% show
                        if phantomInputTF
                            subplot(521), plot(plotx,known_positionOffs(:,2), plotx,accumOffs(:,2)); axis tight;title(['x-axis offset, iter=',num2str(iter)]);
                            subplot(522), plot(plotx,accumOffs(:,2)-known_positionOffs(:,2)); axis tight;title('x-axis offset diff');
                        else
                            subplot(5,2,[1,2]), plot(plotx,accumOffs(:,2)); axis tight; title(['x-axis offset, iter=',num2str(iter)]);
                        end

                        if realTF
                            subplot(5,2,[3,4]), imagesc(squeeze( b(mcoor(size(b,1)),:,:) ));xlabel('rotAngle'); ylabel('x-axis');
                        else
                            subplot(5,2,[3,4]), imagesc(squeeze( angle(b(mcoor(size(b,1)),:,:)) ));xlabel('rotAngle'); ylabel('x-axis');
                        end

                        if phantomInputTF
                            subplot(525), plot(plotx,known_positionOffs(:,1), plotx,accumOffs(:,1));axis tight;   title(['y-axis offset, iter=',num2str(iter)]);
                            subplot(526), plot(plotx,accumOffs(:,1)-known_positionOffs(:,1)); axis tight; title('y-axis offset diff');
                        else
                            subplot(5,2,[5,6]), plot(plotx,accumOffs(:,1));axis tight;   title(['y-axis offset, iter=',num2str(iter)]);
                        end

                        if realTF
                            subplot(5,2,[7,8]), imagesc(squeeze( b(:,mcoor(size(b,2)),:) )); xlabel('rotAngle'); ylabel('y-axis')
                        else
                            subplot(5,2,[7,8]), imagesc(squeeze( angle(b(:,mcoor(size(b,2)),:)) )); xlabel('rotAngle'); ylabel('y-axis')
                        end

                        if phantomInputTF
                            subplot(5,2,9),semilogy(1:totalIter,dx_RMS(1:totalIter));  axis tight; xlabel('iteration #'); title(['difference from the previous iteration RMS (pixels): ', num2str(dx_RMS(totalIter),'%0.4f')])

                            errorRMS_temp = accumOffs-known_positionOffs;
                            errorRMS_temp = errorRMS_temp-mean(errorRMS_temp,1);
                            errorRMS(totalIter) = rms(errorRMS_temp(:));

                            subplot(5,2,10), plot(1:totalIter,errorRMS(1:totalIter)); axis tight; xlabel('iteration #'); title( ['position error RMS (pixels): ', num2str(errorRMS(totalIter),'%0.4f')] )
                        else
                            subplot(5,2,[9,10]), semilogy(1:totalIter,dx_RMS(1:totalIter));  axis tight; xlabel('iteration #'); title(['difference from the previous iteration RMS (pixels): ', num2str(dx_RMS(totalIter),'%0.4f')])
                        end

                        drawnow;
                        set(gcf,'Color','w');
                        %     pause;
                    end
                    fprintf('iter: %02d / %02d .. %.1f sec\n',iter,iterN,toc)
                end % iter ends

                %% apply the results
                obj.positionOffs = obj.positionOffs + measured_positionOffs;

                if reloadTF
                    b = obj.getSinogram(); % reload the sinogram
                else
                    b = b0;
                    zeroMask  = all( abs(b)~=0 ,3);
                    % figure,imagesc(zeroMask),axis image
                    b = fft2(b);
                    for zi = 1:1:Nangle
                        b(:,:,zi)  = circshift_KR( b(:,:,zi), -obj.positionOffs(zi,:) , 'fourierInputTF', true );
                    end
                    
                    if realTF
                        b = ifft2(b,'symmetric');
                    else
                        b = ifft2(b);
                    end

                    refMask = refMaskCalc( zeroMask, obj.positionOffs, [1,1] );
                    b = b.*refMask;
                    % figure,imagesc(refMask ),axis image
                    %   tomoHandle.show(angle(b.*refMask)); 
                    if phantomInputTF
                        b(b<0) = 0;
                    end
                end

            end % binIter ends

            %% nested funcions
            function refMask = refMaskCalc(zeroMask, measured_positionOffs, binN)                
                min_positionOffs = min(measured_positionOffs,[],1); 
                max_positionOffs = max(measured_positionOffs,[],1);    

                ylCut = floor( min_positionOffs(1) / binN(1) ); 
                xlCut = floor( min_positionOffs(2) / binN(2) );
                yuCut = ceil ( max_positionOffs(1) / binN(1) );
                xuCut = ceil ( max_positionOffs(2) / binN(2) );

                zeroMask = ...
                    circshift(zeroMask,[-ylCut,0]) &...
                    circshift(zeroMask,[-yuCut,0]) &...
                    circshift(zeroMask,[0,-xlCut]) &...
                    circshift(zeroMask,[0,-xuCut]);

                rectMask = false( size(zeroMask) );

                xlCut = -min(xlCut,0); % >= 0
                ylCut = -min(ylCut,0); % >= 0
                xuCut = max(xuCut,0);  % >= 0
                yuCut = max(yuCut,0);  % >= 0

                rectMask( (ylCut+1):(end-yuCut) , (xlCut+1):(end-xuCut) ) = true;
                refMask = zeroMask & rectMask;
            end
        end

        function [imgOut,weight] = refimgFromSino(obj,bField,refAngleInd,weight,sparseAngleLib)
            %% parsing
            if nargin < 4
                weight = [];
            end

            if nargin < 5
                sparseAngleLib = [];
            end

            %%
            [yy,xx,zz] = size(bField);

            rotAngleMax_local = obj.rotAngleMax;
            rotAnglelib_local = obj.rotAnglelib;
            drotAngle = rotAnglelib_local(2) - rotAnglelib_local(1);

            rotAnglelib_local = round(rotAnglelib_local / drotAngle); % make this integer

            %  make the frame of interest zero
            %             test = bField(:,:,refAngleInd);
            bField(:,:,refAngleInd) = 0;


            % remove 360 deg if it is duplicated.
            if rotAngleMax_local == 360
                rotAnglelib_local(end) = []; % erase the 360 deg.

                zz = zz-1;
                rotWeight    = ones(1,zz);

                if refAngleInd == 1
                    bField(:,:,1) = bField(:,:,end);
                elseif refAngleInd == zz

                else
                    bField(:,:,1)  = bField(:,:,1) + bField(:,:,end); % merge 0 and 360 degs.
                    rotWeight(1) = 2;
                end
                bField(:,:,end)  = []; % erase the frame of 360 deg.

            else % rotAngleMax_local == 180, no duplicated angles
                rotWeight    = ones(1,zz);
            end

            %% make the refAngle 0 deg
            bField    = circshift( bField,    [0,0,-refAngleInd+1] ); % make refAngleInd the first
            rotWeight = circshift( rotWeight, -refAngleInd+1 );

            if strcmpi(obj.fieldReconMethod,'KK')
                assert(rotAngleMax_local == 360, 'rotAngleMax must be 360 for KK')

                Nind90  = floor(90/drotAngle) + 1;
                part1Ind = 1:1:Nind90;
                part2Ind = [1, zz:-1:(zz-Nind90+2)]; % duplicate the first one

                bField    = bField(:,:,part1Ind) + bField(:,:,part2Ind);
                rotWeight = rotWeight(part1Ind) + rotWeight(part2Ind);

                bField = bField ./ reshape(rotWeight,1,1,[]);
                rotAnglelib_local = rotAnglelib_local(part1Ind); % [0,90)

                if (refAngleInd ~= 1) && (refAngleInd ~= zz)
                    rotWeight(1) = 0;   % first frame is empty unless refAngle == 0 or 360
                end
            else
                switch rotAngleMax_local
                    case 180

                        % start from 2 because there is no duplicated angle in 180 deg case
                        part1Ind = 2:zz-refAngleInd+1;
                        part2Ind = zz:-1:(zz-refAngleInd+2);
                        rotAnglelib_local = rotAnglelib_local(2:end); % (0,180]

                        Npart1 = length(part1Ind);
                        Npart2 = length(part2Ind);


                        if Npart1 > Npart2
                            bField_temp = bField(:,:,part1Ind);
                            bField_temp(:,:,1:Npart2) = bField_temp(:,:,1:Npart2) + bField(:,:,part2Ind);

                            rotWeight_temp = rotWeight(part1Ind);
                            rotWeight_temp(1:Npart2) = rotWeight_temp(1:Npart2) + rotWeight(part2Ind);

                            bField = bField_temp ./ reshape(rotWeight_temp,1,1,[]); % size(bField,3) = Npart1
                            bField    = cat(3, bField, zeros(yy,xx,zz-1-Npart1,'like', bField(1) ) ); % size(bField,3) = zz-1
                            rotWeight = [rotWeight_temp, zeros(1,zz-1-Npart1)];

                        else
                            bField_temp = bField(:,:,part2Ind);
                            bField_temp(:,:,1:Npart1) = bField_temp(:,:,1:Npart1) + bField(:,:,part1Ind);

                            rotWeight_temp = rotWeight(part2Ind);
                            rotWeight_temp(1:Npart1) = rotWeight_temp(1:Npart1) + rotWeight(part1Ind);

                            bField = bField_temp ./ reshape(rotWeight_temp,1,1,[]);% size(bField,3) = Npart2
                            bField = cat(3, bField, zeros(yy,xx,zz-1-Npart2,'like', bField(1) ) ); % size(bField,3) = zz-1
                            rotWeight = [rotWeight_temp, zeros(1,zz-1-Npart2)];
                        end


                    case 360
                        Nind180  = floor(180/drotAngle) + 1;
                        part1Ind = 1:1:Nind180;
                        part2Ind = [1,zz:-1:Nind180];  % share the 180 deg

                        bField    = bField(:,:,part1Ind) + bField(:,:,part2Ind);
                        rotWeight = rotWeight(part1Ind) + rotWeight(part2Ind);

                        bField = bField ./ reshape(rotWeight,1,1,[]);
                        rotAnglelib_local = rotAnglelib_local(part1Ind);
                        if (refAngleInd ~= 1) && (refAngleInd ~= zz)
                            rotWeight(1) = 0;   % first frame is empty unless refAngle == 0 or 360
                        end
                end
            end

            %% angle selection for ref. generation
            if ~isempty(sparseAngleLib)
                [rotAnglelib_local, ia] = intersect(rotAnglelib_local, round(sparseAngleLib/drotAngle));
                bField    = bField(:,:,ia);
                rotWeight = rotWeight(ia);
            end
            rotAnglelib_local = rotAnglelib_local.*drotAngle;
            Nxz       = length(rotAnglelib_local);

            %% jinc convolution matrix calculation
            if isempty(weight)
                %                 assert(mod(xx,2) == 1, 'horizontal size must be odd')

                kxx       = (1:xx)' - mcoor(xx);
                kxcoor    = kxx * cosd( rotAnglelib_local );
                kzcoor_sq = ( kxx * sind( rotAnglelib_local ) ).^2;

                weight = zeros( xx, xx*Nxz,'like', real(bField(1)));

                for ii = 1:xx
                    krho_sq = (kxcoor - kxx(ii)).^2 + kzcoor_sq;
                    closeInd = krho_sq < 0.9;
                    %     figure,imagesc(closeInd),axis image

                    if ~isempty(closeInd)
                        % jinc calc.
                        x = pi * sqrt( krho_sq (closeInd) ); % unitless
                        zeroportion = (abs(x)==0);
                        jincVal =  2 * besselj(1,x) ./ x;
                        jincVal(zeroportion)  = 1;
                        % plot(x,jincVal)

                        % weight alloc.
                        weight(ii, closeInd) = jincVal; % upper half
                        %                         weight(xx-ii+1, flipud(closeInd) ) = jincVal; % lower half
                    end
                end
                weight = reshape(weight,[xx,xx,Nxz]); % [xx,xx,Nxz];
            end
            % imagesc(weight)

            %             imagesc(weight3d(:,:,199));

            %% convolution operation
            realTF = isreal(bField);
            bField = permute(bField, [2,3,1]); % [yy,xx,zz] -> [xx,zz,yy]
            FbField = fftshift(fft(ifftshift( bField , 1), [], 1), 1); % x-direction fft

            Ftest = reshape(weight, xx, []) *  reshape(FbField, [], yy);  % kxx, yy

            %%% accumulated weight
            nonEmptyAngles = rotWeight > 0.5; % erase the weight from the empty frames.
            cumWeight = sum( abs( weight(:,:, nonEmptyAngles ) ) ,[2,3]);

            validkXind = cumWeight > 0;
            % plot(cumWeight)

            Ftest(validkXind,:) = Ftest(validkXind,:)./ cumWeight(validkXind);
            %             Ftest(mcoor(xx)+1:end,:) = conj( Ftest(mcoor(xx)-1:-1:1,:));
            % subplot(121),imagesc(abs(Ftest )),axis image

            if realTF
                imgOut = fftshift(ifft(ifftshift( Ftest , 1), [], 1,'symmetric'), 1); % x-direction ifft
            else
                imgOut = fftshift(ifft(ifftshift( Ftest , 1), [], 1), 1); % x-direction ifft
            end
            imgOut = permute(imgOut,[2,1]); % [xx,yy] -> [yy,xx]

            % subplot(121),imagesc(real(imgOut)),axis image;colorbar; title(refAngleInd)
            % subplot(121),imagesc(angle(imgOut)),axis image;colorbar; title(refAngleInd)
            % subplot(122),imagesc(real(imgOut/norm(imgOut(:))-test/norm(test(:)))),axis image;colorbar; title(refAngleInd)


        end
        function Sinogram = convSinogram(obj,Sinogram, zigzagTF)
            if nargin < 3
                zigzagTF = false;
            end
            objNA  = obj.setup.detector.objectiveLens.NA;
            objMag = obj.setup.detector.objectiveLens.magnification;
            visWl  = obj.setup.detector.scintillator.wavelength;
            camBin = obj.setup.detector.camera.binning;
            camPix = obj.setup.detector.camera.pixelSize;
            try
                PSF    = obj.setup.detector.PSF;
                PSF = fullfile( fileparts(obj.setup.jsonPath), PSF);
            catch
                PSF=[]
                warning('no scyntilator psf !!!')
            end
            %             PSFpix = obj.setup.detector.PSF.pixelSize;
            gpuIndlocal   = obj.gpuInd;

            %%% params def
            scintPix = camPix*camBin/objMag;

            if zigzagTF
                % zigzag boundary: to supress the boundary effect in sub-pixel shifting.
                Sinogram = cat(1,Sinogram,flipud(Sinogram));
                Sinogram = cat(2,Sinogram,fliplr(Sinogram));
            end

            [sy,sx,sz] = size(Sinogram);

            
            cropSize1=[sy,sx];
            scintOTF=make_bright_field_OTFv2(cropSize1,scintPix.*[1 1],obj.setup.detector.objectiveLens.NA,obj.setup.detector.scintillator.wavelength,obj.setup.detector.scintillator.thickness,obj.setup.detector.scintillator.RI,0,1.5,1);%0 um aberation
            scintPix2=scintPix;
           
            fprintf('FourierCrop: [%d, %d] -> [%d, %d], effBin: %.2f\n', sy, sx, cropSize1, scintPix2/scintPix )
            %%% lowpass filter
            Sinogram = mcrop(fftshift(fft2(ifftshift(Sinogram))), cropSize1);
            Sinogram = fftshift(ifft2(ifftshift(Sinogram),'symmetric')) * prod(cropSize1)/sy/sx ;            
            % figure, contrastVis(Sinogram);

            %%% size update
            signalkR = 2 * objNA /visWl * scintPix2;   % pixel^-1, < 1
            signalkR = min( signalkR, 1/2 );
            
            %%% deconvolution
            if gpuIndlocal
                if isgpuarray(Sinogram)
                    gpu=gpuDevice;
                else
                    gpu=gpuDevice(gpuIndlocal);
                    size(Sinogram)
                    Sinogram = gpuArray(Sinogram);
                end

                bytesPerAngle = (sy*sx) * 8; % 8 bytes for complex single
                maxChunkSize = floor( gpu.AvailableMemory / bytesPerAngle / 10 ); % last factor is a heuristic safety value.
                iterN = ceil( sz / maxChunkSize);
                chunks = round( linspace(1, sz+1, iterN+1) );

                meanSNRMap = zeros(size(Sinogram,[1,2]),'like',real(Sinogram(1)));
                for ii = 1:iterN
                    thisChunk = chunks(ii) : chunks(ii+1)-1;

                    Sinogram(:,:,thisChunk) = mfft2(gpuArray(Sinogram(:,:,thisChunk)));
                    [Sinogram(:,:,thisChunk),meanSNRMap_temp] =...
                        convOTF( Sinogram(:,:,thisChunk), scintOTF, signalkR );

                    meanSNRMap = meanSNRMap + meanSNRMap_temp;
                    fprintf('Deconv recon. finished: %d / %d \n', ii,iterN);
                end

            else % no gpu
                [Sinogram,meanSNRMap] = convOTF( mfft2( Sinogram ), scintOTF, signalkR );
            end

            Sinogram = abs(Sinogram);
            if zigzagTF
                %%% undo zigzag boundary
                Sinogram = Sinogram(1:size(Sinogram,1)/2, 1:size(Sinogram,2)/2, :);
            end

            %%% meanSNRMap def.
            meanSNRMap = meanSNRMap / iterN;
            obj.SNRMap = gather(meanSNRMap);
            %             figure, imagesc(log10(abs(meanSNRMap))),axis image

        end
        %%% GPU is used if posssible
        function Sinogram = deconvSinogram(obj,Sinogram, zigzagTF)
            if nargin < 3
                zigzagTF = false;
            end
            objNA  = obj.setup.detector.objectiveLens.NA;
            objMag = obj.setup.detector.objectiveLens.magnification;
            visWl  = obj.setup.detector.scintillator.wavelength;
            camBin = obj.setup.detector.camera.binning;
            camPix = obj.setup.detector.camera.pixelSize;
            try
                PSF    = obj.setup.detector.PSF;
                PSF = fullfile( fileparts(obj.setup.jsonPath), PSF);
            catch
                PSF=[]
                warning('no scyntilator psf !!!')
            end
            %             PSFpix = obj.setup.detector.PSF.pixelSize;
            gpuIndlocal   = obj.gpuInd;

            %%% params def
            scintPix = camPix*camBin/objMag;

            if zigzagTF
                % zigzag boundary: to supress the boundary effect in sub-pixel shifting.
                Sinogram = cat(1,Sinogram,flipud(Sinogram));
                Sinogram = cat(2,Sinogram,fliplr(Sinogram));
            end

            [sy,sx,sz] = size(Sinogram);

            %NAr = objNA / visWl * [sy,sx] * scintPix; % pixel
            %scintOTF = scintOTFgen(NAr, 4, [sy,sx] ); % factor 4 is a heuristic value
            %p=struct;
            %p.calOTF=ones(2000,2000,'single');
            %p.imgPix=scintPix;

            % imagesc(abs(loadPSF.calOTF))

            %scintPix_target = max(scintPix, p.imgPix);
            %cropSize0 = [sy,sx] * (scintPix/scintPix_target);
            %cropSize1 = floor(cropSize0);
            %if mod(cropSize1,2) == 0 % make it odd
            %    cropSize1 = cropSize1 - 1;
            %end
            %addMag = cropSize0./cropSize1; % very close to 1 but >1
            %scintPix2 = mean([sy,sx]./cropSize1 * scintPix );
            % should be similar to scintPix_target, but can be a bit smaller due to floor()

            %scintOTF = imMagRot(p.calOTF, addMag, 0, size(p.calOTF));
            %scintPSF = msize( fftshift(ifft2(ifftshift(scintOTF))), cropSize1 );
            %scintOTF = fftshift(fft2(ifftshift(scintPSF)));
            
            %scintOTF=make_bright_field_OTF(cropSize1,2*4.5/20.*[1 1],0.4,0.532);
            %scintOTF=make_bright_field_OTFv2(cropSize1,0.470.*[1 1],0.4,0.532,20,1.9,500,1.5,1);%500um aberation
            %scintOTF=make_bright_field_OTFv2(cropSize1,0.470.*[1 1],0.4,0.532,20,1.9,0,1.5,1);%0 um aberation
            cropSize1=[sy,sx];
            scintOTF=make_bright_field_OTFv2(cropSize1,scintPix.*[1 1],obj.setup.detector.objectiveLens.NA,obj.setup.detector.scintillator.wavelength,obj.setup.detector.scintillator.thickness,obj.setup.detector.scintillator.RI,0,1.5,1);%0 um aberation
            scintPix2=scintPix;
            %scintPix
            %error('stop')
            
            %error('use the setup info to create')
%            figure; imagesc(abs(scintOTF)); axis image;
            
            fprintf('FourierCrop: [%d, %d] -> [%d, %d], effBin: %.2f\n', sy, sx, cropSize1, scintPix2/scintPix )
            %%% lowpass filter
            Sinogram = mcrop(fftshift(fft2(ifftshift(Sinogram))), cropSize1);
            Sinogram = fftshift(ifft2(ifftshift(Sinogram),'symmetric')) * prod(cropSize1)/sy/sx ;            
            % figure, contrastVis(Sinogram);

            %%% size update
            signalkR = 2 * objNA /visWl * scintPix2;   % pixel^-1, < 1
            signalkR = min( signalkR, 1/2 );
            
            %%% deconvolution
            if gpuIndlocal
                if isgpuarray(Sinogram)
                    gpu=gpuDevice;
                else
                    gpu=gpuDevice(gpuIndlocal);
                    size(Sinogram)
                    Sinogram = gpuArray(Sinogram);
                end

                bytesPerAngle = (sy*sx) * 8; % 8 bytes for complex single
                maxChunkSize = floor( gpu.AvailableMemory / bytesPerAngle / 10 ); % last factor is a heuristic safety value.
                iterN = ceil( sz / maxChunkSize);
                chunks = round( linspace(1, sz+1, iterN+1) );

                meanSNRMap = zeros(size(Sinogram,[1,2]),'like',real(Sinogram(1)));
                for ii = 1:iterN
                    thisChunk = chunks(ii) : chunks(ii+1)-1;

                    Sinogram(:,:,thisChunk) = mfft2(gpuArray(Sinogram(:,:,thisChunk)));
                    [Sinogram(:,:,thisChunk),meanSNRMap_temp] =...
                        deconvOTF( Sinogram(:,:,thisChunk), scintOTF, signalkR );

                    meanSNRMap = meanSNRMap + meanSNRMap_temp;
                    fprintf('Deconv recon. finished: %d / %d \n', ii,iterN);
                end

            else % no gpu
                [Sinogram,meanSNRMap] = deconvOTF( mfft2( Sinogram ), scintOTF, signalkR );
            end

            Sinogram = abs(Sinogram);
            if zigzagTF
                %%% undo zigzag boundary
                Sinogram = Sinogram(1:size(Sinogram,1)/2, 1:size(Sinogram,2)/2, :);
            end

            %%% meanSNRMap def.
            meanSNRMap = meanSNRMap / iterN;
            obj.SNRMap = gather(meanSNRMap);
            %             figure, imagesc(log10(abs(meanSNRMap))),axis image

        end

        function Sinogram = sinoNorm(obj,Sinogram, varargin)
            p = inputParser();
            addParameter(p, 'logInput', true);
            addParameter(p, 'phaseRampCorrection', false);
            addParameter(p, 'phaseCorrDims', 2);
            % phase correlation dimensions
            % put 2 (x-axis) for KK
            %    --> KK does not provide y-axis global phase
            % put [1,2] (xy-plane) for SSM

            parse(p, varargin{:});
            logInput             = p.Results.logInput;
            phaseRampCorrection  = p.Results.phaseRampCorrection;
            phaseCorrDims        = p.Results.phaseCorrDims;

            %%% bgMask init.
            [yy,xx,zz] = size(Sinogram);
            bgMask = false(1,xx);
            bgMask([1:obj.FOV.bgAddXFOV, end-obj.FOV.bgAddXFOV+1:end]) = true;
            % plot(bgMask)

            %%% zero input consideration
            if logInput
                zeroMask = (imag(Sinogram) == Inf);
            else
                zeroMask = abs(Sinogram) == 0;
            end
            
            sinogramHasZeros = any( zeroMask ,"all" );

            if sinogramHasZeros                
                nonZeroMask = all(~zeroMask,3);                
                bgMask = nonZeroMask & bgMask;
                %                 clear zeroMask;
                %   tomoHandle.show(zeroMask);
            end
            clear zeroMask;

            %%% phaseRampCorrection
            if phaseRampCorrection
                outGridCell = ndgrid_matSizeIn(size(Sinogram,[1,2]),1,'centerZero');
                YY = outGridCell{1};
                XX = outGridCell{2};

                %%% relative phase ramp correction: 
                for zi = 1:zz-1
                    testRamp = exp(1i*angle(  conj(Sinogram(:,:,zi)) .* Sinogram(:,:,zi+1) ));
                    [~,centerkvec]= PhiShiftKR( testRamp, 'bgMask', bgMask );

                    phaseRamp = exp(-1i*2*pi* ( YY*centerkvec(1) + XX*centerkvec(2) ) );
                    meanPhase = testRamp .* phaseRamp;
                    meanPhase = angle(mean(meanPhase(bgMask)));

                    Sinogram(:,:,zi+1) = Sinogram(:,:,zi+1) .* phaseRamp .* exp(-1i*meanPhase);
                    %                     imagesc(angle(Sinogram(:,:,zi+1))),axis image;title(zi);
                    %                     drawnow;
                end

                %%% mean phaseRamp calc. --> phaseRamp global apply
                meanSino = mean(exp(1i*angle(Sinogram)),3);
                [~,centerkvec]= PhiShiftKR( meanSino, 'bgMask', bgMask );
                phaseRamp = exp(-1i*2*pi* ( YY*centerkvec(1) + XX*centerkvec(2) ) );
                Sinogram = Sinogram.* phaseRamp;
                % imagesc(angle(meanSino)),axis image;
                % imagesc(angle(goodImg)),axis image;


                %%% relative phase ramp correction:
%                 for zi = 1:zz
%                    Sinogram(:,:,zi) = PhiShiftKR( Sinogram(:,:,zi), 'bgMask', bgMask );
%                 end
            end

            
            %%% bgE calc
            if sinogramHasZeros
                if all( phaseCorrDims == [1,2] )
                    Sinogram = reshape(Sinogram, [], zz);
                    bgE = zeros(1,1,zz,'like',Sinogram);
                    for zi = 1:zz
                        bgE(zi) = mean( Sinogram(bgMask,zi) );
                    end
                    Sinogram = reshape(Sinogram, yy,xx,zz);
                else
                    error('need to build')
                end

            else
                bgE = mean( Sinogram(:,bgMask,:), phaseCorrDims); % get BG mean along the phaseCorrDims.                %             plot(imag(bgE))
            end

            %%% bg normalizatoin
            if logInput
                Sinogram = Sinogram - bgE;
            else
                Sinogram = Sinogram ./ bgE;
            end
            %   tomoHandle.show(real(Sinogram));
            %   tomoHandle.show(abs(Sinogram+1));
        end

        function Sinogram = unwrap2_Lp_FieldIn(obj, Sinogram, p)
            if nargin < 3
                p = 0; % L0-norm default
            end
            
            imagPart = -log(abs(Sinogram));
            realPart = angle(Sinogram);
            zz = size(realPart,3);
            parfor zi = 1:zz
                realPart(:,:,zi) = unwrap2_Lp( realPart(:,:,zi), p ); %Lp-norm unwrap2
                fprintf('unwrapping... %05d / %05d\n',zi,zz)
            end
            % figure,tomoHandle.show(imagPart);
            Sinogram = realPart + 1i*imagPart;            

            %%% normalization
            Sinogram = obj.sinoNorm(Sinogram,'logInput',1, 'phaseRampCorrection',0 , 'phaseCorrDims', [1,2]);            
            Sinogram(imagPart == Inf) = 0;

        end

        function [Sinogram, optionalOut] = getOptField(obj, Sinogram, varargin)
            p = inputParser();
            addParameter(p, 'angleOutTF', true);
            addParameter(p, 'figShow', false);

            parse(p, varargin{:});
            angleOutTF         = p.Results.angleOutTF;
            figShow            = p.Results.figShow;

            %%% deconv.
            if obj.deconvTF
                Sinogram = obj.deconvSinogram(Sinogram);
                fprintf(' * Deconvolution done.\n');
            end
            % tomoHandle.show(Sinogram,95,'refFrame',1);  colormap(turbo)

            %%% field reconstruction
            switch obj.fieldReconMethod
                case {'proj', 'projection', 'PC'}
                    optionalOut = NaN;
                    if obj.gpuInd
                        if ~isgpuarray(Sinogram)
                            gpuDevice(obj.gpuInd);
                            Sinogram = gpuArray(Sinogram);
                        end
                    end

                case 'SSM'
                    %%% input
                    [~,~,ext] = fileparts(obj.samFiles(1).name);
                    inputExtIsMat = strcmp(ext,'.mat');
                    
                    %%% rotation axis angle correction
                    tiltAngle = -obj.setup.detector.camera.tiltAngle;
                    if isempty(tiltAngle)
                        tiltAngle = 0;
                    end

                    if ~inputExtIsMat
                        assert(ischar(Sinogram),'put paramsFile for sinogram for SSM');
                        d  = obj.SSMELibGenerator(Sinogram,figShow,ext,tiltAngle);

                    else % if inputExtIsMat
                        topDir = obj.samFiles(1).folder;
                        ELibGenForTomoFuncOutStr = 'ELibForTomo_*.mat';
                        ELibList = dir(fullfile(topDir,ELibGenForTomoFuncOutStr));
                        if isempty(ELibList) % only SSM done, ELib does not exist
                            d  = obj.SSMELibGenerator('',figShow,ext,tiltAngle);

                        elseif length(ELibList) == 1 %  SSM done --> ELib done
                            d = load(fullfile(topDir,ELibList(1).name));

                        else % multi ELibs exist
                            readFile = uigetfile(fullfile(topDir,'*.mat'),'Select ELib files');
                            d = load(fullfile(topDir,readFile));
                        end
                    end

                    %%% for SSM, return here
                    Sinogram = d.EsamLib;
                    d = rmfield(d,'EsamLib');
                    optionalOut = d;

                    %%% pixelsize def.
                    obj.setup.imagePixelSize = optionalOut.propaXPix;
                    

                    %%% gpu activation
                    if obj.gpuInd
                        Sinogram = gpuArray(Sinogram);
                    end
                    return

                case 'KK'
                    %%% default angleOutTF = false for KK
                    if ismember('angleOutTF',p.UsingDefaults)
                        angleOutTF = false;
                    end
                    ZPdia  = obj.setup.afterSample.zonePlate.diameter;
                    ZPdrn  = obj.setup.afterSample.zonePlate.outermostWidth;
                    imgPix = obj.setup.imagePixelSize; % nm
                    OCD    = obj.setup.afterSample.zonePlate.offAxis.distance;      % nm
                    offAngle = obj.setup.afterSample.zonePlate.offAxis.angleDeg;      % nm
                    bgAddXFOV =  obj.FOV.bgAddXFOV;

                    NAfrac  = 1/2* 1/ZPdrn * imgPix;       % unitless freq. (-1/2, 1/2]
                    OCDfrac = (OCD/ZPdia)*1/ZPdrn* imgPix; % unitless freq. (-1/2, 1/2]

                    [Sinogram,Fsupport] = obj.KKfieldRecon( ...
                        Sinogram,NAfrac,OCDfrac,offAngle,bgAddXFOV,obj.gpuInd,figShow,'ori');
                    optionalOut = Fsupport;
            end  % field reconstruction end

            %%% rough y-axis registration: disabled 22.03.23
            %%%  becuase of the introduction of
            %%%  sinoRegistration function before the field reconstruction.
            %             sucessFlag = 'R'; % R for retry
            %             while strcmpi(sucessFlag,'R')
            %                 [Sinogram,sucessFlag] = obj.yaxisRegistration(Sinogram);
            %             end

            %%% DC term normalization: beam intensity fluctuation + global phase matching
            DCval = mean(Sinogram, [1,2]);
            Sinogram = Sinogram ./ DCval;

            %%% angle (-1i*log) conversion of Fields
            if angleOutTF
                switch obj.fieldReconMethod
                    case {'proj', 'projection'}
                        % attnuation only
                        Sinogram = -real( log(abs(Sinogram)) ); % imag part only, positive ( beta > 0 )

                    case 'PC'
                        % for PC, this value becomes etimated phase
                        Sinogram = 1/2 * real( log(abs(Sinogram)) ); % real part only, negative (-delta < 0)

                    otherwise % KK
                        Sinogram = angle(Sinogram) - 1i*log(abs(Sinogram)); % -1i*RytovField (or compelx angle) -delta+1i*beta
                end

                %%% bg normalization of each measurements
                Sinogram = obj.sinoNorm(Sinogram); % make bg = 0; normalization
            end
            %             Fbfield = fftshift(fft(ifftshift(bField,2),[],2),2); % pix
        end

        function bField = KKAngleMerge(obj, bField)
            %%% warning('This function is only for ''KK'' fieldReconMethod')
            if ~strcmp(obj.fieldReconMethod,'KK')
                return;
            end

            [yy,xx,zz] = size(bField);
            xc = mcoor(xx);

            %%% assert
            offAngle = obj.setup.afterSample.zonePlate.offAxis.angleDeg;
            assert( offAngle == 0 || offAngle == 180, 'offAxis angle must be horizontal (0 or 180)');

            %%% merge
            bField = fftshift(fft(ifftshift( bField,2 ),[],2),2);
            tomoAngles =  obj.rotAnglelib;
            Fweight    =  zeros(1,xx,zz);
            measuredHalf = false(1,xx);
            if offAngle == 180
                measuredHalf(1:xc) = true;
            elseif offAngle == 0
                measuredHalf(xc:xx) = true;
            end
            Fweight(1,measuredHalf,:) = 1;
            %             tomoHandle.show(abs(bField));
            %             imagesc(squeeze(abs(Fweight))); axis image

            %%% add 360
            addInd = find( tomoAngles == 360 );
            if ~isempty(addInd)
                counterPartInd  = addInd - addInd(1) + 1;
                Fweight(1,measuredHalf,counterPartInd) = Fweight(1,measuredHalf,counterPartInd) + 1;
                bField(:,:,counterPartInd) = bField(:,:,counterPartInd) + bField(:,:,addInd);
            end

            %%% flipped add
            flipInd = find( tomoAngles>=180 & tomoAngles<360 );
            if ~isempty(flipInd)
                counterPartInd  = flipInd - flipInd(1) + 1;
                Fweight(1,~measuredHalf,counterPartInd) = Fweight(1,~measuredHalf,counterPartInd) + 1;
                Fweight(1,xc,counterPartInd) = Fweight(1,xc,counterPartInd) + 1;
                bField(:,:,counterPartInd) = bField(:,:,counterPartInd) + flip(bField(:,:,flipInd),2);
            end

            %%% average and truncation
            invalidAngle = tomoAngles >= 180;
            bField(:,:,invalidAngle) = [];
            Fweight(:,:,invalidAngle) = [];

            validInd = (Fweight>0);
            bField  = reshape(bField,yy,[]);
            Fweight = reshape(Fweight,1,[]);

            bField(:,validInd) = bField(:,validInd) ./ Fweight(validInd);
            bField = reshape(bField,yy,xx,[]);

            bField = fftshift(ifft(ifftshift( bField,2 ),[],2),2);
            %             tomoHandle.show(abs(bField))
            %             tomoHandle.show(angle(bField))

            %                         bField = obj.sinoNorm(bField) + 1; % v6.1: make the DC term 1
            %%% unwrap2 + log 
%             bField = obj.unwrap2_Lp_FieldIn(bField);

            %%% angle
            bField = angle(bField) - 1i*log(abs(bField)); % -1i*RytovField (or compelx angle) -delta+1i*beta
            %                 figure,  tomoHandle.show(imag(bField),99.99); colorbar
            %                 figure,  tomoHandle.show(real(bField),99.99); colorbar
            
            
            %%% bg nomarlization of each measurements
            bField = obj.sinoNorm(bField); % make bg = 0; normalization
        end

        function [d,samList,bgList]  = SSMELibGenerator(obj,paramsFile,showTF,ext,tiltAngle)
            %%% CZPNA calc.
            XrayWl  = Etowl( (obj.setup.source.energy_eV)*1e-3 );
            ZPdrn   = obj.setup.beforeSample.zonePlate.outermostWidth;
            OSADia     = obj.setup.beforeSample.pinhole.diameter;
            OSAdefocus = obj.setup.beforeSample.pinhole.defocus;

            CZPNA_original =  XrayWl/ZPdrn/2;
            CZPNA_OSAcut   = OSADia/OSAdefocus/2;
            CZPNA = min( CZPNA_original, CZPNA_OSAcut );
            
            %%% numerialBeamStop calc.
            numerialBeamStopRadius = 2.2331 / OSADia; % nm^-1, 2.2331 = 2nd zero of jinc

            %%% fileList struct
            fileList = struct('samList',obj.samFiles,'bgList',obj.bgFiles,'ext',ext);

            %%% run
            [d,samList,bgList] = ELibGenForTomo('paramsFile',paramsFile, ...
                'gpuInd',obj.gpuInd,...
                'showTF',showTF,'saveTF',1,'CZPNA',CZPNA,'fileList',fileList,'tiltAngle',tiltAngle,...
                'instabilityCompensationOn',0,'numerialBeamStopRadius',numerialBeamStopRadius);

            %%% revise the image files to generated .mat files
            obj.samFiles = samList;
            obj.bgFiles  = bgList;
        end

        function rotAxisTiltAngle = rotAxisTiltRegistration(obj,Sinogram,varargin)
            p = inputParser();

            %%% parameters for fminpowellForRegistration
            addParameter(p, 'boundary', []);
            addParameter(p, 'showTF', 0);

            %%% parsing
            parse(p, varargin{:});
            boundary  = p.Results.boundary;
            showTF    = p.Results.showTF;

            if isempty(boundary)
                yScanRange = [-10,10].'; % pixel
                xScanRange = [-10,10].'; % pixel
                AngleScanRange = [-10,10].'; % deg
                boundary = [yScanRange,xScanRange,AngleScanRange];
            end

            %%% load
            if strcmp(obj.fieldReconMethod,'KK')
                error('KK is not working here')
            end
            ind0   = (obj.rotAnglelib == 0);
            ind180 = (obj.rotAnglelib == 180);

            frame0   = Sinogram(:,:,ind0);
            frame180 = fliplr(Sinogram(:,:,ind180));
            % figure,imagesc(angle(shot0)),axis image
            % figure,imagesc(angle(shot180)),axis image

            %%% function def.
            refMask = (frame0 ~= 0) & (frame180 ~= 0);
            fun = @(vecIn) vecInCostOut(vecIn,frame180,frame0,refMask);            
            options = optimset('tolX',1e-2);

            %%% registration (shift & rotation)
            vIter = fminpowell( fun, boundary, options);
            % disp(vIter)

            if showTF
                [~,imgFinal] = vecInCostOut(vIter,frame180,frame0,refMask);
                figure(1141),imagesc(angle(frame0)),axis image
                figure(1142),imagesc(angle(imgFinal)),axis image
            end

            rotAxisTiltAngle = vIter(3)/2; %
            if isempty(obj.setup.detector.camera.tiltAngle)
                obj.setup.detector.camera.tiltAngle = rotAxisTiltAngle;
            end

            %%% nested function
            function [costOut,test] = vecInCostOut(vecIn, imgFloat, imgRef, refMask, varargin)
                if nargin < 4
                    refMask = true(size(imgRef));
                else
                    assert( all(size(refMask) == size(imgRef)), 'size of refMask must be equal to imgRef' )
                end

                %%% vecIn
                shiftVec = vecIn(1:2);
                rotAngle = vecIn(3);

                %%% rotation
                test =imMagRot(imgFloat,1,rotAngle);

                %%% shift
                test = circshift_KR(test, shiftVec, varargin{:} );
                %     tform = rigid2d(eye(2), [shiftVec(2), shiftVec(1)]);
                %     test = imwarp(imgFloat,tform,'OutputView',imref2d( size(imgFloat) ));
                %     figure, imagesc(angle(test)),axis image

                %%% ovelapMask calc.
                ceiledVec = ceil( abs(shiftVec) );
                yzeroArea = mod( sign(shiftVec(1)) * (1:1:ceiledVec(1)), size(imgFloat,1)+1 );
                xzeroArea = mod( sign(shiftVec(2)) * (1:1:ceiledVec(2)), size(imgFloat,2)+1 );

                %     overlapMask = true( size(fft_imgFloat),'like', refMask(1) );
                %     overlapMask(yzeroArea,:) = false;
                %     overlapMask(:,xzeroArea) = false;

                %%% costFunction calc: normalized correlation
                % refMask = refMask & circshift(refMask,sign(shiftVec).*ceiledVec);
                refMask(yzeroArea,:) = false;
                refMask(:,xzeroArea) = false;
                costOut = - abs(sum( conj(test(refMask)) .* imgRef(refMask) )) /norm(test(refMask)) /norm(imgRef(refMask));
            end

        end

        function bField = filterSinogram(obj,bField,outType,filterName)
            %%% default filter
            if nargin < 4
                filterName = 'wiener';
            end

            %%% real output only?
            if strcmp(outType,'real')
                bField = real(bField);
            elseif strcmp(outType,'imag')
                bField = -imag(bField);
            end
            %                 figure,  tomoHandle.show(real(test),99.99); colorbar

            %%% Sinogram filtering
            switch obj.tomoReconMethod
                case 'FBP'
                    imgPix      = obj.setup.imagePixelSize; % nm
                    xWindowSize = min( imgPix/obj.tomoRes(2) , 1 ); % NAcicle ratio

                    %%% SNRMap consideration
                    if isempty(obj.SNRMap)
                        SNR_kx = [];
                    else
                        % pre-filtering using SNR
                        SNR_kx = mean(obj.SNRMap,[1,3]);
                        % figure,plot(log10(abs(SNR_kx)))
                    end

                    %%% do filtering
                    %  figure,   tomoHandle.show(real(bField),99); axis image;
                    bField = iradon_filterPart(bField,filterName, xWindowSize,2,SNR_kx);
                    %                     bField = iradon_filterPart(bField,'wiener', xWindowSize,2);
                    %                     test = iradon_filterPart(bField,'ram-lak', xWindowSize,2);
                    %                     bField = iradon_filterPart(bField,'hamming', xWindowSize,2);

                    %                     test = iradon_filterPart(bField,'wiener', xWindowSize,2,SNR_kx);
                    %                     test = iradon_filterPart(bField,'wiener', xWindowSize,2);
                    %                     test = iradon_filterPart(bField,'ram-lak', xWindowSize,2);
                    %                     test = iradon_filterPart(bField,'hamming', xWindowSize,2);
                otherwise
                    error('NEED TO BUILD')
            end
            %           figure,  tomoHandle.show(real(test));
            %           figure,  imagesc(squeeze(test(1,:,:))),axis image
            %           figure,  imagesc(squeeze(bField(1,:,:))),axis image

            %%% bg nomarlization of each measurements
            %             bField = obj.sinoNorm(bField); % make bg = 0; normalization
        end

        function axisOffset = autoRotAxisScan(obj,bFieldSelectedY,varargin)
            [yy,xx,~] = size(bFieldSelectedY);
            assert(yy == 1,'autoRotAxisScan works for a single y-slice only')

            %% inputParsing
            p = inputParser();
            addParameter(p, 'filterName','wiener')
            addParameter(p, 'scanRange', [], @(v) isempty(v) || isvector(v))
            parse(p, varargin{:});
            scanRange          = p.Results.scanRange;     % pixel; empty input --> auto finder
            filterName         = p.Results.filterName;     
            %% init. Setting
            if isempty(scanRange)
                scanAmp    = min( round(xx/10) , 10 );
                scanRange  = [-scanAmp, scanAmp];
            else
                scanRange = sort(scanRange);
            end

            xR = (xx-1)/2;
            tomoMask = ~mk_ellipse(xR,xR,xx,xx); %% xx is odd

            %% scan using fminbnd ('golden section search, parabolic interpolation', or Brent's method)
            options = optimset('tolX',1e-2);
            axisOffset = fminbnd( @axisInEntropyOut ,scanRange(1), scanRange(2), options);
            axisOffset = gather(axisOffset);
            minImg = axisInTomoSliceOut( axisOffset );

            figure(15151),
            subplot(121),tomoHandle.show(minImg,99); colormap(turbo); axis image;
            title(['axisOffset: ',num2str(axisOffset)]);  colorbar;
            %             subplot(122),plot(scanPoints(readVec),imgEnt(readVec), scanPoints(minInd), imgEnt(minInd),'R*')
            drawnow;

            %% manual check
            QQ = 0;
            while true
                testImg = axisInTomoSliceOut(axisOffset);
                subplot(122),tomoHandle.show(testImg,99); axis image;
                title(['lastInput: ',num2str(QQ,'%+.2f'), ', axisOffset: ',num2str(axisOffset)]); colorbar;
                drawnow;
                QQbefore = QQ;
                QQ = input('rotation axis move digit (type 0 to stop): ');

                if isempty(QQ)
                    QQ = QQbefore;
                elseif QQ == 0
                    break;
                end

                axisOffset = axisOffset + QQ;
            end


            %% nested functions
            function [tomoSlice_out,histStep] = axisInTomoSliceOut(axisOffset_in)
                % rotation center
                test = obj.rotAxisShift(bFieldSelectedY, axisOffset_in, obj.FOV.bgAddXFOV);
                %     tomoHandle.show(real(test));  colormap(turbo)
                %     tomoHandle.show(angle(test));  colormap(turbo)

                % merge if KK
                test = obj.KKAngleMerge(test);
                
                ringRemovePrmLocal = obj.ringRemovePrm;
                if ~isempty(ringRemovePrmLocal)
                    test = obj.removeStripeOfSino(test, 'wavelet-fft', ...
                        ringRemovePrmLocal.waveletN, ringRemovePrmLocal.wname, ringRemovePrmLocal.sigma);
                end

                % filter
                test = obj.filterSinogram(test,'real',filterName);

                % outputs
                histStep = std(abs(test),0,'all')*6/256;
                tomoSlice_out = obj.getTomogram(test,'physUnitOutTF',false);
            end

            function [entOut,tomoSlice_out] = axisInEntropyOut(axisOffset_in)
                [tomoSlice_out,histStep] = axisInTomoSliceOut(axisOffset_in);

                tomoSlice_out = round( tomoSlice_out./histStep );
                tomoSlice_out = uint8(tomoSlice_out - min(tomoSlice_out(:)));
                entOut = entropy( tomoSlice_out(tomoMask) );
            end

        end
      
        function [bFieldSelectedY, ringRemovePrm] = autoRingPatRmvScan(obj, bFieldSelectedY,varargin )
            %% inputParsing
            p = inputParser();
            addParameter(p, 'filterName', 'wiener') % large waveletN removes more stripes, may lower resolution.
            addParameter(p, 'waveletNLib', 0:4, @isvector) % large waveletN removes more stripes, may lower resolution.
            addParameter(p, 'wname', 'db25', @ischar) % fix for default. Not significantly different results.
            addParameter(p, 'sigma', 2.4, @isscalar)    % fix for default. higher value makes it smoother
            addParameter(p, 'figTF', false, @islogical)

            parse(p, varargin{:});
            filterName    = p.Results.filterName;     % pixel; empty input --> auto finder
            waveletNLib   = p.Results.waveletNLib;     % pixel; empty input --> auto finder
            wname         = p.Results.wname;     % pixel; empty input --> auto finder
            sigma         = p.Results.sigma;
            figTF         = p.Results.figTF;

            %% run
            imgEnt = zeros(length(waveletNLib),1);
            for ii = 1:length(waveletNLib)
                % ring pattern removal
                test    = obj.removeStripeOfSino(bFieldSelectedY,'wavelet-fft',waveletNLib(ii),wname,sigma);

                % filter
                test = obj.filterSinogram(test,'real',filterName);

                % outputs
                histStep = std(abs(test),0,'all')*6/256;
                testImg = obj.getTomogram(test,'physUnitOutTF',false);

                % entropy calc.
                testImg = round(testImg./histStep);
                testImg = uint8(testImg - min(testImg(:)));
                imgEnt(ii) = entropy(testImg);

                [~, minInd] = min(imgEnt(1:ii));
                if ii == minInd
                    minImg   = testImg;
                    minSino  = test;
                    waveletN = waveletNLib(minInd);
                end

                if figTF
                    subplot(131),tomoHandle.show(testImg);axis image; title(waveletN);colorbar;
                    subplot(132),tomoHandle.show(test); axis image; title(wname); colorbar;
                    subplot(133),plot(waveletNLib(1:ii),imgEnt(1:ii),waveletNLib(minInd),imgEnt(minInd),'r*')
                    drawnow;
                    pause;
                end

            end
            figure(15151),
            subplot(131),tomoHandle.show(minImg);axis image; colorbar;
            title(['waveletN: ',num2str(waveletN),', sigma: ', num2str(sigma)]);
            subplot(132),tomoHandle.show(minSino); axis image; title(wname); colorbar;
            subplot(133),plot(waveletNLib,imgEnt,waveletNLib(minInd),imgEnt(minInd),'r*')
            drawnow;

            ringRemovePrm = struct;
            ringRemovePrm.waveletN    = waveletN;
            ringRemovePrm.wname       = wname;
            ringRemovePrm.sigma       = sigma;
            fprintf('*** Auto ringRemove result =  waveletN: %d, wname: %s, sigma: %d\n',waveletN,wname,sigma);


            %%% manual check
            while true
                test    = obj.removeStripeOfSino(bFieldSelectedY,'wavelet-fft', waveletN ,wname,sigma);
                test    = obj.filterSinogram(test,'real',filterName);
                testImg = obj.getTomogram(test,'physUnitOutTF',false);

                subplot(131),tomoHandle.show(testImg);axis image;
                title(['waveletN: ',num2str(waveletN),', sigma: ', num2str(sigma)]);
                colorbar;
                subplot(132),tomoHandle.show(test); axis image; title(wname); colorbar;
                drawnow;

                QQ = input('Type empty to stop. type # = waveletN (>=0); db##= wname; s##= sigma : ','s');

                % break if empty
                if isempty(QQ)
                    break;
                end


                if ~isnan(str2double(QQ)) % waveletN control
                    waveletN = str2double(QQ);
                elseif startsWith(QQ,'db','IgnoreCase',true) % wname control
                    wname    = QQ;
                elseif startsWith(QQ,'s','IgnoreCase',true)  % sigma control
                    sigma    = str2double(QQ(2:end));
                end
            end

            ringRemovePrm.waveletN    = waveletN;
            ringRemovePrm.wname       = wname;
            ringRemovePrm.sigma       = sigma;
            close(gcf);

        end

        function bField  = prepareSinogram( obj,bField,filterName )
            if nargin < 3
                filterName = 'wiener';
            end

            axisOffsetLocal    = obj.axisOffset;
            ringRemovePrmLocal = obj.ringRemovePrm;
            %             outType            = p.Results.outType;

            %%% yind selection
            if isempty(ringRemovePrmLocal) || isempty(axisOffsetLocal)
                th = tomoHandle.show(real(bField));
                drawnow;

                yy = size(bField,1);
                str = sprintf('*** Put y-axis index [center = %d]: ',mcoor(yy));
                yselect = input(str);
                if isempty(yselect)
                    yselect = mcoor(yy);
                end
                close(th)

                bFieldSelectedY = bField(yselect,:,:);
                %     tomoHandle.show(bFieldSelectedY);

            end

            %% auto rotation axis determination
            if isempty(axisOffsetLocal)
                axisOffsetLocal = obj.autoRotAxisScan( bFieldSelectedY,'filterName',filterName);
                % axisOffset = obj.autoRotAxisScan( bFieldSelectedY,histStep );
            end

            %% auto ring pattern removal
            if isempty(ringRemovePrmLocal)
                % axis shift
                bFieldSelectedY = obj.rotAxisShift(bFieldSelectedY, axisOffsetLocal, obj.FOV.bgAddXFOV);

                % merge if KK
                bFieldSelectedY = obj.KKAngleMerge(bFieldSelectedY);

                [~, ringRemovePrmLocal] = obj.autoRingPatRmvScan( bFieldSelectedY, 'filterName',filterName );
            end

            %% outType check
            if ismember( obj.fieldReconMethod, {'proj','projection','PC'})
                outType = 'real';
            else
                outType = 'complex';
            end

            %% bField rdy
            % rotation center.
            bField = obj.rotAxisShift(bField, axisOffsetLocal, obj.FOV.bgAddXFOV);

            % merge if KK -> log
            bField = obj.KKAngleMerge(bField);            
            
            % ring pattern rmv.
            bField = obj.removeStripeOfSino(bField,'wavelet-fft', ...
                ringRemovePrmLocal.waveletN, ringRemovePrmLocal.wname, ringRemovePrmLocal.sigma);
            
            % filter
            bField = obj.filterSinogram(bField,outType,filterName);                      % filter

            %  figure,   tomoHandle.show(real(bField),99); axis image;
            
            %             test = obj.removeStripeOfSino(bField(159,:,:), ...
            %                 3, ringRemovePrm.wname, 1); % ring pattern rmv.

            %             tomoHandle.show(-real(obj.getTomogram(test)),99.99); axis image;% delta
            %             tomoHandle.show(real(bField(159,:,:)),99.99); axis image;% delta
            %             tomoHandle.show(real(bField),99.99); axis image;% delta
            %

            %% out
            obj.axisOffset    = axisOffsetLocal;
            obj.ringRemovePrm = ringRemovePrmLocal;

        end

        function tomoOut = getTomogram(obj,RytovF,varargin)
            [yy,xx,zz] = size(RytovF);

            p = inputParser();
            addParameter(p, 'physUnitOutTF', true, @islogical);

            parse(p, varargin{:});
            physUnitOutTF    = p.Results.physUnitOutTF;


            %% duplicated angle consideration
            tomoAngles =  obj.rotAnglelib;
            tomoAngles( zz+1 : end ) = [];   % erase invalid angles (possibly merged before)
            tomoAngles = cast(tomoAngles,'like',real(RytovF));

            flipInd = find( tomoAngles>=180 & tomoAngles<360 ); % 0 and 180 will be overlapped during iradon
            if ~isempty(flipInd)
                counterPartInd  = flipInd - flipInd(1) + 1;
                RytovF(:,:,[counterPartInd, flipInd]) = 1/2 * RytovF(:,:,[counterPartInd, flipInd]); % 1/2 for overlapped angles.
            end
            NAngles = length(tomoAngles) - length(flipInd); % effectuve number of angles considering overlapping
            %   tomoHandle.show(imag(RytovF));

            %% Tomo recon
            switch obj.tomoReconMethod
                case 'FBP'
                    assert(mod(xx,2) == 1,'size(RytovF,2) must be odd')
                    realTF = isreal(RytovF);

                    %%% gpu check
                    if obj.gpuInd
                        if isgpuarray(RytovF)
                            gpu=gpuDevice;
                        else
                            gpu = gpuDevice(obj.gpuInd);
                            RytovF = gpuArray(RytovF);
                        end
                    end
                    tomoOut  = zeros(yy,xx,xx,'like',RytovF);
                    RytovF = RytovF * length(tomoAngles) / NAngles; % in order to compensate 1/length(tomoAngles) in iradon.

                    if obj.gpuInd
                        %%% chunk def
                        if realTF
                            bytesPerAngle = (xx*xx) * 4; % 4 bytes for single
                        else
                            bytesPerAngle = (xx*xx) * 8; % 8 bytes for complex single
                        end

                        maxChunkSize = floor( gpu.AvailableMemory / bytesPerAngle / 10 ); % last factor is a heuristic safety value.
                        iterN = ceil( yy / maxChunkSize);
                        chunks = round( linspace(1, yy+1, iterN+1) );

                        %%% run
                        for ii = 1:iterN
                            thisChunk = chunks(ii) : chunks(ii+1)-1;
                            parfor yi = thisChunk
                                RytovFy = squeeze(RytovF(yi,:,:));

                                %%% tomoGeneration
                                if realTF
                                    tomoOut(yi,:,:) = iradon_vBPpart(RytovFy,tomoAngles,'linear',xx);
                                    %                                     scatteringPotential(yi,:,:) = iradon(RytovFy,tomoAngles,'none',1,xx);
                                else
                                    tomoOut(yi,:,:) = ...
                                        iradon_vBPpart(real(RytovFy), tomoAngles,'linear',xx)...
                                        +1i* iradon_vBPpart(imag(RytovFy), tomoAngles,'linear',xx);
                                end
                            end
                            %                 fprintf('Tomo recon. finished: %d / %d \n', ii,iterN);
                        end

                    else % no GPU

                        vec = (1:xx) - mcoor(xx); %% this is identical with floor(xx/2) + 1 for odd.
                        [XX,YY] = meshgrid(vec,-vec);

                        parfor yi = 1:yy
                            RytovFy = squeeze(RytovF(yi,:,:));

                            %%% tomoGeneration
                            if realTF
                                tomoOut(yi,:,:) = iradon_vBPpart(RytovFy,tomoAngles,'linear',xx,XX,YY);
                            else
                                tomoOut(yi,:,:) = ...
                                    iradon_vBPpart(real(RytovFy), tomoAngles,'linear',1,xx,XX,YY)...
                                    +1i* iradon_vBPpart(imag(RytovFy), tomoAngles,'linear',1,xx,XX,YY);
                            end
                        end
                    end

                otherwise
                    error('NEED TO BUILD')
            end


            %% masking
            if isempty(obj.axisOffset)
                axisOffsetLocal=0;
            else
                axisOffsetLocal = ceil( abs(obj.axisOffset) );
            end
            xR = (xx-1)/2 - axisOffsetLocal; %% xx is odd
            tomoMask = ~mk_ellipse(xR,xR,xx,xx); %% xx is odd

            %%% make bg = 0
            bgR = xR - obj.FOV.bgAddXFOV;
            bgMask  = mk_ellipse(bgR,bgR,xx,xx) & tomoMask;
            %             figure,imagesc(bgMask),axis image
            tomoOut = reshape(tomoOut, yy, []);
            bgVal   = mean( tomoOut(:,bgMask) , 2);
            tomoOut = tomoOut - bgVal; % make bg zero for every y-slice
            tomoOut = reshape(tomoOut, yy, xx, xx);

            %%% masking
            tomoOut = tomoOut.* reshape(tomoMask,1,xx,xx); % make outside of FOV zero
            % figure,tomoHandle.show(-real( tomoOut ),99.99);
            % figure,tomoHandle.show(-real( tomoOut .* reshape(tomoMask,1,xx,xx) ),99.99);

            %%% At this point, the unit of tomoOut is (pixel)^-1
            %%% = pi/wl * (n^2 -1), n = 1 -delta +i*beta;
            %%% = 2*pi/wl*( -delta +1i*beta ) (if delta, beta << 1)

            if physUnitOutTF
                wl     = Etowl(obj.setup.source.energy_eV/1000); % nm
                imgPix = obj.setup.imagePixelSize;  % nm
                tomoOut = tomoOut * (wl/2/pi /imgPix); % unitless
            end

            %%% At this point, the unit of tomoOut is unitless
            %%% = 1/2 * (n^2 -1), n = 1 -delta +i*beta;
            %%% = -delta +1i*beta; (if delta, beta << 1)

            % tomoHandle.show(real(tomoOut));
        end

        function [sino, groundTruth, positionOffs] = phantomMaker(obj,varargin)
            %%% 3D RI = Modified Shepp-Logan phantom multiplied by (complex RI - 1) of given material for given energy.

            %% def.
            imgPix      = obj.setup.imagePixelSize; % nm
            photonE     = obj.setup.source.energy_eV/1000; % keV

            %% parse
            p = inputParser();
            addOptional(p,  'samSize_nm', [], @isscalar); % nm
            addParameter(p, 'material', 'Si', @(v) ismember(v , {'Au','Si','Ni','Al','Diamond'}));
            addParameter(p, 'padRatio', 3, @isscalar);
            addParameter(p, 'maxPhotonCount', 1000, @isscalar);

            addParameter(p, 'instableStd', [0,0], @isvector);   % pixel
            addParameter(p, 'sampleOffCenter', 0, @isscalar);   % pixel
            addParameter(p, 'axisOffset', 0, @isscalar);   % pixel

            addParameter(p, 'detailParam', struct([]), @isstruct);

            parse(p, varargin{:});
            samSize_nm      = p.Results.samSize_nm;
            material        = p.Results.material;
            padRatio        = p.Results.padRatio;
            instableStd     = p.Results.instableStd;
            sampleOffCenter = p.Results.sampleOffCenter;
            axisOffsetLocal      = p.Results.axisOffset;
            maxPhotonCount  = p.Results.maxPhotonCount;
            detailParam     = p.Results.detailParam;

            v = [padRatio,maxPhotonCount];
            mustBeInteger(v);
            mustBePositive(v);

            %% samSize def.
            if isempty(samSize_nm)
                samSize    = 127;
                samSize_nm = samSize * imgPix;
            else
                samSize    = 2*round( samSize_nm/imgPix/2 ) + 1; % odd
            end

            %% rotAnglelib def.
            pracRes = max( [obj.tomoRes; imgPix*ones(1,2)] );

            if isempty( obj.rotAnglelib)
                Nangle = ceil( pi /2 * samSize_nm / pracRes(2) );
                obj.rotAnglelib = linspace(0,180,Nangle); % last angle duplication mimicking the experiments.

                if obj.rotAngleMax == 360
                    obj.rotAnglelib = [ obj.rotAnglelib,  obj.rotAnglelib(2:end) + 180];
                    % even for rotAngleMax == 360, must include deg = 180 to satisfy symmetry
                    %        Nangle = length(rotAnglelib);
                end
            end

            %% gpuSetting
            if obj.gpuInd
                gpuDevice(obj.gpuInd);
                varType = zeros(1,'single','gpuArray');
            else
                varType = zeros(1,'single');
            end

            %% phantom gen
            cRI = XrayComplexRI(photonE, material)-1;
            wl  = Etowl(photonE);
            phase3D = abs(phantom3dAniso(varType,'Modified Shepp-Logan',samSize)) * cRI * (2*pi/wl * imgPix); % complex
            phase3D = permute(phase3D, [3,2,1]); % make z axis parallel to rotation axis (y, vertical);
            % figure, tomoHandle.show(real(phase3D));  colorbar

            %%% axisOffset application
            assert(round(sampleOffCenter) == sampleOffCenter, 'axisOffset for a phantom must be an integer')
            [yy,xx,zz] = size(phase3D);
            padXsize = round( sampleOffCenter/2 );
            phase3D = cat(2, zeros(yy,padXsize,zz), phase3D, zeros(yy,padXsize,zz));
            phase3D = circshift(phase3D, [0,sampleOffCenter,0]);


            %% window def.
            NAfrac  = 1/2*imgPix /obj.tomoRes(1);
            OCDfrac = 1/2*imgPix * (1/obj.tomoRes(2) - 1/obj.tomoRes(1));

            %% solution gen.
            %%% 3D window making
            vecfunc = @(v) ((1:v) - mcoor(v))/v;
            xvec = vecfunc(xx);
            yvec = vecfunc(yy);
            zvec = vecfunc(zz);

            [XX,YY,ZZ] = meshgrid(xvec,-yvec,zvec);

            rho = sqrt(XX.^2 + ZZ.^2);
            clear XX ZZ
            window3D = (rho-OCDfrac).^2 + YY.^2 <= NAfrac^2;

            %%% application
            window3D = ifftshift(window3D);
            groundTruth  = fftn(phase3D);
            groundTruth(~window3D) = 0;
            groundTruth  = ifftn(groundTruth);

            %%% release from the GPU
            groundTruth  = gather( groundTruth / (2*pi/wl * imgPix) );
            %   tomoHandle.show(window3D,99.99);  colorbar;
            %   figure,tomoHandle.show(real(groundTruth),99.99);  colorbar;

            %% sino gen
            rotAnglelib_forRadon = mod(obj.rotAnglelib,360);
            for yi = 1:samSize
                thisSlice = squeeze( phase3D(yi,:,:) );

                if yi == 1
                    sino = radon(real(thisSlice) ,rotAnglelib_forRadon)...
                        + 1i * radon(imag(thisSlice) ,rotAnglelib_forRadon);
                    xx   = size(sino,1);
                    sino = reshape(sino,1,xx,[]);
                    sino(samSize,:,:) = 0; % initialization

                else
                    sino(yi,:,:) = reshape(...
                        radon(real(thisSlice) ,rotAnglelib_forRadon)...
                        +1i * radon(imag(thisSlice) ,rotAnglelib_forRadon) ,1,xx,[]);
                end
            end
            sino = exp(1i*sino); % E-field conversion.
            clear phase3D
            %   tomoHandle.show(abs(sino),99.99);  colorbar;
            %   tomoHandle.show(angle(sino),99.99); colorbar;


            %% padding : done to reduce the boundary effects during the windowing.
            [yy0,xx0] = size(sino,[1,2]);
            sino = mpad(sino - 1, padRatio*[yy0,xx0]) + 1; % padding with 1

            [yy,xx,zz] = size(sino);
            xvec = vecfunc(xx);
            yvec = vecfunc(yy);
            [XX,YY] = meshgrid(xvec,-yvec);

            %% 2D windowing
            sino = fft2( sino );
            sino = reshape(sino,[],zz);

            switch obj.fieldReconMethod
                case 'KK'

                    %%% offAngle settting: can be inaccurate in the experiments
                    offAngle = obj.setup.afterSample.zonePlate.offAxis.angleDeg;      % nm

                    if isfield(detailParam,'offAngle') %overide
                        if ~isempty(detailParam.offAngle)
                            offAngle = detailParam.offAngle;
                        end
                    end

                    %%% maskOffset settting: can be inaccurate in the experiments
                    maskOffset = 0;
                    if isfield(detailParam,'maskOffset') %overide
                        if ~isempty(detailParam.maskOffset)
                            maskOffset = detailParam.maskOffset / max(yy0,xx0); % 1pixel for 1/FOV shift (before pad)
                        end
                    end
                    mustBeNonnegative(maskOffset);

                    %%% window2D def.
                    yc= sind(offAngle)*OCDfrac;
                    xc= cosd(offAngle)*OCDfrac;

                    window2D = (XX-xc).^2 + (YY-yc).^2 <= NAfrac^2; % NAmask
                    window2D = (cosd(offAngle)*XX + sind(offAngle)*YY >= -maskOffset) & window2D;
                    %                     imagesc(window2D)

                otherwise
                    window2D = XX.^2 + YY.^2 <= NAfrac^2; % NAmask
            end
            % imagesc(window2D); axis image;
            window2D = ifftshift(window2D);
            sino(~window2D,:) = 0;

            %%% phase plate of PC
            if strcmpi(obj.fieldReconMethod, 'PC')
                phDia     = obj.setup.afterSample.phasePlate.diameter;
                ZPdia     = obj.setup.afterSample.zonePlate.diameter;
                ZPdrn     = obj.setup.afterSample.zonePlate.outermostWidth;

                phR = 1/2 * phDia/ZPdia * imgPix/ZPdrn ;
                phMask = XX.^2 + YY.^2 <= phR^2; % NAmask
                phMask = ifftshift(phMask);

                sino(phMask,:) = 1i*sino(phMask,:); % pi/2 shift
                % imagesc(phMask); axis image;
            end

            sino = reshape(sino,yy,xx,zz);
            sino = ifft2( sino );
            % figure,tomoHandle.show(abs(sino).^2,99.99); colorbar;
            % tomoHandle.show(log10(abs(fftshift(fft2(sino)))),99.99); colorbar;


            %% sample position off definition
            positionOffs = zeros( length(obj.rotAnglelib) ,2);

            %%% axisOffset
            positionOffs(:,2) = positionOffs(:,2) + axisOffsetLocal;
            if axisOffsetLocal ~= 0
                sino = circshift_KR( sino, [0, axisOffsetLocal, 0] );
            end

            %%% instablity
            instableStd = reshape(instableStd,1,[]); % make row vector
            insta_positionOffs = instableStd .* randn(zz,2);
            positionOffs = positionOffs + insta_positionOffs;

            if any( instableStd ~= 0 )
                for zi = 1:zz
                    sino(:,:,zi) = circshift_KR( sino(:,:,zi), [insta_positionOffs(zi,1),insta_positionOffs(zi,2),0] );
                end
            end


            %% measurement: modulus square, crop, noise addition
            sino = mcrop(sino,[yy0,xx0]);
            sino = abs(sino).^2;

            mFactor = maxPhotonCount / max(sino(:));
            sino = uint16( sino * mFactor );
            sino = imnoise(sino,'poisson');
            % figure,tomoHandle.show(real(sino),99.99); colorbar;
            sino = gather( single(sino) / mFactor );

            %% FOV def.
            obj.FOV = struct;
            obj.FOV.xlim = [1,xx0];
            obj.FOV.bgAddXFOV = floor( xx0 * 0.03 );
            obj.FOV.ylim = [1,yy0];
            obj.FOV.size = [yy0,xx0];
            obj.FOV.originalSize = [yy0,xx0];
            obj.positionOffs = 0;

        end

        function saveStr = save(obj,varargin)
            %%% parsing
            p = inputParser();
            addOptional(p,'sinogram',[])
            addOptional(p,'scatteringPotential',[])

            %             addParameter(p,'outType','', @(a) ischar(a) || isstring(a) )
            %             addParameter(p,'axisOffset','', @isscalar )
            %             addParameter(p,'ringRemovePrm',struct([]), @isstruct )
            addParameter(p,'saveTag','', @(a) ischar(a) || isstring(a) )

            p.parse(varargin{:});
            sinogram = gather( p.Results.sinogram );
            scatteringPotential = gather( p.Results.scatteringPotential );

            %             outType = p.Results.outType;
            %             axisOffset = p.Results.axisOffset;
            %             ringRemovePrm = p.Results.ringRemovePrm;

            saveTag    = p.Results.saveTag;

            %%% save string
            fprintf('Saving ... '); tic;
            topDir   = obj.samFiles(1).folder;
            if isempty(saveTag)
                saveFileName = ['tomoHadleOut_v',num2str(obj.version),'.mat'];
            else
                saveFileName = ['tomoHadleOut_v',num2str(obj.version),'_',saveTag,'.mat'];
            end
            saveStr = fullfile(topDir,saveFileName);

            %%% save
            saveStruct.tH = obj;

            %             saveStruct.outType = outType;
            %             saveStruct.axisOffset = axisOffset;
            %             saveStruct.ringRemovePrm = ringRemovePrm;

            saveStruct.sinogram = sinogram;
            saveStruct.scatteringPotential = scatteringPotential;

            save(saveStr,'-struct','saveStruct','-v7.3');
            fprintf('done .. %.1f sec\n',toc)
        end

    end

end

%%% used functions
function out = mfft2(in,varargin)
out = fftshift(fft2(ifftshift( in ),varargin{:}));
end
% function out = mifft2(in,varargin)
%     out = fftshift(ifft2(ifftshift( in ),varargin{:}));
% end

function [bout,fileDir,outType] = matParser(fileDir,showType)
mustBeMember(showType,{'auto','sino','delta','beta'})

%%% empty b
if isempty(fileDir)
    [readFile, readPath] = uigetfile({'*.mat', 'MATLAB Data (*.mat)'}, 'Select a tomoHandle savefile');
    fileDir = fullfile(readPath, readFile);
end

%%% show
if ~exist(fileDir,'file')
    error('invalid input')
end

matObj = matfile(fileDir);
sinoTF = all( size(matObj,'sinogram') ~= 0 );
tomoTF = all( size(matObj,'scatteringPotential') ~= 0 );

outType = '';
switch showType
    case 'auto'
        if tomoTF % if tomo is availalbe, view delta.
            bout = -real( matObj.scatteringPotential );
            outType = 'delta';
        elseif sinoTF
            bout = matObj.sinogram;
            outType = 'sino';
        else
            error('empty matfile');
        end
    case 'sino'
        if sinoTF
            bout = matObj.sinogram;
            outType = 'sino';
        else
            error('empty matfile');
        end
    case 'delta' % >0
        if tomoTF
            bout = -real( matObj.scatteringPotential );
            outType = 'delta';
        else
            error('empty matfile');
        end
    case 'beta'  % >0
        if tomoTF
            bout = imag( matObj.scatteringPotential );
            outType = 'beta';
        else
            error('empty matfile');
        end
end
end

function img = imreadAndRotation(imgDir,rotAngle, darkImg)
if isempty(darkImg)
    img = single(imread(imgDir));
else
    img = single(imread(imgDir)) - darkImg;
end
img = imrotate(img,rotAngle,'bilinear');
end