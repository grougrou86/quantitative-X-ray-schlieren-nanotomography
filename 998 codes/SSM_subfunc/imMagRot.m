function output = imMagRot(input,mag,rotAngle,varargin)
% rotation (in deg) in a counter-clockwise direction (as in imrotate).
p = inputParser();

addOptional(p, 'outSize',[], @(s) isvector(s) || isempty(s));    % y,x
addOptional(p, 'interpMethod','linear',@(s) ismember(s,{'linear','nearest','cubic','makima','spline'}));
addParameter(p, 'magMethod','',@(s) ismember(s,{'legacy',''}));

parse(p, varargin{:});
outSize = p.Results.outSize;
magMethod = p.Results.magMethod;
interpMethod = p.Results.interpMethod;

%% make mag vector
if isscalar(mag)
    mag = [mag,mag];
end

%%
if any(mag <= 0)
    error('Magnification factor must be > 0');
end
%%% trivial solution
if all(mag == 1) && (mod(rotAngle,360) == 0)
    output = input;
    if ~isempty(outSize)
        output = mcrop(input,outSize);
    end
    return
end


%%    
if  isempty(magMethod) && any(mag ~= 1) % default magMethod    
    if any(mag < 1)  
        integerOSratio = ceil( 1./mag );
        mag0 = mag;
        mag = mag0 .* integerOSratio;    % make mag >= 1     
    else 
        integerOSratio = 1;
    end
        
    inSize = size(input);    
    padSize = floor(mag.*inSize);
    magForInterp =  mag ./ ( padSize./inSize );    % mag >= 1 but very close to 1.
    
    Finput = fftshift(fft2(ifftshift(input)));
    Finput = mpad(Finput, padSize);
    inputForInterp = fftshift(ifft2(ifftshift(Finput))) * prod(padSize)/prod(inSize); %intensity norm term.
else
    inputForInterp = input;
    magForInterp   = mag;
    integerOSratio = 1;
end

[yy,xx] = size(inputForInterp);

xvec = (1:1:xx) - mcoor(xx);
yvec = (1:1:yy) - mcoor(yy);
[XX,YY] = meshgrid(xvec,yvec);

if isempty(outSize)
    outSize = [yy,xx];
end

outSizeForInterp  = outSize .* integerOSratio;
yvec2 = (1:1:outSizeForInterp(1)) - mcoor( outSizeForInterp(1) );
xvec2 = (1:1:outSizeForInterp(2)) - mcoor( outSizeForInterp(2) );
[XX2,YY2] = meshgrid(xvec2,yvec2);

XX2 = double(XX2);
YY2 = double(YY2);
XX2rot = ( +XX2*cosd(rotAngle) - YY2*sind(rotAngle) ) /magForInterp(2);
YY2rot = ( +XX2*sind(rotAngle) + YY2*cosd(rotAngle) ) /magForInterp(1);

output = interp2(XX,YY,inputForInterp,XX2rot,YY2rot,interpMethod,0);



% Foutput = fftshift(fft2(ifftshift(output)));
% imagesc(real(output)),axis image
% subplot(121),imagesc(log10(abs(Foutput)),[3,6]),axis image
% subplot(122),imagesc(log10(abs(Finput)),[3,6]),axis image

output = downsample2d(output,integerOSratio,0,'centeredTF',1);
 % imagesc(output),axis image
end
