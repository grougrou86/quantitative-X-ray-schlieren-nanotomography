function imgout = mcenter(img, centerPosition)

inputSize = size(img);
inputDim  = length(inputSize);
centerDim = length(centerPosition);

if inputDim < centerDim
    error('Dimensions of input images must be larger than given dimesion')
else
    conservingDim = centerDim+1:inputDim;
end

if  any(centerPosition > inputSize(1:centerDim))
    error('CenterPosition must be smaller than the image size')
end

%%%
displacement = floor(inputSize(1:centerDim)/2)+1 - centerPosition;
displacement(conservingDim) = 0;

imgout  = circshift_KR(img, displacement);
