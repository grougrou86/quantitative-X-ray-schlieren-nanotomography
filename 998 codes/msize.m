function output = msize(input,outputSize)

inputSize = size(input);
inputDim = length(inputSize);
outputDim = length(outputSize);

if inputDim < outputDim
    error('Dimensions of input images must be larger than given dimesion')
       
elseif inputDim > outputDim    
    conservingDim = (outputDim+1):inputDim;
    outputSize(conservingDim) = inputSize(conservingDim);
end

cropInd = inputSize >= outputSize;

cropSize = inputSize;
cropSize(cropInd) = outputSize(cropInd);
cropOut = mcrop(input, cropSize);

padInd = inputSize < outputSize;
padSize = cropSize;
padSize(padInd) = outputSize(padInd);

output = mpad(cropOut, padSize);

