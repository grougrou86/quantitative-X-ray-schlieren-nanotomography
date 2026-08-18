% cimg = mcrop(img,[yy,xx,zz])

function cimg = mcrop(img, outputSize)

inputSize = size(img);
inputDim = length(inputSize);
outputDim = length(outputSize);

if inputDim < outputDim
    error('Dimensions of input images must be larger than given dimesion')
       
elseif inputDim > outputDim    
    conservingDim = (outputDim+1):inputDim;
    outputSize(conservingDim) = inputSize(conservingDim);
end
if  any(outputSize > inputSize)
    error('Final image size must be smaller than initial image size')
end

if  inputDim >3
    error('This function covers upto 3-D')
end

dimvectors = zeros(2,inputDim);
for kk = 1:1:inputDim
    mp = floor(inputSize(kk)/2)+1;
    mside = round((outputSize(kk)-1)/2);
    dimvectors(:,kk) = [ mp-mside; mp+outputSize(kk)-mside-1 ];
end


switch inputDim
    case 1
        cimg = img(dimvectors(1,1):dimvectors(2,1));
    case 2
        cimg = img(dimvectors(1,1):dimvectors(2,1),...
            dimvectors(1,2):dimvectors(2,2));
    case 3
        cimg = img(dimvectors(1,1):dimvectors(2,1),...
            dimvectors(1,2):dimvectors(2,2),...
            dimvectors(1,3):dimvectors(2,3));
end
