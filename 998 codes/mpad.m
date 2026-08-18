function output = mpad(input,finalSize)

% if ndims(input) ~= length(finalSize)
%     error( 'The dimension of matrix must be equal to the length of vector' )
% end

if ~isvector(finalSize)
    error( 'Input size must be a vector' )
end

currentDim = length(finalSize);
originalSize = size(input);
originalSize = originalSize(1:currentDim);

if  any(finalSize < originalSize)
    error( 'The final dimension must be larger than initial dimension' )
end

originalSize = cast(originalSize,'like',finalSize);
mvec = floor(finalSize/2)+1;
mMsz = floor(originalSize/2)+1;
output = circshift(padarray(input,finalSize - originalSize,'post'),mvec-mMsz);
