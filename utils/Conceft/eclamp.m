function outMat = eclamp(inMat, ptg)
%ECLAMP Energy-based thresholding
%   set values in inMat that do not contribute to ptg of the total energy
%   to 0
%
%   Tingran Gao (tingrangao@galton.uchicago.edu)
%   last modified: Sep 13, 2017
%

if ptg < 0 || ptg > 1
    warning('percentage parameter should lie in [0,1]');
end

sortedElements = sort(abs(inMat(:)).^2,'descend');
energyThreshold = sum(sortedElements) * ptg;
thresholdIdx = find(cumsum(sortedElements) >= energyThreshold-1e-4, 1);
thresholdVal = sortedElements(thresholdIdx);

outMat = inMat;

if ~isempty(thresholdVal)
    sqInMat = abs(inMat).^2;
    mask = sqInMat < thresholdVal;
    outMat(mask) = 0;
end

end

