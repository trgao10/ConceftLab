function [outMat, lowQVal, highQVal] = qclamp(inMat, lowQ, highQ)
%QCLAMP Summary of this function goes here
%   Detailed explanation goes here

if nargin < 3
    highQ = 1-lowQ;
end

lowQVal = quantile(inMat(:), lowQ);
highQVal = quantile(inMat(:), highQ);

outMat = inMat;
outMat(inMat < lowQVal) = min(inMat(:));
outMat(inMat > highQVal) = max(inMat(:));
% outMat(inMat < lowQVal) = lowQVal;
% outMat(inMat > highQVal) = highQVal;

end

