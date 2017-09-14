function [outMat, lowQVal, highQVal] = qclamp(inMat, lowQ, highQ)
%QCLAMP Naive thresholding
%   simply set values in inMat beyond [lowQ,highQ] to 0
%
%   Tingran Gao (tingrangao@galton.uchicago.edu)
%   last modified: Sep 13, 2017
%

if nargin < 3
    highQ = 1-lowQ;
end

lowQVal = quantile(inMat(:), lowQ);
highQVal = quantile(inMat(:), highQ);

outMat = inMat;
% outMat(inMat < lowQVal) = min(inMat(:));
% outMat(inMat > highQVal) = max(inMat(:));
outMat(inMat < lowQVal) = lowQVal;
outMat(inMat > highQVal) = highQVal;

end

