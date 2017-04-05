function [sstmt,f,mtcell,cwtcfscell,phasetfcell,cwtfreqcell] = fsstmt(x,mt,varargin)
%FSSTMT Multitapered Short-Time Fourier Transform Synchrosqueezed Transform
%
% Tingran Gao (trgao10@math.duke.edu)
% last modified: Feb 15, 2017
%

mtcell = cell(1,mt);
fcell = cell(1,mt);
vararginCell = cell(1,mt);
cwtcfscell = cell(1,mt);
phasetfcell = cell(1,mt);
cwtfreqcell = cell(1,mt);

winparamsflag = find(strncmpi('WindowParameters',varargin,1));

for k=1:mt
    vararginCell{k} = varargin;
    vararginCell{k}{winparamsflag+1}.Order = k;
end

% parfor k=1:mt
for k=1:mt
    [mtcell{k}, fcell{k}, cwtcfscell{k}, phasetfcell{k}, cwtfreqcell{k}] = fsstgao(x,vararginCell{k}{:});
    mtcell{k} = mtcell{k}/sum(mtcell{k}(:));
end

sstmt = sum(cat(3,mtcell{:}),3) / mt;
f = fcell{1};

end

