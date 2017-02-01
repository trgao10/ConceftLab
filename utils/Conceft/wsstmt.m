function [sstmt,f,mtcell,cwtcfscell,cwtfreqscell] = wsstmt(x,mt,varargin)
%WSSTMT Multitapered Wavelet Synchrosqueezed Transform
%
% Tingran Gao (trgao10@math.duke.edu)
% last modified: Jan 30, 2017
%

mtcell = cell(1,mt);
fcell = cell(1,mt);
vararginCell = cell(1,mt);
cwtcfscell = cell(1,mt);
cwtfreqscell = cell(1,mt);

wavparamsflag = find(strncmpi('WaveletParameters',varargin,1));

for k=1:mt
    vararginCell{k} = varargin;
    vararginCell{k}{wavparamsflag+1}.k = k-1;
end

parfor k=1:mt
    [mtcell{k}, fcell{k}, cwtcfscell{k}, cwtfreqscell{k}] = wsstgao(x,vararginCell{k}{:});
    mtcell{k} = mtcell{k}/sum(mtcell{k}(:));
end

sstmt = sum(cat(3,mtcell{:}),3) / mt;
f = fcell{1};

end

