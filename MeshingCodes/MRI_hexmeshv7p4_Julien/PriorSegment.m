function [regs,noreg,dim]=PriorSegment(usedef,MagIm)
regdef='none';
if(~usedef)
    d=dir('*.mat');
    for ii=1:length(d)
        disp(['File ' int2str(ii) ' :: ' d(ii).name])
    end
    n=input(['Name of Soft Prior Segmentation File (Default ' regdef '):  >>']);
    if(~isempty(n))
        regdef=d(n).name;
    end
end
if ~exist('regf','var')||isempty(regf)
    regf=regdef;
end

% Load segmentation file if supplied
dim=size(MagIm);
if(strcmp('none',regf))
    disp('Using no soft prior segmentation')
    regs=zeros(dim);
    noreg=true;
else
    disp(['Loading segmentation file ' regf])
    load(regf)
    noreg=false;
end
end