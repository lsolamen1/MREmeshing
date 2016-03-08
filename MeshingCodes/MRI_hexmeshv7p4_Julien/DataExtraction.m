function [A,vox,mridim,MagIm,P,freqHz,msk,mask]=DataExtraction(usedef,default)
inputtypedef='D';
inputtype=input(['MR Data Source <D = Dartmouth, I = Illinois> (default is ' inputtypedef '):  >>'],'s');
if isempty(inputtype)
    inputtype=inputtypedef;
end
if inputtype=='D' % Dartmouth Data
    % Define name of File data.
    
    % MRE data files:
    if(~usedef)
        MRIfile=input(['MRI motion file name? (default is ' default.mridef '):  >>'],'s');
    end
    if ~exist('MRIfile','var')||isempty(MRIfile)
        MRIfile=default.mridef;
    end
    
    % Header Data file
    if(~usedef)
        Hdrfile=input(['Header file name? (default is ' default.Hdrdef '):  >>'],'s');
    end
    if ~exist('Hdrfile','var')||isempty(Hdrfile)
        Hdrfile=default.Hdrdef;
    end
    
    % Mask
    if(~usedef)
        msk=input(['Name of mask to use for mesh (Default ' default.mskdef '):  >>'],'s');
    end
    if ~exist('msk','var')||isempty(msk)
        msk=default.mskdef;
    end 
    
    % Load files
    load(MRIfile);
    if ~exist('A');
        if exist('A1');
            num=input('Multi-Frequency Data.  Enter number for frequency set: ','s');
            eval(['A = A',num,';'])
            eval(['P = P',num,';'])
        end
    end
    load(Hdrfile);
    load(msk);
    mridim=size(MagIm);
    vox=DirIndex(4,1:3).*1e-3;
    
elseif inputtype=='I' %Illinois Data
    % MRE data files:
    d=dir('*.mat');
    for ii=1:length(d)
        disp(['File ' int2str(ii) ' :: ' d(ii).name])
    end
    n=input('Number of file to use (default = 1)  >> ');
    if(isempty(n))
        n=1;
    end
    MRIfile=d(n).name;
    load(MRIfile);
    mskdef='Mask.mat';
    if ~exist('msk','var')||isempty(msk)
        msk=mskdef;
    end
    
    mridim=size(t2stack);
    vox=([mreParams.FOVx,mreParams.FOVy,mreParams.FOVz]./[mreParams.nx,mreParams.ny,mreParams.nz])*1e-3;
    freqHz=mreParams.freq;
    MagIm=t2stack;
    A=zeros([size(MagIm) 3]);
    P=zeros([size(MagIm) 3]);
    A(:,:,:,1)=abs(Xmotion);
    A(:,:,:,2)=abs(Ymotion);
    A(:,:,:,3)=abs(Zmotion);
    P(:,:,:,1)=angle(Xmotion);
    P(:,:,:,2)=angle(Ymotion);
    P(:,:,:,3)=angle(Zmotion);
    clear Xmotion Ymotion Zmotion
else
    error('No Suitable Input Data Format Specified');
end %End Data Type Condition
end
