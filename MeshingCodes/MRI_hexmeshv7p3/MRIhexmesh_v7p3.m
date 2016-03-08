function MRIhexmesh_v7p3(default)
%MRIhexmesh_interp Generates required files for MREv7 from MRE data. 
%   Interpolates data to an approximate number of nodes per wavlength,
%   based on prior estiamtes of the stiffness. Builds Hex27 mesh manually,
%   i.e. does not require a 'template mesh'.
%   If called using MRIhexmesh('default'), uses the default values without
%   prompting for inputs. 
%   INPUTS:
%    default(optional): Set to 'default' to disable propmts for inputs. Any
%                       other value, or if not supplied will propmt for
%                       inputs. Using default values allows automated
%                       meshing.
%
%   IMPORTANT FEATURES
%
%   - DISPLACEMENT DATA IS SCALED!! I chose to scale everything so that the
%     average displacement amplitude is the same for all datasets
%     interpolated with this code. 
%   - Interpolation is an option. Cubic spline interpoaltion is used, a
%     custom function is used which does not use values outside of the mask
%     for the interpolated values.
%   - Different meshing strategies are possible: (MRE-Zonev7.05 now takes
%         zone dimensions rather than dividing the extents into an integer
%         number of zones.)
%      meshstrat=1 uses one FE node for each MR voxel, as has been done in
%                  the past.
%      meshstrat=2 allows specification of the number of FE nodes per
%                  wavelength, based on an estiamte of the shear modulus. A
%                  warning will be produced if the data is being
%                  interpolated to a resolutionlower than the data.
%      meshstrat=3 allows specification of a target resolution.
%      Both meshstrat=2 and meshstrat=3 make slight modifications to the
%      resolution to fit an integer number of 27 node elements across the
%      domain defined by the extents of the mask, to avoid wasting planes
%      of data near the edges.
%   - Different zone sizing strategies are possible: 
%      zonestrat=1 Matches a number of FE nodes per zone, as has been done 
%                  in the past.
%      meshstrat=2 matches a number of wavelengths per zone, based on an 
%                  estiamte of the shear modulus 
%      Both of these strategies are only approximate, the way that the
%      zoning proscess works (dividing the mesh extents into an integer
%      number of zones) makes it difficult to exactly match specified
%      values.
%   - Input data which may be useful when analyzing results or looking at
%     a reconstruction later on is saved as outstm.InterpLocations.mat,
%     and outstm.InterpData.mat outstm.meshinput.mat. The idea of saving
%     this data is to attach it to reconstruction results so it is clear
%     what parameters were used to produce it. 
%   
disp('MRE-Zone v7.3 Data Converter')
par=false; % If you have the matlab parallel processing toolbox set this to true, otherwise false.
nlabs=7; % Number of labs for matlab to use (max 7)
if(par)
    disp(['Using Parallelized version with ' int2str(nlabs) ' labs'])
else
    disp('Using non-Parallelized version - modify value of ''par'' to change')
end

if(nargin==0) % Default value is to prompt for inputs.
    default='no';
end
usedef=strcmp(default,'default');
if(usedef) 
    disp('Using Default values, no input prompts')
end

%% Get Inputs

% MRE data files:
mridef='MRE_3DMotionData.mat';
if(~usedef) 
    MRIfile=input(['MRI motion file name? (default is ' mridef '):  >>'],'s');
end
if ~exist('MRIfile','var')||isempty(MRIfile)
    MRIfile=mridef;
end

Hdrdef='HeaderData.mat';
if(~usedef) 
    Hdrfile=input(['Header file name? (default is ' Hdrdef '):  >>'],'s');
end
if ~exist('Hdrfile','var')||isempty(Hdrfile)
    Hdrfile=Hdrdef;
end

mskdef='Mask.mat';
if(~usedef) 
    msk=input(['Name of mask to use for mesh (Default ' mskdef '):  >>'],'s');
end 
if ~exist('msk','var')||isempty(msk)
    msk=mskdef;
end

if(exist('SPRregs.mat','file'))
    regdef='SPRregs.mat';
else
    regdef='none';
end
if(~usedef) 
    regf=input(['Name of Segmented Stack File (Default ' regdef '):  >>'],'s');
end 
if ~exist('regf','var')||isempty(regf)
    regf=regdef;
end


load(MRIfile);
load(Hdrfile);
load(msk);

if(ndims(MagIm)==4);
    MagIm=mean(MagIm,4);
end

vox=DirIndex(4,1:3).*1e-3;
rowswap=[0 1 0; 1 0 0 ; 0 0 1]; % This makes [1,:,:] = x, [:,1,:]=y, [:,:,1]=z
MPSto123=rowswap*DirIndex(1:3,4:6);

mridim=size(mask);

%Meshing Approach: one node per voxel, or interpolate to give nodes per wavelength
meshstratdef=1;
if(~usedef) 
    disp(['Meshing Strategies: 1=one node per MR voxel, 2=Target #nodes per wavelength, 3=Target Resolution'])
    meshstrat=input(['Meshing Strategy, (Default ' int2str(meshstratdef) ') >>']);
end
if ~exist('meshstrat','var')||isempty(meshstrat)
    meshstrat=meshstratdef;
end

%Zone Sizing Approach: nodes per zone or wavelngths per zone
zonestratdef=2;
if(~usedef) 
    disp(['Zone sizing Strategies: 1=fit nodes per zone, 2=fit wavelengths per zone, 3 = fixed size'])
    zonestrat=input(['Zone sizing Strategy, (Default ' int2str(zonestratdef) ') >>']);
end
if ~exist('zonestrat','var')||isempty(zonestrat)
    zonestrat=zonestratdef;
end

mudef=3300;
%mudef=20000;
if(meshstrat==2)||(zonestrat==2)
    % Estimated properties.    
    if(~usedef) 
        muest=input(['Estimate of Shear Modulus (Default ' num2str(mudef) ') >>']);
    end
    if ~exist('muest','var')||(isempty(muest))
        muest=mudef;
    end    

    rhoest=1000;
    wl=1/freqHz*sqrt(muest/rhoest);
end

if(meshstrat==2)
    % Desired nodes per wavelength.
    npwdef=12;
    if(~usedef) 
        npw=input(['Desired Nodes per wavelength: (Default ' num2str(npwdef) ') >>']);
    end
    if ~exist('npw','var')||(isempty(npw))
        npw=npwdef;
    end
end

if(zonestrat==1)
    % Nodes per zone
    npzdef=3500;
    if(~usedef) 
        npz=input(['Target Number of Nodes per zone (Default ' num2str(npzdef) ') >>']);
    end
    if ~exist('npz','var')||(isempty(npz))
        npz=npzdef;
    end     
elseif(zonestrat==2)
    %Wavelengths per zone    
    wlperzonedef=0.7;
    if(~usedef) 
        wlperzone=input(['Target Number of wavelengths per zone (Default ' num2str(wlperzonedef) ') >>']);
    end
    if ~exist('wlperzone','var')||(isempty(wlperzone))
        wlperzone=wlperzonedef;
    end 
elseif(zonestrat==3)
    zldef=1.95633e-02;
    if(~usedef) 
        zlength=input(['Zone size (m) (Default ' num2str(zldef) ') >>']);
    end
    if ~exist('wlperzone','var')||(isempty(zlength))
        zlength=zldef;
    end 
end
     
if(meshstrat==3)
    resdef=[0.002 0.002 0.002];
    if(~usedef) 
        targetres=input(['Target Resolution (Default ' num2str(resdef) ') >>']);
    end
    if ~exist('targetres','var')||(isempty(targetres))
        targetres=resdef;
    end
    
end

% output file stem
%outstmdef='MREv7mesh';
junk=pwd;
dirloc=find(junk=='/');
outstmdef=junk(dirloc(end)+1:end);
%outstmdef='MREv7mesh';

if(meshstrat==1)
    outstmdef=[outstmdef '_voxelmesh'];
elseif(meshstrat==2)
    outstmdef=[outstmdef '_npw' int2str(npw)];
elseif(meshstrat==3)
    outstmdef=[outstmdef '_res' sprintf('%2.1f-%2.1f-%2.1f',targetres.*1000)];
end
    
if(~usedef) 
    outstm=input(['output file stem (default is ' outstmdef '):  >>'],'s');
end
if ~exist('outstm','var')||isempty(outstm)
    outstm=outstmdef;
end

if(par)
    matlabpool('open',nlabs);
end

% Displacement Scaling - Many MRE regularization techniques are sensitive
% to the size of the displacements. Either the regularization weights can
% be altered for each case, or the displacements can be scaled so that they
% are always almost the same size.
dispscale = true; % Switch to turn on displacement scaling
dispscalar = 1e-3; % Average displacement amplitude is scaled to this size.
t0=tic;
%% Scale Data if appropriate
meanA = mean(A(repmat(mask,[1 1 1 3])==1));
disp(['Mean Displacement Amplitude all directions ' sprintf('%10.3e',meanA)])
if(dispscale)
    disp(['Scaling Displacements to an average size of ' sprintf('%10.3e',dispscalar) 'm'])
    A=A./meanA.*dispscalar;
end
%disp(['Check: mean(A(mask)) = ' num2str( mean(A(repmat(mask,[1 1 1 3])==1)) )])

% Generate Real and Imag Displacements and apply mask
Ur=A.*cos(P);
Ui=A.*sin(P);
clear A P


% % Apply a masked selective median filter to clear up outliers
% thresh=0.2; % Thereshold to apply median filter
% disp(['Median filter with threshold ' num2str(thresh)])
% for ii=1:3
%     Ur(:,:,:,ii)=selectivemedianfilter(Ur(:,:,:,ii),mask,thresh);
%     Ui(:,:,:,ii)=selectivemedianfilter(Ui(:,:,:,ii),mask,thresh);
% end

Ur(repmat(mask,[1 1 1 3])==0)=nan;
Ui(repmat(mask,[1 1 1 3])==0)=nan;

t1=tic;
dim=size(MagIm);

% Load segmentation file if supplied
if(strcmp('none',regf))
    disp('Using no soft prior segmentation')
    regs=zeros(dim);
    noreg=true;
else
    disp(['Loading segmentation file ' regf])
    load(regf)
    noreg=false;
end


if(meshstrat==1) % Each MR voxel is a FE node
    maskint=mask;
    Urint=Ur;
    clear Ur
    Uiint=Ui;
    clear Ui
    MagIm_int=MagIm;
    clear MagIm
    xin=0:vox(1):vox(1)*(dim(1)-1);
    yin=0:vox(2):vox(2)*(dim(2)-1);
    zin=0:vox(3):vox(3)*(dim(3)-1);
    xout=xin;
    yout=yin;
    zout=zin;   
    regs_int=regs;
    meshres=DirIndex(4,1:3)*1e-3;
elseif(meshstrat==2)||(meshstrat==3)  %  Process geometry and interpolate to get desired nodes per wavelength if meshstrat = 2.
    
    if(meshstrat)==2
        targetres(1:3) = wl/npw;
    elseif(meshstrat==3)
       % targetres already defined from inputs;
    end
    disp(['Target Resolution: ' num2str(targetres)])
        
    
    % Allow some deviation from the target resolution to fit an integer
    % number of Hex27 elements across the geometry. This will eliminate
    % wasting planes of data when there is an even number of slices which
    % occurs using a voxel based Hex27 meshing scheme.
    
    Imask=find(mask);
    [mi,mj,mk]=ind2sub(dim,Imask);
    ext(1)=(max(mi)-min(mi))*vox(1);
    ext(2)=(max(mj)-min(mj))*vox(2);
    ext(3)=(max(mk)-min(mk))*vox(3);
    
    % Number of targetres fitting across ext = ext(ii)/targetres. Want an
    % EVEN number of resolutions across the domain to fit 27 node elements
    % (Even number of resolutions = odd number of nodes) 
    meshres=zeros(1,3);
    for ii=1:3        
        nres=floor(ext(ii)/targetres(ii));
        if(mod(nres,2)==1)
            nres=nres+1;
        end
        meshres(ii) = ext(ii)/nres;
    end
    disp(['Actual New Resolution(m): ' num2str(meshres)])
    disp(['Data Resolution(m): ' num2str(vox)])
    
    % Warning if mesh is being interpolated to a lower resolution than the
    % data, this is not really a good idea.
    if(meshres(1)>vox(1)||meshres(2)>vox(2)||meshres(3)>vox(3))
        disp([' >> Warning: Mesh is being interpolated to a lower resolution than the data'])        
    end
    
    % Perform Interpolation
    xin=0:vox(1):vox(1)*(dim(1)-1);
    yin=0:vox(2):vox(2)*(dim(2)-1);
    zin=0:vox(3):vox(3)*(dim(3)-1);
    
    
    xout_tmp= (min(mi)-1)*vox(1):meshres(1):(max(mi)-1)*vox(1);
    yout_tmp= (min(mj)-1)*vox(2):meshres(2):(max(mj)-1)*vox(2);
    zout_tmp= (min(mk)-1)*vox(3):meshres(3):(max(mk)-1)*vox(3);
    
    % Include a buffer around the extents of the mask so that MagIm can be
    % interpolated, and the dataset extends past the boundaries of the
    % mask. e.g. for brain data, the buffer means the interpolated MagIm
    % shows the skull etc, not just the masked region of the brain.
    bufsiz=[4 4 0]; % Size of buffer in each direction (interpolated MR voxels) MUST BE EVEN
    if(sum(mod(bufsiz,2)>0))
        error('bufsiz must be even to fit hex27 elements efficiently inside mask extents!!')
    end
    
    xout=zeros(1,2*bufsiz(1)+length(xout_tmp));
    xout(bufsiz(1)+1:bufsiz(1)+length(xout_tmp))=xout_tmp;
    for ii=1:bufsiz(1)
        xout(bufsiz(1)+1-ii)=xout(bufsiz(1)+1)-meshres(1)*ii;
        xout(bufsiz(1)+length(xout_tmp)+ii)=xout(bufsiz(1)+length(xout_tmp))+meshres(1)*ii;
    end
    
    yout=zeros(1,2*bufsiz(2)+length(yout_tmp));
    yout(bufsiz(2)+1:bufsiz(2)+length(yout_tmp))=yout_tmp;
    for ii=1:bufsiz(2)
        yout(bufsiz(2)+1-ii)=yout(bufsiz(2)+1)-meshres(2)*ii;
        yout(bufsiz(2)+length(yout_tmp)+ii)=yout(bufsiz(2)+length(yout_tmp))+meshres(2)*ii;
    end
    
    zout=zeros(1,2*bufsiz(3)+length(zout_tmp));
    zout(bufsiz(3)+1:bufsiz(3)+length(zout_tmp))=zout_tmp;
    for ii=1:bufsiz(3)
        zout(bufsiz(3)+1-ii)=zout(bufsiz(3)+1)-meshres(3)*ii;
        zout(bufsiz(3)+length(zout_tmp)+ii)=zout(bufsiz(3)+length(zout_tmp))+meshres(3)*ii;
    end
        
    if(0==1) % Check locations of xin and xout
        figure
        vjunk=zeros(size(xin));vjunk(min(mi):max(mi))=1;
        plot(xin,vjunk,'r-',xin,zeros(size(xin)),'ro',xout,zeros(size(xout)),'b.')

        xout(1)
        xin(min(mi))

        xout(end)
        xin(max(mi))
        xout
    end
    
    % Interpolate Ur, Ui and MagIm with the function SplineInterp3D_withnans. 
    Urint=zeros([length(xout) length(yout) length(zout) 3]);
    Uiint=zeros([length(xout) length(yout) length(zout) 3]);
    %MagIm_int=zeros([length(xout) length(yout) length(zout) 1]);
    if(par)
        n_int=8;
        needsint=zeros([size(MagIm) n_int]);
        needsint(:,:,:,1:3)=Ur;
        needsint(:,:,:,4:6)=Ui;
        needsint(:,:,:,7)=MagIm;
        needsint(:,:,:,8)=regs;
        
        intdata=zeros([length(xout) length(yout) length(zout) n_int]);
               
        parfor ii=1:n_int
            intdata(:,:,:,ii)=SplineInterp3D_withnans(xin,yin,zin,needsint(:,:,:,ii),xout,yout,zout,2,'no');
        end
        Urint(:,:,:,:)=intdata(:,:,:,1:3);
        Uiint(:,:,:,:)=intdata(:,:,:,4:6);
        MagIm_int=intdata(:,:,:,7);
        regs_int=intdata(:,:,:,8);
        clear intdata needsint        
    else
        for ii=1:3
            Urint(:,:,:,ii)=SplineInterp3D_withnans(xin,yin,zin,Ur(:,:,:,ii),xout,yout,zout,2,'no');        
        end
        for ii=1:3
            Uiint(:,:,:,ii)=SplineInterp3D_withnans(xin,yin,zin,Ui(:,:,:,ii),xout,yout,zout,2,'no');
        end
        MagIm_int=SplineInterp3D_withnans(xin,yin,zin,MagIm,xout,yout,zout,2,'no');
        
        if(noreg)
            regs_int=zeros(size(MagIm_int));
        else
            regs_int=SplineInterp3D_withnans(xin,yin,zin,regs,xout,yout,zout,2,'no');
        end
    end

%     if(par)
%         parfor ii=1:3
%             Urint(:,:,:,ii)=SplineInterp3D_withnans(xin,yin,zin,Ur(:,:,:,ii),xout,yout,zout,2,'no');        
%         end
%     else
%         for ii=1:3
%             Urint(:,:,:,ii)=SplineInterp3D_withnans(xin,yin,zin,Ur(:,:,:,ii),xout,yout,zout,2,'no');        
%         end    
%     end
    clear Ur
    disp('Real Displacements Interpolated')
    % Interpolate Ui with the function SplineInterp3D_withnans. 
%    Uiint=zeros([length(xout) length(yout) length(zout) 3]);
%     if(par)
%         parfor ii=1:3
%             Uiint(:,:,:,ii)=SplineInterp3D_withnans(xin,yin,zin,Ui(:,:,:,ii),xout,yout,zout,2,'no');
%         end
%     else
%         for ii=1:3
%             Uiint(:,:,:,ii)=SplineInterp3D_withnans(xin,yin,zin,Ui(:,:,:,ii),xout,yout,zout,2,'no');
%         end
%     end
    clear Ui
    disp('Imag Displacements Interpolated')
    maskint=~isnan(Urint(:,:,:,1));
    clear MagIm
    disp('MagIm Interpolated')
    if(0==1)
        for ii=1:3
            %montagestack(Ur(:,:,:,ii))
            montagestack(Urint(:,:,:,ii))
        end
        montagestack(MagIm_int)
    end    
    
end

if(par)
    matlabpool('close');
end

disp('Displacement Processing Complete, Beginning FE Meshing Process')
tdisp=toc(t1);
tmesh=tic;

%Start the meshing process:
dim=size(maskint);

% Assign node numbers to each interpolated voxel
nodnum=1:numel(Urint(:,:,:,1));
nodnum=reshape(nodnum,size(Urint(:,:,:,1)));

% Start in bottom corner, march elements along mesh, inserting them if they
% are inside the mask
intmp=zeros(prod(floor(dim/2)),27);
nel=0;
nodin=false(size(nodnum));
for ii=1:2:dim(1)-2
    for jj=1:2:dim(2)-2
        for kk=1:2:dim(3)-2
            if(sum(sum(sum(maskint(ii:ii+2,jj:jj+2,kk:kk+2))))==27) % Element is inside mask
                nel=nel+1;
                [intmp(nel,:)]=Hex27incidencelist(nodnum(ii:ii+2,jj:jj+2,kk:kk+2)); % Add element to mesh                
                nodin(ii:ii+2,jj:jj+2,kk:kk+2)=true; % Tag these nodes as included                
            end
        end
    end
end

intmp=intmp(1:nel,:);
tin=toc(tmesh);
tnodes=tic;
% Renumber Nodes (exclude unused nodes)
nod=zeros(prod(dim),3);
uvwr=zeros(prod(dim),3);
uvwi=zeros(prod(dim),3);
magimnod=zeros(prod(dim),1);
regnod=zeros(prod(dim),1);
idx=zeros(prod(dim),3);
old2newnod=zeros(size(nodnum));
%in=zeros(size(intmp));
nn=0;
nbnod=0;
bnod=zeros(prod(dim),1);

for ii=1:dim(1)
    for jj=1:dim(2)
        for kk=1:dim(3)
            if(nodin(ii,jj,kk))
                nn=nn+1;
                old2newnod(ii,jj,kk)=nn;
                %in(intmp==nodnum(ii,jj,kk))=nn;
                nod(nn,:)=[xout(ii) yout(jj) zout(kk)];
                uvwr(nn,:)=[Urint(ii,jj,kk,1) Urint(ii,jj,kk,2) Urint(ii,jj,kk,3)];
                uvwi(nn,:)=[Uiint(ii,jj,kk,1) Uiint(ii,jj,kk,2) Uiint(ii,jj,kk,3)];
                idx(nn,:)=[ii jj kk];
                magimnod(nn)=MagIm_int(ii,jj,kk);
                regnod(nn)=regs_int(ii,jj,kk);
            end
        end
    end    
end
in=old2newnod(intmp); % 330 times faster than in(intmp==nodnum(ii,jj,kk))=nn; inserted in nested loop above. Checked output files are the same with 'diff' command.
clear intmp
tnodrenum=toc(tnodes);

% Resize arrays to correct size.
nod=nod(1:nn,:);
uvwr=uvwr(1:nn,:);
uvwi=uvwi(1:nn,:);
magimnod=magimnod(1:nn);
regnod=regnod(1:nn);
idx=idx(1:nn,:);
bnod=bnod(1:nbnod);

tbni=tic;
disp('Building Element Connetivity')
% Build element connectivity number to determine boundary nodes
elcon=zeros(nn,1);
for ii=1:nel
    for jj=1:27
        elcon(in(ii,jj))=elcon(in(ii,jj))+1;
    end
end
disp('Finding Boundary nodes')
disp(['nn = ' int2str(nn)])
% Build Boundary node array based on nodal connectivity % THis step accounts for ~99% of the meshing time 
bnum=1;


%New bnode finding - 2800x faster than old version below - bnodes come out
%in differet order but get all the same nodes.
bnodchecked=false([nn 1]);
for ii=1:nel
    for jj=1:27
       if(~bnodchecked(in(ii,jj)))
            if jj<=4 || (jj>=10 && jj<=13) % corner nodes =1,2,3,4,10,11,12,13
                if elcon(in(ii,jj))<8 %corner node is boundary node
                    bnod(bnum)=in(ii,jj);
                    bnum=bnum+1;
                end
            elseif (jj>=5 && jj<=8) || (jj>=14 && jj<=17) || (jj>=19 && jj<=22) % mid-edge nodes = 5,6,7,8,14,15,16,17,19,20,21,22
                if elcon(in(ii,jj))<4 %mid-edge node is boundary node
                    bnod(bnum)=in(ii,jj);
                    bnum=bnum+1;
                end
            elseif  jj==9 || jj==18 || (jj>=23 && jj<=26) %mid-fanodce nodes are 9,18,23,24,25,26
                if elcon(in(ii,jj))<2 %mid-face node is boundary node
                    bnod(bnum)=in(ii,jj);
                    bnum=bnum+1;
                end
            end
            bnodchecked(in(ii,jj))=true;
       end               
    end
end
                
    
% Old slow bnod finding
% for ii=1:nn
%     %find out what sort of node they are (corner,mid-edge,midface,center)
%     I=find(in==ii,1); %find the linear index of the first occurance of this node
%     [iel,jnod]=ind2sub(size(in),I);
%     if jnod<=4 || (jnod>=10 && jnod<=13) % corner nodes =1,2,3,4,10,11,12,13
%         if elcon(ii)<8 %corner node is boundary node
%             bnod(bnum)=ii;
%             bnum=bnum+1;
%         end
%     elseif (jnod>=5 && jnod<=8) || (jnod>=14 && jnod<=17) || (jnod>=19 && jnod<=22) % mid-edge nodes = 5,6,7,8,14,15,16,17,19,20,21,22
%         if elcon(ii)<4 %mid-edge node is boundary node
%             bnod(bnum)=ii;
%             bnum=bnum+1;
%         end
%     elseif  jnod==9 || jnod==18 || (jnod>=23 && jnod<=26) %mid-fanodce nodes are 9,18,23,24,25,26
%         if elcon(ii)<2 %mid-face node is boundary node
%             bnod(bnum)=ii;
%             bnum=bnum+1;
%         end
%     elseif jnod~=27 %something has gone wrong
%         disp(['ERROR: type of node not recognised, jnod =' int2str(jnod)])
%     end
% end


nbnod=bnum-1;
tbnod=toc(tbni);


% Perform coordinate transformation
uvwr=uvwr*MPSto123';
uvwi=uvwi*MPSto123';  % Nnx3*3x3=Nn*3, Transpose IS required

tfemesh=toc(tmesh); % Time for full meshing process

%% Calcualte Zone sizes
touti=tic;
% Zonestrat=1: target number of nodes per zone
% Zonestrat=2: target number of wavelengths per zone
znovlp=[0.15 0.15 0.15]; % Zone overlap
znovlp1=1+znovlp;
%Calculate edge lenth factors
rangedim=max(nod)-min(nod);
disp('Mesh Size Ratio (x,y,z):')
disp(rangedim./min(rangedim))

if(zonestrat==1)
    
    % -> Aim for approximatly cubic zones, with a specified number of nodes per zone
    % npz = Target number of nodes per zone (Note that this strategy usually ends up with zones ~70% of target size).
    disp(['zonestrat=1, aiming for ' int2str(npz) ' nodes per zone'])
    zfac=rangedim./min(rangedim); % Ratio of mesh dimensions
   
    % Nn/Nz=Npz/(znovlp+1)^3
    % Nz=Nn*(znovrlp+1)^3/Npz
    % zx*zy*zz = Nz
    % zfacx*F * zfacy*F * zfacz*F = Nz
    % F = (Nz/(zfacx*zfacy*zfacz))^(1/3)
    % zx =zfacx*F, zy=zfacy*F, zz=zfacz*Ffprintf(fid,'%c',direc);fprintf(fid,'/');fprintf(fid,'%c',elmf);fprintf(fid,'/n')
    
    Nz=nn*znovlp1(1)*znovlp1(2)*znovlp1(3)/npz;
    F=(Nz/(zfac(1)*zfac(2)*zfac(3)))^(1/3);
    znedge=zfac.*F;

    znedgeint=round(znedge);
    for ii=1:3
        if znedgeint(ii)==0
            znedgeint(ii)=1;
        end
    end
    disp(['v <= 7.04 Zone edge factors calculated = ' int2str(znedgeint)])
    disp(['Estimated nodes per zone = ' num2str(nn*prod(znovlp1)/(prod(znedgeint)))])
    
    % v7.05 = direct specification of zone sizes, [L1 L2 L3]
    % Aim for cubic zones, NPZ =(L*(1+2*ovlp(1)))/meshres(1))*(L*(1+2*ovlp(2)))/meshres(2))*(L*(1+2*ovlp(3)))/meshres(3))
    % L^3=NPZ*mesres(1)*meshres(2)*meshres(3)/(1+2*ovlp(1))*(1+2*ovlp(2))*1+2*ovlp(3)))
    znedgelength(1:3)=( npz*meshres(1)*meshres(2)*meshres(3) / ((1+2*znovlp(1))*(1+2*znovlp(2))*(1+2*znovlp(3))) )^(1/3);
    disp(['v7.05 zone edge length ' num2str(znedgelength)])
    disp(['Estimated nodes per zone = ' int2str(prod(znedgelength.*(1+2*znovlp)./meshres))])
elseif(zonestrat==2)
    disp(['zonestrat=2, aiming for ' int2str(wlperzone) ' wavelengths per zone'])    
     
    znesz=wlperzone*wl;
    znedgeint=round(rangedim./znesz.*(1+2*znovlp));
    for ii=1:3       
        if(znedgeint(ii)<3) % Mostly Edge zones, only one SZ overlap
            znedgeint(ii)=round(rangedim(ii)./znesz.*(1+znovlp(ii)));
        end
        if znedgeint(ii)==0
            znedgeint(ii)=1;
        end
    end
    
    disp(['v <= 7.04 Zone edge factors calculated = ' int2str(znedgeint)])  
    disp(['Estimated nodes per zone = ' num2str(nn*prod(znovlp1)/(prod(znedgeint)))])
    
    % v7.05 = direct specification of zone sizes, [L1 L2 L3]
    % Actual SZ length = (1+2*znovlp(i))*L(i)
    for ii=1:3
        znedgelength(ii)=znesz/(1+2*znovlp(ii));
    end
    disp(['v7.05 zone edge length ' num2str(znedgelength)])
    disp(['Estimated nodes per zone = ' num2str(prod(znesz./meshres)) ])
elseif(zonestrat==3)
    
    znedgeint=round(rangedim./zlength.*(1+2*znovlp));
    for ii=1:3       
        if(znedgeint(ii)<3) % Mostly Edge zones, only one SZ overlap
            znedgeint(ii)=round(rangedim(ii)./zlength.*(1+znovlp(ii)));
        end
        if znedgeint(ii)==0
            znedgeint(ii)=1;
        end
    end
    disp(['v <= 7.04 Zone edge factors calculated = ' int2str(znedgeint)])  

    for ii=1:3
        znedgelength(ii)=zlength;
    end
    disp(['v7.05 zone edge length ' num2str(znedgelength)])
end

disp(['Meshing Complete, ' int2str(nn) ' Nodes and ' int2str(nel) ' elements'])


%% Output mesh

% Output everything to files

%mask file
maskoutf=[outstm '.mask.mat'];
save(maskoutf,'mask');

%nod files
ind=(1:nn)';
nodoutf=[outstm '.nod'];
nodhmgoutf=[outstm '.hom.nod'];
nodregoutf=[outstm '.reg.nod'];


fid=fopen(nodoutf,'w');
fprintf(fid,'%7i %15.8e %15.8e %15.8e %12.5e \n',[ind nod ones(nn,1)]');
fclose(fid);

fid=fopen(nodhmgoutf,'w');
fprintf(fid,'%7i %15.8e %15.8e %15.8e  1 \n',[ind nod]');
fclose(fid);

fid=fopen(nodregoutf,'w');
fprintf(fid,'%7i %15.8e %15.8e %15.8e %3i \n',[ind nod regnod]');
fclose(fid);


%elm files
elmind=(1:nel)';
elmf=[outstm '.elm'];
fid=fopen(elmf,'w');
fprintf(fid,'%7i %8i %8i %8i %8i %8i %8i %8i %8i %8i %8i %8i %8i %8i %8i %8i %8i %8i %8i %8i %8i %8i %8i %8i %8i %8i %8i %8i \n',[elmind in]');
fclose(fid);

% Element file for MRE-Zone-v6 (extra column of ones to indicate
% homogeneous elementally defined properties
elmhmgf=[outstm '.hom.elm'];
fid=fopen(elmhmgf,'w');
fprintf(fid,'%7i %8i %8i %8i %8i %8i %8i %8i %8i %8i %8i %8i %8i %8i %8i %8i %8i %8i %8i %8i %8i %8i %8i %8i %8i %8i %8i %8i  1 \n',[elmind in]');
fclose(fid);


%index file -> Note this is interpolated index. Need to back-interpolate
%using xin and xout to get values at original MR voxels.
idxf=[outstm '.idx'];
fid=fopen(idxf,'w');
fprintf(fid,'%7i %6i %6i %6i\n',[ind idx]');
fclose(fid);

%disp file
dspoutf=[outstm '.dsp'];
fid=fopen(dspoutf,'w');
fprintf(fid,'%7i %15.8e %15.8e %15.8e %15.8e %15.8e %15.8e \n',[ind uvwr(:,1) uvwi(:,1) uvwr(:,2) uvwi(:,2) uvwr(:,3) uvwi(:,3)]');
fclose(fid);

%magim.mtrin file
magoutf=[outstm 'magim.mtrin'];
fid=fopen(magoutf,'w');
fprintf(fid,'%7i %15.8e\n',[ind magimnod]');
fclose(fid);


%bnod file
bcoutf=[outstm '.bnd'];
fid=fopen(bcoutf,'w');
fprintf(fid,'%7i  %7i\n',[(1:nbnod)' bnod']');
fclose(fid);


% region file for soft prior -> Output the region file at the original data
% resolution so that the regions can be interpolated to the material
% property meshes inside the recon code.
regoutf=[outstm '.regstack'];
% Format:
% StackSize
% x coords
% y coords
% z coords
% slice1
% slice2
% ..
% Slice_end
fid=fopen(regoutf,'w');
fprintf(fid,'Array Dimensions \n'); 
fprintf(fid,'%7i %7i %7i \n',mridim);
fprintf(fid,'x coordinates \n');
fclose(fid);
dlmwrite(regoutf,xin,'-append','delimiter',' '); 
fid=fopen(regoutf,'a'); fprintf(fid,'y coordniates \n'); fclose(fid);
dlmwrite(regoutf,yin,'-append','delimiter',' ');
fid=fopen(regoutf,'a'); fprintf(fid,'z coordniates \n'); fclose(fid);
dlmwrite(regoutf,zin,'-append','delimiter',' ');
for ii=1:mridim(3)
    dlmwrite(regoutf,regs(:,:,ii),'-append','delimiter',' ');
end
  

  
  
  



%save data required to get back to original MR voxels (including mask filename)
if isempty(msk)
    msk=maskoutf;
end
save([outstm '.InterpLocations.mat'],'maskint','nodin','xout','yout','zout','xin','yin','zin','msk') 

%save interpolated motion data for possible later use
save([outstm '.InterpData.mat'],'maskint','Urint','Uiint','MagIm_int','MPSto123')

%save inputs to meshing code to link with mesh
save([outstm '.meshinput.mat'],'meshstrat','zonestrat','dispscale','dispscalar')

if((meshstrat==2)||(zonestrat==2))
    save([outstm '.meshinput.mat'],'muest','rhoest','-append')
end
if(meshstrat==2)||(meshstrat==3)
    save([outstm '.meshinput.mat'],'targetres','meshres','-append')
end
if(meshstrat==2)
    save([outstm '.meshinput.mat'],'npw','-append')
end

if(zonestrat==2)
    save([outstm '.meshinput.mat'],'wlperzone','-append')
end

if(zonestrat==2)
    save([outstm '.meshinput.mat'],'wlperzone','-append')
end


%% Output Runfiles and submit files

% Make sure mu and rho estiamtes are there
if(meshstrat~=2)&&(zonestrat~=2)
    muest=mudef;
    rhoest=1000;
end
DR=0.18; % Default damping ratio
%DR=0.10; % Default damping ratio
%DR=0.025;
vincomp=0.499; % Default Incompressible Damping Ratio
vcomp=0.4; % Default compressible damping ratio


% Create Directories if they dont exist
D=dir; % Check to see if hex exists
hexflg=0;
for ii=1:length(D)
    if length(D(ii).name)==3
        if (strcmp(D(ii).name,'hex')&&(D(ii).isdir==1)) 
            hexflg=1;
        end
    end
end
if (hexflg==0) % Make Hex
    disp('Creating hex directory')
    mkdir('hex')
end
junk=['hex/' outstm]; % Make subdirectiories
mkdir('hex',outstm);
mkdir(junk,'inv');


direc=pwd;
inpath=[direc '/hex/' outstm '/'];
outpath=['inv/'];

% Initialize list of files to move
mvfiles(1).name=[outstm '*']; % Mesh files


%% MREv7.3 runfile
% Iso incompressible run file - viscoelastic - Soft Prior On

if(~noreg) % do not output soft prior runfiles without a supplied segmentation.

viscoutputfv7p3isoincomp=[outstm  '_G' sprintf('%4.0f',muest) '.v7.3.inv.iso.incomp.visc_SPon'];
viscrfilev7p3='runfile-v7p3visc.dat_SPon';
mvfiles(length(mvfiles)+1).name=viscrfilev7p3;
fid=fopen(viscrfilev7p3,'w');
fprintf(fid,'Subzone Reconstruction Data File \n');
fprintf(fid,'Problem Type <0=Forward Problem Only> <1+=Inverse Problem> \n');
fprintf(fid,'1 \n');
fprintf(fid,'Node File: \n');
fprintf(fid,'%c',nodhmgoutf);fprintf(fid,'\n');
fprintf(fid,'Element File: \n');
fprintf(fid,'%c',elmhmgf);fprintf(fid,'\n');
fprintf(fid,'Boundary Node File: \n');
fprintf(fid,'%c',bcoutf);fprintf(fid,'\n');
fprintf(fid,'Region Stack File: \n');
fprintf(fid,'%c',regoutf);fprintf(fid,'\n');
fprintf(fid,'Measured Displacement File: \n');
fprintf(fid,'1 \n');
fprintf(fid,'%10.4fd0',freqHz);fprintf(fid,'\n');
fprintf(fid,'%c',dspoutf);fprintf(fid,'\n');
fprintf(fid,'Initial Solution File: \n');
fprintf(fid,'0 \n');
fprintf(fid,'Output File Stem: \n');
fprintf(fid,'%c',outpath);fprintf(fid,'%c',viscoutputfv7p3isoincomp);fprintf(fid,'\n');
fprintf(fid,'Print Detailed Runtime Debugging and Execution Comments <verb> <file>:\n');
fprintf(fid,'.false.,.false. \n');
fprintf(fid,'Material Model <1=isotropic> <2=orthotropic> <3=iso compress>: \n');
fprintf(fid,'1 \n');
fprintf(fid,'Number of Material Properties: \n');
fprintf(fid,'3 \n');
fprintf(fid,'Material Description Style <1=nodal> <2=element> [<shear modulus> <density> <bulk modulus>]: \n');
fprintf(fid,'1,1,1 \n');
fprintf(fid,'Number of Parameters per Material Point: \n');
fprintf(fid,'1,1,1 \n');
fprintf(fid,'Property Scalars: \n');
fprintf(fid,'%5.0f.d0,%5.0f.d0,%5.0f.d0,%8.5fd0,%10.0f.d0,%7.0f.d0',[muest 2*DR*muest rhoest -1d-1 2*muest*(1+vincomp)/(3*(1-2*vincomp)) 0]);fprintf(fid,'\n');
fprintf(fid,'Number of Material Property Mesh Resolutions: \n');
fprintf(fid,'3 \n');
fprintf(fid,'Material Property Mesh Resolutions (x,y,z): \n');
fprintf(fid,'%9.7f , %9.7f , %9.7f \n',DirIndex(4,1:3)*1e-3);
fprintf(fid,'%9.7f , %9.7f , %9.7f \n',5.*DirIndex(4,1:3)*1e-3);
fprintf(fid,'%9.7f , %9.7f , %9.7f \n',10.*DirIndex(4,1:3)*1e-3);
fprintf(fid,'Material Mesh Index (1st line real part, 2nd line imag) \n ');
fprintf(fid,'1 2 2 \n');
fprintf(fid,'1 2 2 \n');
fprintf(fid,'Reconstruction Indicators: \n');
fprintf(fid,'.true.,.true. \n');
fprintf(fid,'.false.,.false. \n');
fprintf(fid,'.false.,.false. \n');
fprintf(fid,'Property Estimate Variance Calculation: \n');
fprintf(fid,'.false. \n');
fprintf(fid,'Zone Sizes (not including overlap factor, actual size = (1+2*ovlp)*siz) [x y z]: \n');
fprintf(fid,'%12.5e,%12.5e,%12.5e',znedgelength);fprintf(fid,'\n');
fprintf(fid,'Zone Overlap Percentage [x y z]: \n');
fprintf(fid,'%4.3fd0,%4.3fd0,%4.3fd0',znovlp);fprintf(fid,'\n');
fprintf(fid,'Iteration Limits [global zone]: \n');
fprintf(fid,'100,1 \n');
fprintf(fid,'Minimum Parameter Update Size [global zone line]: \n');
fprintf(fid,'%c','0.1d0, 0.1d0, 0.1d0');fprintf(fid,'\n');
fprintf(fid,'Number of Zone Iteration Structures: \n');
fprintf(fid,'3 \n');
fprintf(fid,'Iteration Limits for Zone Iteration Structures (NOTE:  provide one less limit than # of structures!!!): \n');
fprintf(fid,'10,150 \n');
fprintf(fid,'Zone Iteration Structures [<# of CG iters> <# of GN iters> <!!! QN CURRENTLY UNAVAILABLE !!!> <# of line search iters>]: \n');
fprintf(fid,'1,0,0,1 \n');
fprintf(fid,'2,0,0,2 \n');
fprintf(fid,'3,0,0,2 \n');
fprintf(fid,'Number of Processors per MUMPS Communicator \n');
fprintf(fid,'1 \n');
fprintf(fid,'Maximum Amount of RAM per Processor [MB] \n');
fprintf(fid,'999999999 \n');
fprintf(fid,'ooooooooooo  REGULARIZATION DESCRIPTORS  ooooooooooooooooooooooooooooooooooooooooooooooooooooooooooo \n');
fprintf(fid,'Regularization Indicators [<TV> <SF> <TK> <MQ> <JW> <constraints> <CG residual> <VH> <McG> <soft prior>]: \n');
fprintf(fid,'.false.,.true.,.false.,.false.,.false.,.false.,.true.,.true.,.true.,.true. \n');
fprintf(fid,'Number of constant regularization iterations (final N iterations use final regularization weights) \n');
fprintf(fid,'30 \n');
fprintf(fid,'Number of Parameters to Treat with Total Variation: \n');
fprintf(fid,'3 \n');
fprintf(fid,'TV Parameter List: \n');
fprintf(fid,'1,2,3 \n');
fprintf(fid,'TV Delta Values: <1st line real, 2nd line imag> \n');
fprintf(fid,'1.d-19,1.d-19,1.d-19 \n');
fprintf(fid,'1.d-19,1.d-19,1.d-19 \n');
fprintf(fid,'TV Initial Weights: <1st line real, 2nd line imag> \n');
fprintf(fid,'5.d-16,5.d-16,5.d-16 \n');
fprintf(fid,'5.d-16,5.d-16,5.d-16 \n');
fprintf(fid,'TV Final Weights: <1st line real, 2nd line imag> \n');
fprintf(fid,'5.d-16,5.d-16,5.d-16 \n');
fprintf(fid,'5.d-16,5.d-16,5.d-16 \n');
fprintf(fid,'Number of Parameters to Treat with Spatial Filtering: \n');
fprintf(fid,'3 \n');
fprintf(fid,'SF Parameter List: \n');
fprintf(fid,'1,2,3 \n');
fprintf(fid,'SF Gaussian filter initial Widths: <1st line real, 2nd line imag> \n');
fprintf(fid,'0.003d0,0.003d0,0.003d0 \n');
fprintf(fid,'0.003d0,0.003d0,0.003d0 \n');
fprintf(fid,'SF Final Gaussian width: <1st line real, 2nd line imag> \n');
fprintf(fid,'0.0015d0,0.0015d0,0.0015d0 \n');
fprintf(fid,'0.0015d0,0.0015d0,0.0015d0 \n');
fprintf(fid,'Number of Parameters to Treat with Tikhonov Regularization: \n');
fprintf(fid,'3 \n');
fprintf(fid,'TK Parameter List: \n');
fprintf(fid,'1,2,3 \n');
fprintf(fid,'TK Initial Weights: <1st line real, 2nd line imag> \n');
fprintf(fid,'1.d-18,1.d-18,1.d-18 \n');
fprintf(fid,'1.d-18,1.d-18,1.d-18 \n');
fprintf(fid,'TK Final Weights: <1st line real, 2nd line imag> \n');
fprintf(fid,'1.d-18,1.d-18,1.d-18 \n');
fprintf(fid,'1.d-18,1.d-18,1.d-18 \n');
fprintf(fid,'Number of Parameters to Treat with Marquardt Regularization: \n');
fprintf(fid,'3 \n');
fprintf(fid,'Distance of alpha from 1 at which MQ weights are adjusted: \n');
fprintf(fid,'0.25d0 \n');
fprintf(fid,'MQ Parameter List: \n');
fprintf(fid,'1,2,3 \n');
fprintf(fid,'MQ Weight Delta: <1st line real, 2nd line imag> \n');
fprintf(fid,'0.5d0,0.5d0,0.5d0 \n');
fprintf(fid,'0.5d0,0.5d0,0.5d0 \n');
fprintf(fid,'MQ Initial Weights: <1st line real, 2nd line imag> \n');
fprintf(fid,'1.d2,1.d2,1.d2 \n');
fprintf(fid,'1.d2,1.d2,1.d2 \n');
fprintf(fid,'MQ Minimum Weights: <1st line real, 2nd line imag> \n');
fprintf(fid,'1.d-11,1.d-11,1.d-11 \n');
fprintf(fid,'1.d-11,1.d-11,1.d-11 \n');
fprintf(fid,'Number of Parameters to Treat with Joachimowicz Regularization: \n');
fprintf(fid,'3 \n');
fprintf(fid,'JW Parameter List: \n');
fprintf(fid,'1,2,3 \n');
fprintf(fid,'JW Initial Weights: <1st line real, 2nd line imag> \n');
fprintf(fid,'5.d0,5.d0,5.d0 \n');
fprintf(fid,'5.d0,5.d0,5.d0 \n');
fprintf(fid,'JW Final Weights: <1st line real, 2nd line imag> \n');
fprintf(fid,'5.d0,5.d0,5.d0 \n');
fprintf(fid,'5.d0,5.d0,5.d0 \n');
fprintf(fid,'Number of Parameters to Treat with Constraints: \n');
fprintf(fid,'3 \n');
fprintf(fid,'Constraint Parameter List: \n');
fprintf(fid,'1,2,3 \n');
fprintf(fid,'Constraint Weights <converges to exact solution as weight --> inf>: \n');
fprintf(fid,'1.d-14,1.d-14,1.d-14 \n');
fprintf(fid,'Constraint Minimums: \n');
fprintf(fid,'200.d0, 1.d0, 1000.d0, -5.d4, 1000.d0, 0.d0 \n');
fprintf(fid,'Constraint Maximums: \n');
fprintf(fid,'5.d5, 5.d5, 1000.d0, -1.d-2, 1.d12, 0.d0 \n');
fprintf(fid,'CG Residual Scaling Initial Weights: <1st line real, 2nd line imag> \n');
fprintf(fid,'0.04d0,0.04d0,0.04d0 \n');
fprintf(fid,'0.04d0,0.04d0,0.04d0 \n');
fprintf(fid,'CG Residual Scaling Final Weights: <1st line real, 2nd line imag> \n');
fprintf(fid,'0.01d0,0.01d0,0.01d0 \n');
fprintf(fid,'0.01d0,0.01d0,0.01d0 \n');
fprintf(fid,'Van Houten Regularization Level: \n');
fprintf(fid,'1.2d0 \n');
fprintf(fid,'Number of parameters to treat with soft prior: \n');
fprintf(fid,'3 \n');
fprintf(fid,'SP Parameter List: \n');
fprintf(fid,'1,2,3 \n');
fprintf(fid,'SP Initial Weights: <1st line real, 2nd line imag> \n');
fprintf(fid,'1.d-9,1.d-9,1.d-9 \n');
fprintf(fid,'1.d-9,1.d-9,1.d-9 \n');
fprintf(fid,'SP Final Weights: <1st line real, 2nd line imag> \n');
fprintf(fid,'1.d-9,1.d-9,1.d-9 \n');
fprintf(fid,'1.d-9,1.d-9,1.d-9 \n');
fprintf(fid,'SP start iteration delay (first N iterations do not use soft prior \n');
fprintf(fid,'0 \n');
fclose(fid);

% Build MRE-submit file
viscsubnamev7p3='MPICH2-Submitv7.3visc_SPon';
fid=fopen(viscsubnamev7p3,'w');
mvfiles(length(mvfiles)+1).name=viscsubnamev7p3;
fprintf(fid,'#!/bin/bash -l \n');
fprintf(fid,'# declare a name for this job \n'); 
fprintf(fid,'#PBS -N ');fprintf(fid,'%c',viscoutputfv7p3isoincomp);fprintf(fid,'\n');
fprintf(fid,'# request the queue (enter the possible names, if omitted, serial is the default) \n');
fprintf(fid,'#PBS -q default  \n');
fprintf(fid,'# request  node  \n');
fprintf(fid,'#PBS -l nodes=4:ppn=8 \n');
fprintf(fid,'#PBS -l feature=amd \n');
fprintf(fid,'# request some hours of wall time  \n');
fprintf(fid,'#PBS -l walltime=36:00:00  \n');
fprintf(fid,'#combine PBS standard output and error files  \n');
fprintf(fid,'#PBS -j oe  \n');
fprintf(fid,'# mail is sent to you when the job starts and when it terminates or aborts  \n');
fprintf(fid,'##PBS -m bea  \n');
fprintf(fid,'# specify your email address  \n');
fprintf(fid,'##PBS -M matthew.d.mcgarry@dartmouth.edu  \n');
fprintf(fid,'# Change to Submission Directory  \n');
fprintf(fid,'cd $PBS_O_WORKDIR  \n');
fprintf(fid,'# run the program \n');
fprintf(fid,'cat $PBS_NODEFILE | uniq > node_file \n');
fprintf(fid,'nnodes=$(cat node_file | wc -l) \n');
fprintf(fid,'nprocs=$(cat $PBS_NODEFILE | wc -l) \n');
fprintf(fid,'export MKL_NUM_THREADS=1 \n');
fprintf(fid,'echo Nodes $nnodes \n');
fprintf(fid,'echo Procs $nprocs \n');
fprintf(fid,'mpirun -np $nprocs -hostfile $PBS_NODEFILE /ihome/mmcgarry/code/MREv7p35/MRE-Zone.v7.35.discov ');fprintf(fid,'%c',[viscrfilev7p3]);fprintf(fid,'\n');
fprintf(fid,'exit 0  \n');
fclose(fid);

end

%% MREv7.3 runfile
% Iso incompressible run file - viscoelastic - No soft Prior
viscoutputfv7p3isoincomp_noSP=[outstm  '_G' sprintf('%4.0f',muest) '.v7.3.inv.iso.incomp.visc_SPoff_SF0p0015_CG2p2'];
viscrfilev7p3_noSP='runfile-v7p3visc.dat';
mvfiles(length(mvfiles)+1).name=viscrfilev7p3_noSP;
fid=fopen(viscrfilev7p3_noSP,'w');
fprintf(fid,'Subzone Reconstruction Data File \n');
fprintf(fid,'Problem Type <0=Forward Problem Only> <1+=Inverse Problem> \n');
fprintf(fid,'1 \n');
fprintf(fid,'Node File: \n');
fprintf(fid,'%c',nodhmgoutf);fprintf(fid,'\n');
fprintf(fid,'Element File: \n');
fprintf(fid,'%c',elmhmgf);fprintf(fid,'\n');
fprintf(fid,'Boundary Node File: \n');
fprintf(fid,'%c',bcoutf);fprintf(fid,'\n');
fprintf(fid,'Region Stack File: \n');
fprintf(fid,'%c',regoutf);fprintf(fid,'\n');
fprintf(fid,'Measured Displacement File: \n');
fprintf(fid,'1 \n');
fprintf(fid,'%10.4fd0',freqHz);fprintf(fid,'\n');
fprintf(fid,'%c',dspoutf);fprintf(fid,'\n');
fprintf(fid,'Initial Solution File: \n');
fprintf(fid,'0 \n');
fprintf(fid,'Output File Stem: \n');
fprintf(fid,'%c',outpath);fprintf(fid,'%c',viscoutputfv7p3isoincomp_noSP);fprintf(fid,'\n');
fprintf(fid,'Print Detailed Runtime Debugging and Execution Comments <verb> <file>:\n');
fprintf(fid,'.false.,.false. \n');
fprintf(fid,'Material Model <1=isotropic> <2=orthotropic> <3=iso compress>: \n');
fprintf(fid,'1 \n');
fprintf(fid,'Number of Material Properties: \n');
fprintf(fid,'3 \n');
fprintf(fid,'Material Description Style <1=nodal> <2=element> [<shear modulus> <density> <bulk modulus>]: \n');
fprintf(fid,'1,1,1 \n');
fprintf(fid,'Number of Parameters per Material Point: \n');
fprintf(fid,'1,1,1 \n');
fprintf(fid,'Property Scalars: \n');
fprintf(fid,'%5.0f.d0,%5.0f.d0,%5.0f.d0,%8.5fd0,%10.0f.d0,%7.0f.d0',[muest 2*DR*muest rhoest -1d-1 2*muest*(1+vincomp)/(3*(1-2*vincomp)) 0]);fprintf(fid,'\n');
fprintf(fid,'Number of Material Property Mesh Resolutions: \n');
fprintf(fid,'3 \n');
fprintf(fid,'Material Property Mesh Resolutions (x,y,z): \n');
fprintf(fid,'%9.7f , %9.7f , %9.7f \n',DirIndex(4,1:3)*1e-3);
fprintf(fid,'%9.7f , %9.7f , %9.7f \n',5.*DirIndex(4,1:3)*1e-3);
fprintf(fid,'%9.7f , %9.7f , %9.7f \n',10.*DirIndex(4,1:3)*1e-3);
fprintf(fid,'Material Mesh Index (1st line real part, 2nd line imag) \n ');
fprintf(fid,'1 2 2 \n');
fprintf(fid,'1 2 2 \n');
fprintf(fid,'Reconstruction Indicators: \n');
fprintf(fid,'.true.,.true. \n');
fprintf(fid,'.false.,.false. \n');
fprintf(fid,'.false.,.false. \n');
fprintf(fid,'Property Estimate Variance Calculation: \n');
fprintf(fid,'.false. \n');
fprintf(fid,'Zone Sizes (not including overlap factor, actual size = (1+2*ovlp)*siz) [x y z]: \n');
fprintf(fid,'%12.5e,%12.5e,%12.5e',znedgelength);fprintf(fid,'\n');
fprintf(fid,'Zone Overlap Percentage [x y z]: \n');
fprintf(fid,'%4.3fd0,%4.3fd0,%4.3fd0',znovlp);fprintf(fid,'\n');
fprintf(fid,'Iteration Limits [global zone]: \n');
fprintf(fid,'100,1 \n');
fprintf(fid,'Minimum Parameter Update Size [global zone line]: \n');
fprintf(fid,'%c','0.1d0, 0.1d0, 0.1d0');fprintf(fid,'\n');
fprintf(fid,'Number of Zone Iteration Structures: \n');
fprintf(fid,'3 \n');
fprintf(fid,'Iteration Limits for Zone Iteration Structures (NOTE:  provide one less limit than # of structures!!!): \n');
fprintf(fid,'10,150 \n');
fprintf(fid,'Zone Iteration Structures [<# of CG iters> <# of GN iters> <!!! QN CURRENTLY UNAVAILABLE !!!> <# of line search iters>]: \n');
fprintf(fid,'1,0,0,1 \n');
fprintf(fid,'2,0,0,2 \n');
fprintf(fid,'3,0,0,2 \n');
fprintf(fid,'Number of Processors per MUMPS Communicator \n');
fprintf(fid,'1 \n');
fprintf(fid,'Maximum Amount of RAM per Processor [MB] \n');
fprintf(fid,'999999999 \n');
fprintf(fid,'ooooooooooo  REGULARIZATION DESCRIPTORS  ooooooooooooooooooooooooooooooooooooooooooooooooooooooooooo \n');
fprintf(fid,'Regularization Indicators [<TV> <SF> <TK> <MQ> <JW> <constraints> <CG residual> <VH> <McG> <soft prior>]: \n');
fprintf(fid,'.false.,.true.,.false.,.false.,.false.,.false.,.true.,.true.,.true.,.false. \n');
fprintf(fid,'Number of constant regularization iterations (final N iterations use final regularization weights) \n');
fprintf(fid,'30 \n');
fprintf(fid,'Number of Parameters to Treat with Total Variation: \n');
fprintf(fid,'3 \n');
fprintf(fid,'TV Parameter List: \n');
fprintf(fid,'1,2,3 \n');
fprintf(fid,'TV Delta Values: <1st line real, 2nd line imag> \n');
fprintf(fid,'1.d-19,1.d-19,1.d-19 \n');
fprintf(fid,'1.d-19,1.d-19,1.d-19 \n');
fprintf(fid,'TV Initial Weights: <1st line real, 2nd line imag> \n');
fprintf(fid,'5.d-16,5.d-16,5.d-16 \n');
fprintf(fid,'5.d-16,5.d-16,5.d-16 \n');
fprintf(fid,'TV Final Weights: <1st line real, 2nd line imag> \n');
fprintf(fid,'5.d-16,5.d-16,5.d-16 \n');
fprintf(fid,'5.d-16,5.d-16,5.d-16 \n');
fprintf(fid,'Number of Parameters to Treat with Spatial Filtering: \n');
fprintf(fid,'3 \n');
fprintf(fid,'SF Parameter List: \n');
fprintf(fid,'1,2,3 \n');
fprintf(fid,'SF Gaussian filter initial Widths: <1st line real, 2nd line imag> \n');
fprintf(fid,'0.003d0,0.003d0,0.003d0 \n');
fprintf(fid,'0.003d0,0.003d0,0.003d0 \n');
fprintf(fid,'SF Final Gaussian width: <1st line real, 2nd line imag> \n');
fprintf(fid,'0.0015d0,0.0015d0,0.0015d0 \n');
fprintf(fid,'0.0015d0,0.0015d0,0.0015d0 \n');
fprintf(fid,'Number of Parameters to Treat with Tikhonov Regularization: \n');
fprintf(fid,'3 \n');
fprintf(fid,'TK Parameter List: \n');
fprintf(fid,'1,2,3 \n');
fprintf(fid,'TK Initial Weights: <1st line real, 2nd line imag> \n');
fprintf(fid,'1.d-18,1.d-18,1.d-18 \n');
fprintf(fid,'1.d-18,1.d-18,1.d-18 \n');
fprintf(fid,'TK Final Weights: <1st line real, 2nd line imag> \n');
fprintf(fid,'1.d-18,1.d-18,1.d-18 \n');
fprintf(fid,'1.d-18,1.d-18,1.d-18 \n');
fprintf(fid,'Number of Parameters to Treat with Marquardt Regularization: \n');
fprintf(fid,'3 \n');
fprintf(fid,'Distance of alpha from 1 at which MQ weights are adjusted: \n');
fprintf(fid,'0.25d0 \n');
fprintf(fid,'MQ Parameter List: \n');
fprintf(fid,'1,2,3 \n');
fprintf(fid,'MQ Weight Delta: <1st line real, 2nd line imag> \n');
fprintf(fid,'0.5d0,0.5d0,0.5d0 \n');
fprintf(fid,'0.5d0,0.5d0,0.5d0 \n');
fprintf(fid,'MQ Initial Weights: <1st line real, 2nd line imag> \n');
fprintf(fid,'1.d2,1.d2,1.d2 \n');
fprintf(fid,'1.d2,1.d2,1.d2 \n');
fprintf(fid,'MQ Minimum Weights: <1st line real, 2nd line imag> \n');
fprintf(fid,'1.d-11,1.d-11,1.d-11 \n');
fprintf(fid,'1.d-11,1.d-11,1.d-11 \n');
fprintf(fid,'Number of Parameters to Treat with Joachimowicz Regularization: \n');
fprintf(fid,'3 \n');
fprintf(fid,'JW Parameter List: \n');
fprintf(fid,'1,2,3 \n');
fprintf(fid,'JW Initial Weights: <1st line real, 2nd line imag> \n');
fprintf(fid,'5.d0,5.d0,5.d0 \n');
fprintf(fid,'5.d0,5.d0,5.d0 \n');
fprintf(fid,'JW Final Weights: <1st line real, 2nd line imag> \n');
fprintf(fid,'5.d0,5.d0,5.d0 \n');
fprintf(fid,'5.d0,5.d0,5.d0 \n');
fprintf(fid,'Number of Parameters to Treat with Constraints: \n');
fprintf(fid,'3 \n');
fprintf(fid,'Constraint Parameter List: \n');
fprintf(fid,'1,2,3 \n');
fprintf(fid,'Constraint Weights <converges to exact solution as weight --> inf>: \n');
fprintf(fid,'1.d-14,1.d-14,1.d-14 \n');
fprintf(fid,'Constraint Minimums: \n');
fprintf(fid,'200.d0, 1.d0, 1000.d0, -5.d4, 1000.d0, 0.d0 \n');
fprintf(fid,'Constraint Maximums: \n');
fprintf(fid,'5.d5, 5.d5, 1000.d0, -1.d-2, 1.d12, 0.d0 \n');
fprintf(fid,'CG Residual Scaling Initial Weights: <1st line real, 2nd line imag> \n');
fprintf(fid,'0.04d0,0.04d0,0.04d0 \n');
fprintf(fid,'0.04d0,0.04d0,0.04d0 \n');
fprintf(fid,'CG Residual Scaling Final Weights: <1st line real, 2nd line imag> \n');
fprintf(fid,'0.01d0,0.01d0,0.01d0 \n');
fprintf(fid,'0.01d0,0.01d0,0.01d0 \n');
fprintf(fid,'Van Houten Regularization Level: \n');
fprintf(fid,'1.2d0 \n');
fprintf(fid,'Number of parameters to treat with soft prior: \n');
fprintf(fid,'3 \n');
fprintf(fid,'SP Parameter List: \n');
fprintf(fid,'1,2,3 \n');
fprintf(fid,'SP Initial Weights: <1st line real, 2nd line imag> \n');
fprintf(fid,'1.d-11,1.d-11,1.d-11 \n');
fprintf(fid,'1.d-11,1.d-11,1.d-11 \n');
fprintf(fid,'SP Final Weights: <1st line real, 2nd line imag> \n');
fprintf(fid,'1.d-11,1.d-11,1.d-11 \n');
fprintf(fid,'1.d-11,1.d-11,1.d-11 \n');
fprintf(fid,'SP start iteration delay (first N iterations do not use soft prior \n');
fprintf(fid,'0 \n');
fclose(fid);


% Build MRE-submit file
viscsubnamev7p3_noSP='MPICH2-Submitv7.3visc';
fid=fopen(viscsubnamev7p3_noSP,'w');
mvfiles(length(mvfiles)+1).name=viscsubnamev7p3_noSP;
fprintf(fid,'#!/bin/bash -l \n');
fprintf(fid,'# declare a name for this job \n'); 
fprintf(fid,'#PBS -N ');fprintf(fid,'%c',viscoutputfv7p3isoincomp_noSP);fprintf(fid,'\n');
fprintf(fid,'# request the queue (enter the possible names, if omitted, serial is the default) \n');
fprintf(fid,'#PBS -q default  \n');
fprintf(fid,'# request  node  \n');
fprintf(fid,'#PBS -l nodes=4:ppn=8 \n');
fprintf(fid,'#PBS -l feature=amd \n');
fprintf(fid,'# request some hours of wall time  \n');
fprintf(fid,'#PBS -l walltime=36:00:00  \n');
fprintf(fid,'#combine PBS standard output and error files  \n');
fprintf(fid,'#PBS -j oe  \n');
fprintf(fid,'# mail is sent to you when the job starts and when it terminates or aborts  \n');
fprintf(fid,'##PBS -m bea  \n');
fprintf(fid,'# specify your email address  \n');
fprintf(fid,'##PBS -M matthew.d.mcgarry@dartmouth.edu  \n');
fprintf(fid,'# Change to Submission Directory  \n');
fprintf(fid,'cd $PBS_O_WORKDIR  \n');
fprintf(fid,'# run the program \n');
fprintf(fid,'cat $PBS_NODEFILE | uniq > node_file \n');
fprintf(fid,'nnodes=$(cat node_file | wc -l) \n');
fprintf(fid,'nprocs=$(cat $PBS_NODEFILE | wc -l) \n');
fprintf(fid,'export MKL_NUM_THREADS=1 \n');
fprintf(fid,'echo Nodes $nnodes \n');
fprintf(fid,'echo Procs $nprocs \n');
fprintf(fid,'mpirun -np $nprocs -hostfile $PBS_NODEFILE /ihome/mmcgarry/code/MREv7p35/MRE-Zone.v7.35.discov ');fprintf(fid,'%c',[viscrfilev7p3_noSP]);fprintf(fid,'\n');
fprintf(fid,'exit 0  \n');
fclose(fid);


%% Move all files to where they should be.
for ii=1:length(mvfiles)
movefile(mvfiles(ii).name,inpath); 
end
tout=toc(touti);
ttotal=toc(t0);
%% Disply time for each part:
disp(['Time for Displacment processing: ' sprintf('%6.2f',tdisp) ' seconds'])
disp(['Time for FE meshing Process: ' sprintf('%6.2f',tfemesh) ' seconds'])
disp(['  Time to build incidence list: ' sprintf('%6.2f',tin) ' seconds'])
disp(['  Time to renumber nodes: ' sprintf('%6.2f',tnodrenum) ' seconds'])
disp(['  Time to find boundary nodes: ' sprintf('%6.2f',tbnod) ' seconds'])
disp(['Time to output files: ' sprintf('%6.2f',tout) ' seconds'])
disp(' ')
disp(['Total processing time: ' sprintf('%6.2f',ttotal) ' seconds'])
end

function [in]=Hex27incidencelist(nodnum)
%Hex27incidencelist:outputs the row of the incidence list from a Hex27
%element built of of a 3x3x3 cube of node numbers, nodnum.

% Nodal coords copied directly from matts masters thesis.
% X(1,:) = [-1,-1,-1 ]; X(10,:) = [-1,-1, 1 ]; X(19,:) = [-1,-1, 0 ];
% X(2,:) = [ 1,-1,-1 ]; X(11,:) = [ 1,-1, 1 ]; X(20,:) = [ 1,-1, 0 ];
% X(3,:) = [ 1, 1,-1 ]; X(12,:) = [ 1, 1, 1 ]; X(21,:) = [ 1, 1, 0 ];
% X(4,:) = [-1, 1,-1 ]; X(13,:) = [-1, 1, 1 ]; X(22,:) = [-1, 1, 0 ];
% X(5,:) = [ 0,-1,-1 ]; X(14,:) = [ 0,-1, 1 ]; X(23,:) = [ 0,-1, 0 ];
% X(6,:) = [ 1, 0,-1 ]; X(15,:) = [ 1, 0, 1 ]; X(24,:) = [ 1, 0, 0 ];
% X(7,:) = [ 0, 1,-1 ]; X(16,:) = [ 0, 1, 1 ]; X(25,:) = [ 0, 1, 0 ];
% X(8,:) = [-1, 0,-1 ]; X(17,:) = [-1, 0, 1 ]; X(26,:) = [-1, 0, 0 ];
% X(9,:) = [ 0, 0,-1 ]; X(18,:) = [ 0, 0, 1 ]; X(27,:) = [ 0, 0, 0 ];
% Plot the nodal coords to check
% for ii=1:27
%     plot3(X(ii,1),X(ii,2),X(ii,3),'r.','markersize',18)
%     title(['node ' int2str(ii) ' added'])
%     grid on
%     hold on
%     pause
% end
% Added 2 to X in matlab to get indices for nodnum 

in=zeros(1,27);

in(1)=nodnum(1,1,1);
in(2)=nodnum(3,1,1);
in(3)=nodnum(3,3,1);
in(4)=nodnum(1,3,1);
in(5)=nodnum(2,1,1);
in(6)=nodnum(3,2,1);
in(7)=nodnum(2,3,1);
in(8)=nodnum(1,2,1);
in(9)=nodnum(2,2,1);
in(10)=nodnum(1,1,3);
in(11)=nodnum(3,1,3);
in(12)=nodnum(3,3,3);
in(13)=nodnum(1,3,3);
in(14)=nodnum(2,1,3);
in(15)=nodnum(3,2,3);
in(16)=nodnum(2,3,3);
in(17)=nodnum(1,2,3);
in(18)=nodnum(2,2,3);
in(19)=nodnum(1,1,2);
in(20)=nodnum(3,1,2);
in(21)=nodnum(3,3,2);
in(22)=nodnum(1,3,2);
in(23)=nodnum(2,1,2);
in(24)=nodnum(3,2,2);
in(25)=nodnum(2,3,2);
in(26)=nodnum(1,2,2);
in(27)=nodnum(2,2,2);
end

%% Thresholded Median filter
function [stackout]=selectivemedianfilter(stackin,mask,thresh)

s=size(stackin);
stackout=stackin;
fsz=1;

for ii=1:s(1)
    %disp(['ii = ' int2str(ii)])
    for jj=1:s(2)
        for kk=1:s(3)
            if(mask(ii,jj,kk)==1)
                nmask=mask( max(ii-fsz,1):min(ii+fsz,s(1)),max(jj-fsz,1):min(jj+fsz,s(2)),max(kk-fsz,1):min(kk+fsz,s(3)) )==1;
                nstack=stackin( max(ii-fsz,1):min(ii+fsz,s(1)),max(jj-fsz,1):min(jj+fsz,s(2)),max(kk-fsz,1):min(kk+fsz,s(3)) );
                medv=median(nstack(nmask));
                if( abs((medv-stackin(ii,jj,kk))/medv)>thresh )
                    stackout(ii,jj,kk)=medv;
                end
            end
        end
    end
end

end
            
            

  

