

function bgmeshobjects(fname); 

close all;


%input is the base of the input file ex: input or meshinginput
%base is the filebase base. (Ex: tumor1)
%stackfn (same thing as mask)

load Mask.mat mask
stackfn=logical(mask);
save stackfn.mat stackfn;

%% Building stackinfo
%load parms.mat parms
load HeaderData DirIndex; 
load parms.mat
% Rows=parms.tags(1,10);
% Columns=parms.tags(1,11);
% PixelSpacing=([parms.tags(1,29),parms.tags(1,30);]);
% SliceThickness=([parms.tags(1,23)+parms.tags(1,24)]);

Rows=size(mask,1); Columns=size(mask,2);
PixelSpacing=[DirIndex(4,1) DirIndex(4,2)];
%SliceThickness=DirIndex(4,3);
SliceThickness=parms.tags(1,23);
SpacingBetweenSlices=parms.tags(1,24); %Slice Gap


stackInfo=struct('Rows',           Rows,...
                  'Columns',        Columns,...
                  'PixelSpacing',   PixelSpacing,...
                  'SliceThickness', SliceThickness,...
                  'SpacingBetweenSlices',SpacingBetweenSlices);

save stackInfo.mat stackInfo;

%% Create Nodes and Elements

load stackfn.mat stackfn
%montagestack(stackfn); drawnow;
load stackInfo.mat stackInfo

isovalue=1;
Nxyz=[2 2 1];

[hiso,hcap,pb]=visibleCube(stackfn,stackInfo,isovalue,Nxyz,{'box',5});

%pb= get(hiso,'vertices');
eb = get(hiso,'faces');

% plot3(pb(:,1)*1e-3, pb(:,2)*1e-3, pb(:,3)*1e-3, '.'); 
% AXIS = ([0 max(pb(:,1))*1e-3 0 max(pb(:,2))*1e-3 0 max(pb(:,3))*1e-3]);
% axis equal; axis(AXIS);
% grid on;
% title('Surface Mesh Nodes'); hold on;


save surfacemesh.mat eb pb

figure(1);
fig.Position=[200 700 1100 400];
subplot (1,2,1)
%plot3(pb(:,1), pb(:,2), pb(:,3), '.'); 
plot3(pb(:,1)*1e-3, pb(:,2)*1e-3, pb(:,3)*1e-3, '.'); 
AXIS = ([0 max(pb(:,1))*1e-3 0 max(pb(:,2))*1e-3 0 max(pb(:,3))*1e-3]);
%AXIS = ([0 max(pb(:,1)) 0 max(pb(:,2)) 0 max(pb(:,3))]);
axis equal; axis(AXIS);
%'Color',[0.5 0 1]
grid on;
title('Surface Mesh Nodes'); hold on;


%% WRITE SURFACE TO VTK FORMAT

fid=fopen(fname,'rt');
    
todo = fgetl(fid);
disp(['Selection --> ',todo]);
todo = str2num(todo);
fprintf('\nEnter the series name: \n');
%read series name
base = fgetl(fid);
disp(['Series Name --> ',base]);


mesh2vtk(eb,pb,[base, '.vtk']);

%% CONVERT SURFACE TO NML FORMAT
fidvt2nml = fopen('convt2nml.inp','wt');
fprintf(fidvt2nml,'%s\n','1');                    % input file format (*.vtk)
fprintf(fidvt2nml,'%s\n','2');                    % Output file formate (*.nml)
fprintf(fidvt2nml,'%s%s\n',base,'.vtk');      % input file name
fprintf(fidvt2nml,'%s%s\n',base,'.nml');      % Output file name
fclose(fidvt2nml);
eval('!convt < convt2nml.inp');


%% Start of Internal Mesh
default=input('Do you want to customize your criteria? 1=yes, 0=no>>:');
if default==1
    facet_angle=input('Facet angle (default=22)>>: '); 
    facet_size=input('Facet size (default=5)>>:');
    facet_distance=input('Facet distance (default=0.8)>>:');
    cell_radius_edge=input('Cell radius edge (default=3)>>:');
    cell_size=input('Cell size (default=2)>>:');
    tumor_cen=input('Tumor Center (default=[100.0000 100.000 20.00])>>:');
    surf_cen=input('Surface Center (default=[100.0000 100.0000 25.000];)>>:');
    R1=input('R1 (default=30)>>:');
    S1=input('S1 (default=3)>>:');
    R2=input('R2 (default=5)>>:');
    S2=input('S2 (default=3)>>:');
    ref_ratio=input('Reference Ratio (default=1)>>:');
    criteria=struct('facet_angle',facet_angle,'facet_size',facet_size, ...
                    'facet_distance',facet_distance, 'cell_radius_edge',cell_radius_edge, ...
                    'cell_size', cell_size, 'tumor_cen',tumor_cen, 'surf_cen',surf_cen, ...
                    'R1',R1, 'S1',S1, 'R2',R2, 'S2',S2', 'ref_ratio',ref_ratio);
                
    
    if isempty(criteria.facet_angle)
        criteria.facet_angle=22;
    end
    
    if isempty(criteria.facet_size)
        criteria.facet_size=5;
    end
    
    if isempty(criteria.facet_distance)
        criteria.facet_distance=0.8;
    end
    
    
    if isempty(criteria.cell_radius_edge)
        criteria.cell_radius_edge=3;
    end
    
    if isempty(criteria.cell_size)
        criteria.cell_size=2; 
    end
    
    
    if isempty(criteria.tumor_cen)
        %criteria.tumor_cen=[117.1875 93.7500 97.5000];
        criteria.tumor_cen=[100.0000 100.000 20.00];

    end
    
    if isempty(criteria.surf_cen)
        %criteria.surf_cen=[110.1174 55.4880 106.7469];
        criteria.surf_cen=[100.0000 100.0000 25.000];
    end
    
    if isempty(criteria.R1)
        criteria.R1=30;
    end
    
    if isempty(criteria.R2)
        criteria.R2=5;
    end

    if isempty(criteria.S1)
        criteria.S1=3;
    end
    
    if isempty(criteria.S2)
        criteria.S2=3;
    end
    
    if isempty(criteria.ref_ratio)
       criteria.ref_ratio=1; 
    end

    criteria
    tic;    
    
    clear eb pb; load surfacemesh.mat
    
    figure(2); 
    plot3(pb(:,1)*1e-3, pb(:,2)*1e-3, pb(:,3)*1e-3, '.'); 
    AXIS = ([0 max(pb(:,1))*1e-3 0 max(pb(:,2))*1e-3 0 max(pb(:,3))*1e-3]);

    axis equal; axis(AXIS);
    
    [elmtemp,nodtemp]=MakeTetraMesh(eb,pb,criteria);

else 
    figure(2); 
    plot3(pb(:,1)*1e-3, pb(:,2)*1e-3, pb(:,3)*1e-3, '.'); 
    AXIS = ([0 max(pb(:,1))*1e-3 0 max(pb(:,2))*1e-3 0 max(pb(:,3))*1e-3]);

    axis equal; axis(AXIS);
    
    tic
    
    [elmtemp,nodtemp]=MakeTetraMesh(eb,pb);

end


%% Generating Boundary Nodes and Elements

nodtemp(:,1)=nodtemp(:,1)*1e-3; 
nodtemp(:,2)=nodtemp(:,2)*1e-3; 
if min(nodtemp(:,3))<0
    nodtemp(:,3)=nodtemp(:,3)*1e-3+min(nodtemp(:,3))*1e-3+(2/3)*DirIndex(4,3)*1e-3;
else
nodtemp(:,3)=nodtemp(:,3)*1e-3+(2/3)*DirIndex(4,3)*1e-3;
end


[belmtemp bnod]=boundfaces(nodtemp,elmtemp);

save meshingoutput elmtemp nodtemp belmtemp bnod

figure(1);
subplot(1,2,2);
plot3(nodtemp(:,1),nodtemp(:,2),nodtemp(:,3),'.b')
hold on;
plot3(bnod(:,1),bnod(:,2),bnod(:,3),'.r');
title('All Nodes');
xlabel('X'); ylabel('Y'); zlabel('Z')
grid on; AXIS = ([0 max(nodtemp(:,1)) 0 max(nodtemp(:,2)) 0 max(nodtemp(:,3))]);
axis equal; axis(AXIS);
legend('Interior Nodes','Boundary Nodes','Location','Best')


%% Writing everything to files

if ~isempty(elmtemp);
    disp('There are elements!');
    elm=zeros(size(elmtemp,1),6);
    elm(:,1)=1:1:size(elmtemp,1);
    elm(:,2)=elmtemp(:,1); elm(:,3)=elmtemp(:,2); elm(:,4)=elmtemp(:,3); 
    elm(:,5)=elmtemp(:,4); elm(:,6)=elmtemp(:,5);

else 
    elm=elmtemp;
    warning('There are no elements! Remesh using different criteria.');
end


nod=zeros(size(nodtemp,1),4);
nod(:,1)=1:1:size(nodtemp,1);
nod(:,2)=nodtemp(:,1); 
nod(:,3)=nodtemp(:,2); 
nod(:,4)=nodtemp(:,3);

% Writing Nodes to File
nodoutf=[base '.nod'];
fidnodout=fopen(nodoutf,'w');
fprintf(fidnodout, '%7i %15.8e %15.8e %15.8e \n',nod');
fclose(fidnodout);

% Writing Elements to File
elmoutf=[base,'.elm'];
fidelmout=fopen(elmoutf,'w');
fprintf(fidelmout, '%7i %7i %7i %7i %7i %7i \n',elm');
fclose(fidelmout);

if ~isempty(belmtemp);
    disp('There are boundary elements!');
    belm_ones=(ones(size(belmtemp,1),1))';
    belm_zeros=(zeros(size(belmtemp,1),1))';
    belm_nn=(1:1:size(belmtemp,1))';

    bel(:,1)=belm_nn;
    bel(:,2)=belmtemp(:,1);
    bel(:,3)=belmtemp(:,2);
    bel(:,4)=belmtemp(:,3);
    bel(:,5)=belm_ones;
    bel(:,6)=belm_zeros;

    % Writing Boundary Elements to File
    belf=[base '.bel'];
    fidbelf=fopen(belf,'w');
    fprintf(fidbelf,'%7i %7i %7i %7i %7i %7i \n',bel');
    fclose(fidbelf);
else 
    warning('There does not seem to be any boundary elements')
end


if ~isempty(bnod);
    disp('There are boundary nodes!');
    for ii=1:size(bnod,1)
        bnod(ii,4)=find(bnod(ii,1)==nodtemp(:,1) & bnod(ii,2)==nodtemp(:,2) & bnod(ii,3)==nodtemp(:,3));
    end

    bnod_nn=1:1:size(bnod,1);
    bnodfin(:,1)=bnod_nn';
    bnodfin(:,2)=bnod(:,4);

    % Writing Boundary Nodes to File
    bnodf=[base,'.bnod'];
    fidbnodf=fopen(bnodf,'w');
    fprintf(fidbnodf,'%7i %7i \n',bnodfin');
    fclose(fidbnodf)
else 
    warning('There does not seem to be any boundary nodes')
    
end




%% Interpolate displacements from grid to the FE mesh
load MRE_3DMotionData A P MagIm% variables A, P, MagIm
load HeaderData

A(isnan(A))=0;
P(isnan(P))=0;

fid1=fopen('freq.tmp','wt');

%w=create new file for writing
%t=open file in text mode

fprintf(fid1,'%.4f',freqHz);

%gets 4 places after the decimal point from fid1

fclose(fid1);

load HeaderData DirIndex %variable DirIndex
load Mask.mat mask
 
% gen3ddisp(nod, A,P,MagIm,DirIndex,base,mask); %Internal Function
% 
% eval(sprintf(['! awk ''{print $1,0,0}'' < ',base,'.nod > pressure.s1c']));
% 
% %generating file containing initial guesses
% eval(sprintf(['! awk ''{print $1,1,1,1}'' < ',base,'.nod > ',base,'.inv.mtr']));


    % INTERPOLATE DISPLACEMENT DATA ONTO FINITE ELEMENT MESH, GENERATE
    % RECONSTRUCTION ALGORITHIM input FILES, AND LAUNCH MRE
    % RECONSTRUCTION (PARALLEL VERSION)
    
    
    fprintf('Enter Mesh Resolution: \n');
    res = fgetl(fid);
    disp(['Mesh Resolution --> ',res]);
    
    fprintf('\nEnter your email address: \n');
    %read email address
    add = fgetl(fid);
    disp(['Email Address --> ',add]);

    freq = load('freq.tmp');

    fprintf('\nEnter spatial filter weight [0-1]: \n');
    %read spatial filter weight
    spfilt = fgetl(fid);
    disp(['Spatial Filter Weight --> ',spfilt]);

    fprintf('\nEnter maximum number of iterations: \n');
    %read max iterations
    itmax = fgetl(fid);
    disp(['Maximum Iterations --> ',itmax]);

    fprintf('\nEnter regularization parameter: \n');
    %read regularization parameter
    reg = fgetl(fid);
    disp(['Regularization Parameter --> ',reg]);

    fprintf('\nEnter number of processors: \n');
    %read number of processors
    nprocs = fgetl(fid);
    disp(['Number of Processors --> ',nprocs]);

    %create INV directories if they do not already exist
    !mkdir -p INV
    !mkdir -p INVpar

    %ELASTIC RECONSTRUCTION
    fprintf('\nGENERATING 3D_MRE_FILES.DAT\n');
    dir = 'INVpar';
    udopt = 2;
    znopt = 1;
    diag = 0.02;
    ovrlp = 0.10;
    gen_MREpar_reconfile(base,dir,itmax,freq,reg,spfilt,udopt,znopt,diag,ovrlp,0);

    fprintf('\nPERFORMING DATA INTERPOLATION\n');
    %interpolate data onto finite element mesh using MATLAB scripts

    nod = load([base,'.nod']);
    %interpolate measured displacements onto finite element mesh
    gen3ddisp(nod,A,P,MagIm,DirIndex,base,mask);
    
    %generate pressure file
    eval(sprintf(['! awk ''{print $1,0,0}'' < ',base,'.nod > pressure.s1c']));

    %generate complex input files  

    cmplx_input(base);

    
    %generate file containing initial guess
    eval(sprintf(['! awk ''{print $1,1,1,1}'' < ',base,'.nod > ',base,'.inv.mtr']));
    
    % POROELASTIC RECONSTRUCTION
    fprintf('\nGENERATING 3D_MRPE_FILES.DAT\n');
    dir = 'INV';
    udopt = 1;
    znopt = 1;
    diag = 0.02;
    ovrlp = 0.10;
    gen_MRPE_reconfile(base,dir,itmax,freq,reg,spfilt,udopt,znopt,diag,ovrlp,0);
    gen_MRPEkappa_reconfile(base,dir,itmax,freq,reg,spfilt,udopt,znopt,diag,ovrlp,0)
    
    %generate discovery cluster submit-files
    gen_runfiles(add,nprocs,base);
    
    softprior(base,2,0,0)
    
    ! mkdir -p MRPE
    ! mkdir -p MRPE/INV
    ! mv *.ci pressure.s1c 3D_MRPE* MRPE-* MRPE_*  MRPE/ 

    % MAP RECONSTRUCTED ELASTIC PARAMETERS BACK TO MR SPACE
    fprintf('\nGENERATING VOL_INTERP.INP\n');
    %generate vol_interp.inp file
    [mp,np,zp]=size(MagIm);
    fid1=fopen('vol_interp.inp','wt');
    fprintf(fid1,'%i %i %i\n',mp,np,zp+2);
    fprintf(fid1,'%e %e %e\n',DirIndex(4,1:3)'*1e-3);
    fprintf(fid1,[base,'.nod.header\n']);
    fprintf(fid1,[base,'.elm.header\n']);
    fprintf(fid1,'    1\n');
    fprintf(fid1,['2   INV/',base,'.inv.',itmax,'\n']);
    fprintf(fid1,['2   INVpar/',base,'.inv.',itmax,'\n']);
    fprintf(fid1,['2   MRPE/INV/',base,'.inv.',itmax,'\n']);
    fclose(fid1);

    disp(['base = ' base]);
    % Generate strain data.
    fprintf('\nCALCULATING STRAINS:\n');

    %Generate strain runfile
    fidstr=fopen('strainrunfile.inp','w');
    fprintf(fidstr,[base '.nod\n']);
    fprintf(fidstr,[base '.elm\n']);
    fprintf(fidstr,[base '.dsp\n']);
    fclose(fidstr);

    %run strain calculator
    %strex='/ihome/lsolamen/code/misc/Strain_calc_tet.x';  % base of strain calculating executable
    %eval(['!' strex ' strainrunfile.inp']);

    %run strain calculator
    strex='/home/apattiso/bin/Strain_calc_tet_auto.x';  % Name of strain calculating executable
    eval(['!' strex ' strainrunfile.inp']);

    %create auxillary file directories if they do not already exist
    ! mkdir -p INP
    ! mkdir -p LOG
    ! mkdir -p TMP
    ! mkdir -p DIA
    
    %move all auxillary files to their respective directories
    ! mv *.inp INP
    ! mv *.log LOG
    ! mv *.tmp TMP
    ! mv *.dia DIA
    
    %remove files
    ! rm -r -f INP
    ! rm -r -f LOG
    ! rm -r -f TMP
    ! rm -r -f DIA

    eval(sprintf(['! rm ',base,'.vtk']));
    eval(sprintf(['! rm ',base,'.nml'])); 
    
    
%% Generating Region Stack File (Adapted from MRIhexmex_v7p3)

vox=DirIndex(4,1:3).*1e-3;
dim=size(MagIm);
xin=0:vox(1):vox(1)*(dim(1)-1);
yin=0:vox(2):vox(2)*(dim(2)-1);
zin=0:vox(3):vox(3)*(dim(3)-1);

    
mridim=size(mask);

regoutf=[base '.regstack'];
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
    dlmwrite(regoutf,mask(:,:,ii),'-append','delimiter',' ');
end    
   
disp('Mesh Created!');

disp(sprintf('Number of Boundary Nodes: %d \r',size(bnod,1))); 
disp(sprintf('Number of Nodes: %d \r',size(nod,1)));
disp(sprintf('Number of Internal Nodes: %d \r',size(nod,1)-size(bnod,1)))

disp('Generating Run File');
genrunfilev8poro(base,itmax);


!mv runfile_v8poro.dat MRPE/
!mv MPICH2-Submitv8poro MRPE/

disp(sprintf('REMEMBER TO CALCULATE INITIAL BOUNDARY CONDITIONS!')); 


toc
end


%% Internal Function gen3ddisp   
function gen3ddisp(nod,A,P,MagIm,DirIndex,base,mask);
%ooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooo
%   Adapted from MREmesh2.m
%   GEN3DDISP.M
%     This routine is used to compute and interpolate the measured
%     displacement field onto the finite element mesh
%
%   inputS
%     mask        [nx,ny,nz] array of segmented binary images outlining
%                   the region of interest to be meshed
%     stackInfo   a structure containing the pixel-spacing and distance
%                   between slices for the image stack
%     base        series base
%ooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooo

% pixel spacing
dx=DirIndex(4,1)*1e-3;
dy=DirIndex(4,2)*1e-3;
dz=DirIndex(4,3)*1e-3;

% displacements - stack without top&bottom slices
Dr=1e-6.*A.*cos(P);
Di=1e-6.*A.*sin(P);
[nx, ny, nslice, nd]=size(Dr);
clear A P;

% 'pad' displ
pDr=zeros(nx,ny,nslice+2,nd);
pDi=zeros(nx,ny,nslice+2,nd);

pDr(:,:,1,:)=Dr(:,:,1,:);
pDr(:,:,2:nslice+1,:)=Dr;
pDr(:,:,nslice+2,:)=Dr(:,:,end,:);

pDi(:,:,1,:)=Di(:,:,1,:);
pDi(:,:,2:nslice+1,:)=Di;
pDi(:,:,nslice+2,:)=Di(:,:,end,:);

% % Dilate motions around edges a few times to reduce edge interpolation warnings
disp('Dilating to avoid edge warnings')
ndil=3;
mskpad=zeros(nx,ny,nslice+2);
mskpad(:,:,2:nslice+1,:)=mask;

for ii=1:ndil
    for jj=1:3
        maskdisp=((pDr(:,:,:,jj)~=0)|mskpad);
        [pDr(:,:,:,jj)]=dilatearray3D(pDr(:,:,:,jj),maskdisp);
        [pDi(:,:,:,jj)]=dilatearray3D(pDi(:,:,:,jj),maskdisp);         
    end
end


%[xx,yy,zz]=(dx:dx:(nx)*dx,dy:dy:(ny)*dy,0:dz:(nslice)*dz);
[xx,yy,zz]=meshgrid(dx:dx:(nx)*dx,dy:dy:(ny)*dy,0:dz:(nslice+1)*dz);

%% Interpolate displacements from grid to the FE mesh

% real displacements
displ(:,1)=nod(:,1);
displ(:,2)=interp3(xx,yy,zz,pDr(:,:,:,1),nod(:,2),nod(:,3),nod(:,4),'spline');
displ(:,3)=interp3(xx,yy,zz,pDr(:,:,:,2),nod(:,2),nod(:,3),nod(:,4),'spline');
displ(:,4)=interp3(xx,yy,zz,pDr(:,:,:,3),nod(:,2),nod(:,3),nod(:,4),'spline');

% imaginary displacements
cdispl(:,1)=nod(:,1);
cdispl(:,2)=interp3(xx,yy,zz,pDi(:,:,:,1),nod(:,2),nod(:,3),nod(:,4),'spline');
cdispl(:,3)=interp3(xx,yy,zz,pDi(:,:,:,2),nod(:,2),nod(:,3),nod(:,4),'spline');
cdispl(:,4)=interp3(xx,yy,zz,pDi(:,:,:,3),nod(:,2),nod(:,3),nod(:,4),'spline');

%%Transform Nodal Coordinates according with DirIndex directions
nod_transform(nod,displ,cdispl,DirIndex,base);

end

%% Nodal Transform
function nod_transform(nod,displ,cdispl,DirIndex,base)
%Adapted from MREmesh2    
%ooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooo
%   NOD_TRANSFORM.M
%     This routine is used to transform the nodal co-ordinate and
%     displacement data using DirIndex
%
%   inputS
%     nod         [nn,4] array containing nodal coordiantes
%     displ       [nn,3] array containing real-valued displacements
%     cdispl      [nn,6] array containing complex-valued displacements
%     DirIndex    4 X 6 matrix containing information regarding how to
%                   transform displacement data and nodal co-ordinates
%     base        series base
%ooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooo    

% Determine Column Locations for Transformation Matrices
Index_1=find(abs(DirIndex(1,:))==1);
Index_2=find(abs(DirIndex(2,:))==1);
Index_3=find(abs(DirIndex(3,:))==1);

% Get Voxel Spacing
dx = DirIndex(4,1)*1e-3;
dy = DirIndex(4,2)*1e-3;
dz = DirIndex(4,3)*1e-3;

% new_nod(:,1) = nod(:,1);
% new_nod(:,2) = nod(:,Index_1(1)+1)*DirIndex(1,Index_1(1));
% new_nod(:,3) = nod(:,Index_2(1)+1)*DirIndex(2,Index_2(1));
% new_nod(:,4) = nod(:,Index_3(1)+1)*DirIndex(3,Index_3(1));    
   
% Displacement Transformations
% Real Component
new_disp(:,1) = displ(:,1);
new_disp(:,2) = displ(:,Index_1(2)-3+1)*DirIndex(1,Index_1(2));
new_disp(:,3) = displ(:,Index_2(2)-3+1)*DirIndex(2,Index_2(2));
new_disp(:,4) = displ(:,Index_3(2)-3+1)*DirIndex(3,Index_3(2));

% Imaginary Component
new_cdisp(:,1) = cdispl(:,1);
new_cdisp(:,2) = cdispl(:,Index_1(2)-3+1)*DirIndex(1,Index_1(2));
new_cdisp(:,3) = cdispl(:,Index_2(2)-3+1)*DirIndex(2,Index_2(2));
new_cdisp(:,4) = cdispl(:,Index_3(2)-3+1)*DirIndex(3,Index_3(2));

% Write Transformed Displacement Vectors to File
fid=fopen([base,'.v3r'],'wt');
fprintf(fid,'%6d  %24.16E  %24.16E  %24.16E\n',new_disp');
fclose(fid);

temp = [new_disp(:,1:2) new_cdisp(:,2) new_disp(:,3) new_cdisp(:,3) new_disp(:,4) new_cdisp(:,4)];
% Write Transformed Complex Valued Displacements to File
fid=fopen([base,'.dsp'],'wt');
fidv3c=fopen([base,'.v3c'],'wt'); %v3c is the same as dsp, just different tags
fprintf(fid,'%6d  %24.16E  %24.16E  %24.16E  %24.16E  %24.16E  %24.16E\n',temp');
fprintf(fidv3c,'%6d  %24.16E  %24.16E  %24.16E  %24.16E  %24.16E  %24.16E\n',temp');
fclose(fid); fclose(fidv3c);

    
end

%% 
%ooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooo
%   CMPLX_input
%     This routine is used to rewrite the complex displacement and
%     pore-pressure data files such that they can be read into FORTRAN
%     as a complex variables
%
%   inputS
%     base        series base
%
%   Phillip R. Perriñez, Ph.D.
%   Thayer School of Engineering
%   Dartmouth College
%   July 2009
%ooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooo

function cmplx_input(base)

%load .nod file
nod = load([base,'.nod']);
nn = size(nod,1);

%load complex displacement file:
UC = load([base,'.dsp']);
PC = load('pressure.s1c');

disp_output = [base,'.dsp.ci'];
press_output = 'pressure.s1c.ci';

Uxr = UC(:,2);  Uxi = UC(:,3);
Uyr = UC(:,4);  Uyi = UC(:,5);
Uzr = UC(:,6);  Uzi = UC(:,7);

Pr = PC(:,2);   Pi = PC(:,3);

fid=fopen(disp_output,'wt');
d = [Uxr Uxi Uyr Uyi Uzr Uzi];
fprintf(fid,'%6d  (%24.16E,%24.16E)  (%24.16E,%24.16E)  (%24.16E,%24.16E)\n',[1:nn; d']);
fclose(fid);

fid=fopen(press_output,'wt');
p = [Pr Pi];
fprintf(fid,'%6d  (%24.16E,%24.16E)\n',[1:nn; p']);
fclose(fid);
end

%% Dilate Array 3D

function [arrayout]=dilatearray3D(arrayin,mask)
% Dilate the edge values in an array to avoid warnings at the boundary during
% interpolation

s=size(arrayin);

arrayout=arrayin;

mask=double(mask);
if(max(mask(:))~=1)
    warning('Dilate3d: Max(mask) ~=1')
end
    

for ii=1:s(1)
    for jj=1:s(2)
        for kk=1:s(3)
            if(mask(ii,jj,kk)~=1)
                msk=mask(max(1,ii-1):min(s(1),ii+1),max(1,jj-1):min(s(2),jj+1),max(1,kk-1):min(s(3),kk+1));
                if(sum(msk(:)>0))
                    ary=arrayin(max(1,ii-1):min(s(1),ii+1),max(1,jj-1):min(s(2),jj+1),max(1,kk-1):min(s(3),kk+1));
                    arrayout(ii,jj,kk)=sum(msk(:).*ary(:))/sum(msk(:));
                end
            end
        end
    end
end

end
  



