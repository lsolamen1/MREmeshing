%ooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooo
%
%   MREMESH.M
%
%   DESCRIPTION
%     This script is used to automate the MRE preprocessing chain
%
%   INPUTS
%     fname     name of input file
%
%   Phillip R. Perriñez, Ph.D.
%   Thayer School of Engineering
%   Dartmouth College
%   July 2009
%ooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooo
function MREmesh2(fname)

close all;
tic
fid=fopen(fname,'rt');

fprintf('-------------------------------------------------------------\n');
fprintf('\n');
fprintf('\tTime-Harmonic M?.R Elastographic Imaging\n\n');
fprintf('-------------------------------------------------------------\n\n');

fprintf('\nChoose from one of the following options:\n');
fprintf(' [1]  mesh generation only\n');
fprintf(' [2]  data interpolation only\n');
fprintf(' [3]  mesh generation and data interpolation\n');
%read option number
todo = fgetl(fid);
disp(['Selection --> ',todo]);
todo = str2num(todo);
fprintf('\nEnter the series name: \n');
%read series name
name = fgetl(fid);
disp(['Series Name --> ',name]);

fprintf('\n-------------------------------------------------------------\n');
fprintf('USING MATLAB FILES\n');
fname1='HeaderData.mat';
fname2='MRE_3DMotionData.mat';
fname3='Mask.mat';

disp([fname1,'  ',fname2,'  ',fname3]);

%load matlab files
load(fname1);
load(fname2);
load(fname3);

A(isnan(A))=0;
P(isnan(P))=0;


%write MRE frequency to file
fid1=fopen('freq.tmp','wt');
%w=create new file for writing
%t=open file in text mode
fprintf(fid1,'%.4f',freqHz);
%gets 4 places after the decimal point from fid1
fclose(fid1);

%read mesh resolution
fprintf('Enter Mesh Resolution: \n');
res = fgetl(fid);
disp(['Mesh Resolution --> ',res]);

%PERFORM MESH GENERATION
if todo == 1 || todo == 3
    %todo is mesh generation option

    %retain pixel spacing
    stackInfo.PixelSpacing = DirIndex(4,1:2)'*1e-3;
    stackInfo.SpacingBetweenSlices = DirIndex(4,3)'*1e-3;
    fprintf('\nGENERATING SURFACE MESH\n');
    disp('MRE_surfacemesh.m')
    MRE_surfacemesh(mask,stackInfo,name);
    %MRE_surfacemesh is an internal function (look below)

    %generate SPMESH run file
    fid1=fopen('spmesh.inp','wt');
    fprintf(fid1,'1 \n');
    fprintf(fid1,'0 \n');
    fprintf(fid1,[name,'.nml \n']);
    fprintf(fid1,'1 \n');
    fprintf(fid1,[res,'\n']);
    fprintf(fid1,'0 \n');
    fprintf(fid1,'0 \n');
    fprintf(fid1,'1 \n');
    fprintf(fid1,'1 \n');
    fprintf(fid1,'1 \n');
    fprintf(fid1,'0 \n');
    fprintf(fid1,[name,'\n']);
    fprintf(fid1,'-1 \n');
    fclose(fid1);
    
    fprintf('\nGENERATING VOLUME MESH\n');
    ! spmesh < spmesh.inp
  
    eval(['!cp ' name '.elm ' name '.elm.header']);
    eval(['!cp ' name '.nod ' name '.nod.header']);
    
    % generate 3DTRY run file
    fid1=fopen('3dtry.inp','wt');
    fprintf(fid1,'1 \n');
    fprintf(fid1,[name,'.elm.header\n']);
    fprintf(fid1,[name,'.nod.header\n']);
    fprintf(fid1,'0 \n');
    fprintf(fid1,'0 \n');
    fprintf(fid1,'1 \n');
    fprintf(fid1,[name,'.bel\n']);
    fclose(fid1);
    
    fprintf('\nGENERATING BOUNDARY ELEMENTS\n')
    ! 3dtry < 3dtry.inp
    
    fprintf('\n');
    eval(['!cp ' name '.bel ' name '.bel.header']);
    
    %strip file headers
    eval(sprintf(['! sed ''1d'' < ',name,'.elm > ',name,'.elm.tmp']));
    eval(sprintf(['! sed ''1d'' < ',name,'.nod > ',name,'.nod.tmp']));
    eval(sprintf(['! sed ''1d'' < ',name,'.bel > ',name,'.bel.tmp']));
    
    %rename files
    movefile([name,'.elm.tmp'],[name,'.elm']);
    movefile([name,'.nod.tmp'],[name,'.nod']);
    movefile([name,'.bel.tmp'],[name,'.bel']);
    
    %generate BNODGEN run file
    fid1=fopen('bnodgen.inp','wt');
    fprintf(fid1,[name,'.bel\n']);
    fprintf(fid1,[name,'.bnod\n']);
    fclose(fid1);
    
    fprintf('\nGENERATING BOUNDARY NODE INDEX\n');
    ! bnodgen < bnodgen.inp
    
    fprintf('\n****** MESH GENERATION COMPLETE ******\n\n');
end
    
if (todo == 2 || todo == 3)
    % INTERPOLATE DISPLACEMENT DATA ONTO FINITE ELEMENT MESH, GENERATE
    % RECONSTRUCTION ALGORITHIM INPUT FILES, AND LAUNCH MRE
    % RECONSTRUCTION (PARALLEL VERSION)
    
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
    gen_MREpar_reconfile(name,dir,itmax,freq,reg,spfilt,udopt,znopt,diag,ovrlp,0);

    fprintf('\nPERFORMING DATA INTERPOLATION\n');
    %interpolate data onto finite element mesh using MATLAB scripts

    nod = load([name,'.nod']);
    %interpolate measured displacements onto finite element mesh
    gen3ddisp(nod,A,P,MagIm,DirIndex,name);
    
    %generate pressure file
    eval(sprintf(['! awk ''{print $1,0,0}'' < ',name,'.nod > pressure.s1c']));

    %generate complex input files  
    cmplx_input(name);

    
    %generate file containing initial guess
    eval(sprintf(['! awk ''{print $1,1,1,1}'' < ',name,'.nod > ',name,'.inv.mtr']));
    
    % POROELASTIC RECONSTRUCTION
    fprintf('\nGENERATING 3D_MRPE_FILES.DAT\n');
    dir = 'INV';
    udopt = 1;
    znopt = 1;
    diag = 0.02;
    ovrlp = 0.10;
    gen_MRPE_reconfile(name,dir,itmax,freq,reg,spfilt,udopt,znopt,diag,ovrlp,0);
    gen_MRPEkappa_reconfile(name,dir,itmax,freq,reg,spfilt,udopt,znopt,diag,ovrlp,0)
    
    %generate discovery cluster submit-files
    gen_runfiles(add,nprocs,name);
    
    softprior(name,2,0,0)
    
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
    fprintf(fid1,[name,'.nod.header\n']);
    fprintf(fid1,[name,'.elm.header\n']);
    fprintf(fid1,'    1\n');
    fprintf(fid1,['2   INV/',name,'.inv.',itmax,'\n']);
    fprintf(fid1,['2   INVpar/',name,'.inv.',itmax,'\n']);
    fprintf(fid1,['2   MRPE/INV/',name,'.inv.',itmax,'\n']);
    fclose(fid1);

    disp(['NAME = ' name]);
    % Generate strain data.
    fprintf('\nCALCULATING STRAINS:\n');

    %Generate strain runfile
    fidstr=fopen('strainrunfile.inp','w');
    fprintf(fidstr,[name '.nod\n']);
    fprintf(fidstr,[name '.elm\n']);
    fprintf(fidstr,[name '.v3c\n']);
    fclose(fidstr);

    %run strain calculator
    strex='/home/apattiso/bin/Strain_calc_tet_auto.x';  % Name of strain calculating executable
    eval(['!' strex ' strainrunfile.inp']);

end % END IF TODO == 2

if (todo == 1 || todo == 2 || todo == 3)
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
    ! rm *.out *.header
    eval(sprintf(['! rm ',name,'.vtk']));
    eval(sprintf(['! rm ',name,'.nml'])); 
end
fprintf('\n\n\tDONE\n');
toc
end % END MAIN FUNCTION
%ooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooo
%ooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooo

%% INTERNAL FUNCTIONS

function MRE_surfacemesh(mask,stackInfo,name)
%ooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooo
%   MRE_SURFACEMESH.M
%       This routine is used to generate a smooth surfacemesh to be used
%       as a guide when generating the volume mesh
%
%   INPUTS
%       mask        [nx,ny,nz] array of segmented binary images outlining
%                   the region of interest to be meshed
%       stackInfo   a structure containing the pixel-spacing and distance
%                   between slices for the image stack
%       name        series name
%ooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooo

% SCALE INPUT
[rows,cols,slices] = size(mask);

% PAD IMAGE STACK TO OBTAIN CLOSED SURFACE MESH
A = zeros(rows,cols,slices+2);
A(:,:,2:end-1) = mask;

% SMOOTH THE DATA
Ds = smooth3(A);

% MESH
% an isovalue of 0.67 creates a surface mesh within 0.01 distance to the
% first and last images
% an isovalue of 0.68 creates a smaller surface mesh within 0.04 distance
% to the first and last images

isoV = 0.6; %0.68;
figure
hiso = patch(isosurface(Ds,isoV),'FaceColor',[0,0.5,0],'EdgeColor','none');
set(hiso, 'tag', 'isosurface');

lightangle(45,30);
set(gcf,'Renderer','zbuffer'); lighting phong;
isonormals(Ds,hiso);
set(hiso,'SpecularColorReflectance',0,'SpecularExponent',50);

% SCALE DATA TO RETAIN CORRECT SIZE
% shifting Z coordinate of nodes to compensate for the shift introduced by
% the zero padding of the mask
nod = get(hiso,'vertices');
if ~isempty(nod)
    set(hiso,'vertices',[ nod(:,1)*stackInfo.PixelSpacing(1), ...
        nod(:,2)*stackInfo.PixelSpacing(2), ...
        (nod(:,3)-1)*stackInfo.SpacingBetweenSlices ]);
end
view(3); grid;
%view(45,30);  axis tight; axis equal; grid;

nod = get(hiso,'vertices');
elm = get(hiso,'faces');

%% WRITE SURFACE TO VTK FORMAT
mesh2vtk(elm,nod,[name, '.vtk']);

%% CONVERT SURFACE TO NML FORMAT
fid = fopen('convt2nml.inp','wt');
fprintf(fid,'%s\n','1');                    % Input file format (*.vtk)
fprintf(fid,'%s\n','2');                    % Output file formate (*.nml)
fprintf(fid,'%s%s\n',name,'.vtk');      % Input file name
fprintf(fid,'%s%s\n',name,'.nml');      % Output file name
fclose(fid);
eval('!convt < convt2nml.inp');
end

function gen3ddisp(nod,A,P,MagIm,DirIndex,name)
%ooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooo
%   GEN3DDISP.M
%     This routine is used to compute and interpolate the measured
%     displacement field onto the finite element mesh
%
%   INPUTS
%     mask        [nx,ny,nz] array of segmented binary images outlining
%                   the region of interest to be meshed
%     stackInfo   a structure containing the pixel-spacing and distance
%                   between slices for the image stack
%     name        series name
%ooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooo
load('Mask.mat');
% pixel spacing
dx=DirIndex(4,1)*1e-3;
dy=DirIndex(4,2)*1e-3;
dz=DirIndex(4,3)*1e-3;

% displacements - stack without top&bottom slices
Dr=1e-6.*A.*cos(P);
Di=1e-6.*A.*sin(P);
[nx ny nslice nd]=size(Dr);
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

% % Dilate motions around edges a few times to reduce edge interpolation errors
disp('Dilating to avoid edge errors')
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

%% Create mesh grid to interpolate displacements and display
% z dir includes the 'extra' slices in the nslice dimension
[xx,yy,zz]=meshgrid(dx:dx:(nx)*dx,dy:dy:(ny)*dy,0:dz:(nslice+1)*dz);
% [yy,xx,zz]=ndgrid(dx:dx:(nx)*dx,dy:dy:(ny)*dy,0:dz:(nslice+1)*dz);

MINx = min(xx(:));MAXx = max(xx(:));
MINy = min(yy(:));MAXy = max(yy(:));
MINz = min(zz(:));MAXz = max(zz(:));

dim = [MINx MINy MINz; ...
    MINx MAXy MINz; ...
    MAXx MINy MINz; ...
    MAXx MAXy MINz; ...
    MINx MINy MAXz; ...
    MINx MAXy MAXz; ...
    MAXx MINy MAXz; ...
    MAXx MAXy MAXz];

plot3(dim(:,1),dim(:,2),dim(:,3),'.r');
grid on;
hold on;

%% Plot MRE First & Last Images
MagIm = MagIm(:,:,:,1).*mask;

X=[1,nx;1,nx]*dx;
Y=[1,1;ny,ny]*dy;
Z=[1,1;1,1]*dz;
surface(X,Y,Z,'CData',double(MagIm(:,:,1)),'FaceColor','texturemap');
colormap(gray);drawnow;

Z=[1,1;1,1]*dz*(nslice);
surface(X,Y,Z,'CData',double(MagIm(:,:,end)),'FaceColor','texturemap');
colormap(gray);drawnow;

%% Display FE mesh
bel = load([name,'.bel']);
trisurf(bel(:,2:4),nod(:,2),nod(:,3),nod(:,4),'facecolor','g','EdgeColor','none');
axis equal;xlabel('\bfX'),ylabel('\bfY'),zlabel('\bfZ');

%% Interpolate displacements from grid to the FE mesh

% real displacements
displ(:,1)=nod(:,1);
displ(:,2)=interp3(xx,yy,zz,pDr(:,:,:,1),nod(:,2),nod(:,3),nod(:,4),'spline');
displ(:,3)=interp3(xx,yy,zz,pDr(:,:,:,2),nod(:,2),nod(:,3),nod(:,4),'spline');
displ(:,4)=interp3(xx,yy,zz,pDr(:,:,:,3),nod(:,2),nod(:,3),nod(:,4),'spline');

% imaginary displlacements
cdispl(:,1)=nod(:,1);
cdispl(:,2)=interp3(xx,yy,zz,pDi(:,:,:,1),nod(:,2),nod(:,3),nod(:,4),'spline');
cdispl(:,3)=interp3(xx,yy,zz,pDi(:,:,:,2),nod(:,2),nod(:,3),nod(:,4),'spline');
cdispl(:,4)=interp3(xx,yy,zz,pDi(:,:,:,3),nod(:,2),nod(:,3),nod(:,4),'spline');

%% Transform Nodal Coordinates according with DirIndex directions
nod_transform(nod,displ,cdispl,DirIndex,name);
end

function nod_transform(nod,displ,cdispl,DirIndex,name)
%ooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooo
%   NOD_TRANSFORM.M
%     This routine is used to transform the nodal co-ordinate and
%     displacement data using DirIndex
%
%   INPUTS
%     nod         [nn,4] array containing nodal coordiantes
%     displ       [nn,3] array containing real-valued displacements
%     cdispl      [nn,6] array containing complex-valued displacements
%     DirIndex    4 X 6 matrix containing information regarding how to
%                   transform displacement data and nodal co-ordinates
%     name        series name
%ooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooo

% Determine Column Locations for Transformation Matrices
Index_1=find(abs(DirIndex(1,:))==1);
Index_2=find(abs(DirIndex(2,:))==1);
Index_3=find(abs(DirIndex(3,:))==1);

% Get Voxel Spacing
dx = DirIndex(4,1)*1e-3;
dy = DirIndex(4,2)*1e-3;
dz = DirIndex(4,3)*1e-3;

% Nodal Co-ordinate transformation
% new_nod(:,1)=nod(:,1);
% new_nod(:,2)=nod(:,2);
% new_nod(:,3)=nod(:,3);
% new_nod(:,4)=nod(:,4);

new_nod(:,1) = nod(:,1);
new_nod(:,2) = nod(:,Index_1(1)+1)*DirIndex(1,Index_1(1));
new_nod(:,3) = nod(:,Index_2(1)+1)*DirIndex(2,Index_2(1));
new_nod(:,4) = nod(:,Index_3(1)+1)*DirIndex(3,Index_3(1));

% Write Transformed Nodal Coordinates to File
fid=fopen([name,'.nod.new.header'],'wt');
np=length(new_nod);
fprintf(fid,'%i\t%i\n',1,np);
fprintf(fid,'%i\t%e\t%e\t%e\n',new_nod');
fclose(fid);
%-------------------------------------------------------------------------

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
fid=fopen([name,'.v3r'],'wt');
fprintf(fid,'%6d  %24.16E  %24.16E  %24.16E\n',new_disp');
fclose(fid);

temp = [new_disp(:,1:2) new_cdisp(:,2) new_disp(:,3) new_cdisp(:,3) new_disp(:,4) new_cdisp(:,4)];
% Write Transformed Complex Valued Displacements to File
fid=fopen([name,'.v3c'],'wt');
fprintf(fid,'%6d  %24.16E  %24.16E  %24.16E  %24.16E  %24.16E  %24.16E\n',temp');
fclose(fid);

bnod=load([name,'.bnod']);
np=length(bnod);
fid=fopen([name,'.bnod.header'],'wt');
fprintf(fid,'%i %i\n',1,np);
fprintf(fid,'%i %i\n',bnod');
fclose(fid);
end

function cmplx_input(name)
%ooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooo
%   CMPLX_INPUT
%     This routine is used to rewrite the complex displacement and
%     pore-pressure data files such that they can be read into FORTRAN
%     as a complex variables
%
%   INPUTS
%     name        series name
%
%   Phillip R. Perriñez, Ph.D.
%   Thayer School of Engineering
%   Dartmouth College
%   July 2009
%ooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooo

%load .nod file
nod = load([name,'.nod']);
nn = size(nod,1);

%load complex displacement file:
UC = load([name,'.v3c']);
PC = load('pressure.s1c');

disp_output = [name,'.v3c.ci'];
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

function [arrayout]=dilatearray3D(arrayin,mask)
% Dilate the edge values in an array to avoid errors at the boundary during
% interpolation

s=size(arrayin);

arrayout=arrayin;

mask=double(mask);
if(max(mask(:))~=1)
    error('Dilate3d: Max(mask) ~=1')
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





