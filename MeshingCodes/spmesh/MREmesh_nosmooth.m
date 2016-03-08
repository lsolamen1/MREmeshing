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
%  
%   New addition: Optional input 'scalars'
%   1x2 vector [shear lambda]
%   If not supplied, uses default value
%   
%ooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooo
function MREmesh(fname,scalars,bcf,HC)



close all;
fid=fopen(fname,'rt');

fprintf('-------------------------------------------------------------\n');
fprintf('\n');
fprintf('\tTime-Harmonic MR Elastographic Imaging\n\n');
fprintf('-------------------------------------------------------------\n\n');

%fprintf('Enter path to MRE matlab script files:\n');
%read path information
%dir = fgetl(fid);
%disp(['dir --> ',dir]);
%addpath(dir);

fprintf('\nChoose from one of the following options:\n');
fprintf(' [1]  mesh generation only\n');
fprintf(' [2]  mesh generation and data interpolation\n');
%read option number
todo = fgetl(fid);
disp(['Selection --> ',todo]);
todo = str2num(todo);

%PERFORM MESH GENERATION
if todo == 1 || todo == 2
    fprintf('\nEnter the series name: \n');
    %read series name
    name = fgetl(fid);
    disp(['Series Name --> ',name]);
    
    if(nargin==1)
        scalars=[3000 50000];
        bcf=[name,'.inv.bcs\n'];
        HC=1e-7;
    elseif(nargin==2)
        bcf=[fname,'.inv.bcs\n'];
        HC=1e-7;
    elseif(nargin==3)
        HC=1e-7;
    end
    
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
    DirIndex(1:3,:)=round(DirIndex(1:3,:));
    %retain pixel spacing
    stackInfo.PixelSpacing = DirIndex(4,1:2)'*1e-3;
    stackInfo.SpacingBetweenSlices = DirIndex(4,3)'*1e-3;
    fprintf('\nGENERATING SURFACE MESH\n');
    disp('MRE_surfacemesh.m')
    MRE_surfacemesh(mask,stackInfo,name);
    
    %write MRE frequency to file
    fid1=fopen('freq.tmp','wt');
    fprintf(fid1,'%.4f',freqHz);
    fclose(fid1);
    
    %read mesh resolution
    fprintf('Enter Mesh Resolution: \n');
    res = fgetl(fid);
    disp(['Mesh Resolution --> ',res]);
    
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
    
    copyfile([name,'.elm'],[name,'.elm.header']);
    copyfile([name,'.nod'],[name,'.nod.header']);
    
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
    copyfile([name,'.bel'],[name,'.bel.header']);
    
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
    
    if (todo == 2)
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
        
        %ELASTIC RECONSTRCUTION (original)
        fprintf('\nGENERATING 3D_LAME_FILES.DAT\n');
        dir = 'INV';
        gen_MRE_reconfile(name,dir,itmax,freq,reg,spfilt,add,scalars)
        
        %ELASTIC RECONSTRUCTION (PARDISO)
        fprintf('\nGENERATING 3D_MRE_FILES.DAT\n');
        dir = 'INVpar';
        udopt = 2;
        znopt = 1;
        diag = 0.049;
        ovrlp = 0.10;
        gen_MREpar_reconfile(name,dir,itmax,freq,reg,spfilt,udopt,znopt,diag,ovrlp,scalars);
        
        fprintf('\nPERFORMING DATA INTERPOLATION\n');
        %interpolate data onto finite element mesh using MATLAB scripts
        
        nod = load([name,'.nod']);
        %interpolate measured displacements onto finite element mesh
        gen3ddisp(nod,A,P,MagIm,DirIndex,name);
        
        %strip file headers
        eval(sprintf(['! sed ''1d'' < ',name,'.nod.new.header > ',name,'.nod.new']));
        
        %generate file containing initial guess
        eval(sprintf(['! awk ''{print $1,1,1,' num2str(HC) '}'' < ',name,'.nod > ',name,'.inv.mtr']));
        
        %generate discovery cluster submit-files
        gen_runfiles(add,nprocs,name);
        
        % LAUNCH RECONSTRUCTION
        %fprintf('\nLAUNCHING PARALLEL RECONSTRUCTION\n');
        %! qsub MRE-submit
        %! qsub MREpar-submit
        
        % POROELASTIC RECONSTRUCTION
        fprintf('\nGENERATING 3D_MRPE_FILES.DAT\n');
        dir = 'INV';
        udopt = 1;
        znopt = 1;
        diag = 0.049;
        ovrlp = 0.10;
        gen_MRPE_reconfile(name,dir,itmax,freq,reg,spfilt,udopt,znopt,diag,ovrlp,add,scalars,bcf);
        
        %generate pressure file
        eval(sprintf(['! awk ''{print $1,0,0}'' < ',name,'.nod > pressure.s1c']));
        
        %generate complex input files
        cmplx_input(name);
        
        ! mkdir -p MRPE
        ! mkdir -p MRPE/INV
        ! mv *.ci pressure.s1c 3D_MRPE_Files.dat MRPE-submit MRPE
        
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
        fprintf(fid1,['2   INV/MRE-',name,'-inv.',itmax,'\n']);
        fprintf(fid1,['2   INVpar/MRE-',name,'-inv.',itmax,'\n']);
        fprintf(fid1,['2   MRPE/INV/MRPE-',name,'-inv.',itmax,'\n']);
        fclose(fid1);
        
    end % END IF TODO == 2
    %create auxillary file directories if they do not already exist
    ! mkdir -p INP
    ! mkdir -p LOG
    ! mkdir -p TMP
    ! mkdir -p DIA
    
    %move all auxillary files to their respective directories
    ! mv *.inp INP
    ! mv INP/vol_interp.inp .
    ! mv *.log LOG
    ! mv *.tmp TMP
    ! mv *.dia DIA
    
    %remove files
    %! rm -r -f INP
    ! rm -r -f LOG
    ! rm -r -f TMP
    ! rm -r -f DIA
    ! rm *.out
    eval(sprintf(['! rm ',name,'.vtk']));
    eval(sprintf(['! rm ',name,'.nml']));
    
end % END IF TODO == 1 || TODO == 2
fprintf('\n\n\tTHAT''S ALL FOLKS\n');
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
%clear mask;

% SMOOTH THE DATA
Ds = smooth3(A);

% MESH
% an isovalue of 0.67 creates a surface mesh within 0.01 distance to the
% first and last images
% an isovalue of 0.68 creates a smaller surface mesh within 0.04 distance
% to the first and last images

% June 10 MODIFIED HERE!!!
isoV = 0.6; %0.68;
figure;
hiso = patch(isosurface(Ds,isoV),'FaceColor',[0,0.5,0]);%,'EdgeColor','none');
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

function gen_MRE_reconfile(name,dir,itmax,freq,reg,spfilt,add,scalars)
%ooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooo
%   GEN_MRE_RECONFILE.M
%       This routine is used to generate the input file for the original
%       linearly elastic reconstruction
%
%   INPUTS
%       name        series name
%       dir         directory for reconstruction files
%       itmax       maximum number of iterations
%       freq        MRE frequency
%       reg         Marquardt regularization parameter
%       spfilt      spatial filtering level
%       add         email address
%ooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooo

fid = fopen('3D_Lame_Files.dat','wt');
fprintf(fid,'3D Subzone Inversion Run File\n');
fprintf(fid,'Number of Iterations:\n');
fprintf(fid,[itmax,'\n']);
fprintf(fid,'Scaling Factors: <G> <L>\n');
fprintf(fid,[num2str(round(scalars(1))) '.d0, ' num2str(round(scalars(2))) '.d0' '\n']);  % Defaults set at top of code
fprintf(fid,'Lower Constraints: <G> <L>\n');
fprintf(fid,'1.d0,1.d0\n');
fprintf(fid,'Frequency (Hz):\n');
fprintf(fid,[num2str(freq,'%.4f'),'\n']);
fprintf(fid,'Density of Solid:\n');
fprintf(fid,'1020.d0\n');
fprintf(fid,'Zone generation (1:random zonegen; 0:new zonegen)\n');
fprintf(fid,'1\n');
fprintf(fid,'overlap ratio for (option 0 only)\n');
fprintf(fid,'0.2\n');
fprintf(fid,'Subzone size (for option 0)\n');
fprintf(fid,'0.0155\n');
fprintf(fid,'Marqauardt Regularization Parameters: <lambda>,<delta>\n');
fprintf(fid,[reg,',',reg,',',reg,',1.d-1\n']);
fprintf(fid,'Spatial Filtering Weight: <theta>\n');
fprintf(fid,[spfilt,',',spfilt,',',spfilt,'\n']);
fprintf(fid,'Target Number of Zone Nodes\n');
fprintf(fid,'450\n');
fprintf(fid,'Max number of zone nodes & Max bandwidth\n');
fprintf(fid,'4200 8000 3000\n');
fprintf(fid,'Max Number of Zones\n');
fprintf(fid,'12000\n');
fprintf(fid,'Max Number of connectivity\n');
fprintf(fid,'100\n');
fprintf(fid,'Third party for zonegen flag & zone level Recon information\n');
fprintf(fid,'0 0\n');
fprintf(fid,'Notification Address:\n');
fprintf(fid,[add,'\n']);
fprintf(fid,'Update File Name:\n');
fprintf(fid,[name,'.upd\n']);
fprintf(fid,'Data Style: 0:<mag><phase>1:<real><imag>\n');
fprintf(fid,'0\n');
fprintf(fid,'File List:\n');
fprintf(fid,[name,'.nod.new.header\n']);
fprintf(fid,[name,'.elm.header\n']);
fprintf(fid,[name,'.bnod.header\n']);
fprintf(fid,[name,'.v3r\n']);
fprintf(fid,[name,'.inv.mtr\n']);
fprintf(fid,[dir,'/MRE-',name,'-inv\n']);
fclose(fid);
end

function gen_MREpar_reconfile(name,dir,itmax,freq,reg,spfilt,udopt,znopt,diag,ovrlp,scalars)
%ooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooo
%   GEN_MREPAR_RECONFILE.M
%       This routine is used to generate the input file for the PARDISO-
%       based linearly elastic reconstruction
%   INPUTS
%       name        series name
%       dir         directory for reconstruction files
%       itmax       maximum number of iterations
%       freq        MRE frequency
%       reg         Marquardt regularization parameter
%       spfilt      spatial filtering level
%       udopt       update option
%       znopt       zone option
%       diag        cube diagonal
%       ovrlp       percent zone overlap
%ooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooo

fid = fopen('3D_MRE_Files.dat','wt');
fprintf(fid,'3D Subzone Inversion Run File [PARDISO]\n');
fprintf(fid,'Number of OpenMP Processes: <numprocs>\n');
fprintf(fid,'8\n');
fprintf(fid,'Number of Iterations:\n');
fprintf(fid,[itmax,'\n']);
fprintf(fid,'Global Error Estimate: <1> on <2> off\n');
fprintf(fid,'1\n');
fprintf(fid,'Scaling Factors: <G> <L>\n');
fprintf(fid,[num2str(round(scalars(1))) '.d0, ' num2str(round(scalars(2))) '.d0' '\n']);  % Defaults set at top of code
fprintf(fid,'Lower Constraints: <G> <L>\n');
fprintf(fid,'1.d0,1.d0\n');
fprintf(fid,'Frequency (Hz):\n');
fprintf(fid,[num2str(freq,'%.4f'),'\n']);
fprintf(fid,'Density of Solid:\n');
fprintf(fid,'1020.d0\n');
fprintf(fid,'Joachimowicz Regularization Parameter: <alpha>\n');
fprintf(fid,'5.d0\n');
fprintf(fid,'Marqauardt Regularization Parameters: <lambda>,<delta>\n');
fprintf(fid,[reg,',1.d-1\n']);
fprintf(fid,'Spatial Filtering Weight: <theta>\n');
fprintf(fid,[spfilt,'\n']);
fprintf(fid,'Parameter Update Option: <1> global <2> continuous\n');
fprintf(fid,'%d\n',udopt);
fprintf(fid,'Zone Option: <1> cubes <2> rectangular prisms\n');
fprintf(fid,'%d\n',znopt);
fprintf(fid,'Cube Diagonal:\n');
fprintf(fid,'%.4f\n',diag);
fprintf(fid,'Zone Edge Factor: <znedgfac(1:3)>\n');
fprintf(fid,'1 1 1\n');
fprintf(fid,'Zone Overlap Factor: <znovrlp(1:3)>\n');
fprintf(fid,'%.2f %.2f %.2f\n',ovrlp,ovrlp,ovrlp);
fprintf(fid,'Update File Name:\n');
fprintf(fid,'MRE.upd\n');
fprintf(fid,'File List:\n');
fprintf(fid,[name,'.nod\n']);
fprintf(fid,[name,'.elm\n']);
fprintf(fid,[name,'.bel\n']);
fprintf(fid,[name,'.bnod\n']);
fprintf(fid,[name,'.v3r\n']);
fprintf(fid,[name,'.inv.mtr\n']);
fprintf(fid,[dir,'/MRE-',name,'-inv-stab\n']);
fprintf(fid,[dir,'/relerr.dat\n']);
fclose(fid);
end

function gen_MRPE_reconfile(name,dir,itmax,freq,reg,spfilt,udopt,znopt,diag,ovrlp,add,scalars,bcf)
%ooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooo
%   GEN_MRPE_RECONFILE.M
%     This routine is used to generate the input file for the poroelastic 
%     reconstruction
%
%   INPUTS
%     name        series name
%     dir         directory for reconstruction files
%     itmax       maximum number of iterations
%     freq        MRE frequency
%     reg         Marquardt regularization parameter
%     spfilt      spatial filtering level
%     udopt       update option
%     znopt       zone option
%     diag        cube diagonal
%     ovrlp       percent zone overlap
%ooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooo

fid = fopen('3D_MRPE_Files.dat','wt');
fprintf(fid,'3D Subzone Inversion Run File [PARDISO]\n');
fprintf(fid,'Number of OpenMP Processes: <numprocs>\n');
fprintf(fid,'8\n');
fprintf(fid,'Number of Iterations:\n');
fprintf(fid,[itmax,'\n']);
fprintf(fid,'Scaling Factors: <G> <L>\n');
fprintf(fid,[num2str(round(scalars(1))) '.d0, ' num2str(round(scalars(2))) '.d0' '\n']);  % Defaults set at top of code
fprintf(fid,'Lower Constraints: <G> <L>\n');
fprintf(fid,'1.d0,1.d0\n');
fprintf(fid,'Upper Constraints: <G> <L>\n');
fprintf(fid,'1.d5,1.d7\n');
fprintf(fid,'Frequency (Hz):\n');
fprintf(fid,[num2str(freq,'%.4f'),'\n']);
fprintf(fid,'porosity,density of solid,fluid,apparent:\n');
fprintf(fid,'0.2d0,1020.d0,1000.d0,150.d0\n');
fprintf(fid,'Joachimowicz Regularization Parameter: <alpha>\n');
fprintf(fid,'75.d0\n');
fprintf(fid,'Marqauardt Regularization Parameters: <lambda>,<delta>\n');
fprintf(fid,[reg,',1.d-1\n']);
fprintf(fid,'Spatial Filtering Weight: <theta>\n');
fprintf(fid,[spfilt,'\n']);
fprintf(fid,'Parameter Update Option: <1> global <2> continuous\n');
fprintf(fid,'%d\n',udopt);
fprintf(fid,'Zone Option: <1> cubes <2> rectangular prisms\n');
fprintf(fid,'%d\n',znopt);
fprintf(fid,'Cube Diagonal:\n');
fprintf(fid,'%.4f\n',diag);
fprintf(fid,'Zone Edge Factor: <znedgfac(1:3)>\n');
fprintf(fid,'1 1 1\n');
fprintf(fid,'Zone Overlap Factor: <znovrlp(1:3)>\n');
fprintf(fid,'%.2f %.2f %.2f\n',ovrlp,ovrlp,ovrlp);
fprintf(fid,'Update File Name:\n');
fprintf(fid,'MRPE.upd\n');
fprintf(fid,'File List:\n');
fprintf(fid,['../',name,'.nod\n']);
fprintf(fid,['../',name,'.elm\n']);
fprintf(fid,['../',name,'.bel\n']);
fprintf(fid,[bcf '\n']);
fprintf(fid,[name,'.v3c.ci\n']);
fprintf(fid,['../',name,'.inv.mtr\n']);
fprintf(fid,'pressure.s1c.ci\n');
fprintf(fid,[dir,'/MRPE-',name,'-inv\n']);
fprintf(fid,[dir,'/displace.v3c\n']);
fprintf(fid,[dir,'/pressure.s1c\n']);
fprintf(fid,[dir,'/relerr.dat\n']);
fclose(fid);
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

%% Create mesh grid to interpolate displacements and display
% z dir includes the 'extra' slices in the nslice dimension
[yy,xx,zz]=meshgrid(dx:dx:(nx)*dx,dy:dy:(ny)*dy,0:dz:(nslice+1)*dz); % 

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
load('Mask.mat');
MagIm = MagIm.*mask;

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

%saveas(gcf,'mesh-validation.jpg','jpeg');%close;

%% Interpolate displacements from grid to the FE mesh
Dxr=interp3(xx,yy,zz,pDr(:,:,:,1),nod(:,2),nod(:,3),nod(:,4),'*cubic');
Dyr=interp3(xx,yy,zz,pDr(:,:,:,2),nod(:,2),nod(:,3),nod(:,4),'*cubic');
Dzr=interp3(xx,yy,zz,pDr(:,:,:,3),nod(:,2),nod(:,3),nod(:,4),'*cubic');

Dxi=interp3(xx,yy,zz,pDi(:,:,:,1),nod(:,2),nod(:,3),nod(:,4),'*cubic');
Dyi=interp3(xx,yy,zz,pDi(:,:,:,2),nod(:,2),nod(:,3),nod(:,4),'*cubic');
Dzi=interp3(xx,yy,zz,pDi(:,:,:,3),nod(:,2),nod(:,3),nod(:,4),'*cubic');

if(exist('Shear_sinkus.mat','file')) % Output tetrahedral initial guess file
    fprintf(1,'Outputting Initial Guess file from Shear_sinkus.mat');
    load Shear_sinkus.mat  
    IGr2=ones(size(IGr,1),size(IGr,2),size(IGr,3)+2).*mean(IGr(IGr>100));
    IGr2(:,:,2:end-1)=IGr;
    clear IGr
    IGi2=ones(size(IGi,1),size(IGi,2),size(IGi,3)+2).*mean(IGi(IGi>10));
    IGi2(:,:,2:end-1)=IGi;
    clear IGi
        
    IGr_nod=interp3(xx,yy,zz,IGr2(:,:,:),nod(:,2),nod(:,3),nod(:,4),'*cubic');
    IGi_nod=interp3(xx,yy,zz,IGi2(:,:,:),nod(:,2),nod(:,3),nod(:,4),'*cubic');
    
    nu=0.47; % poissons ratio
    mu_sc=3000;
    lam_sc=5d4;
    
    fid=fopen([name '_mu' int2str(mu_sc) '_lam' int2str(lam_sc) '.IGr'],'w')
      fprintf(fid,'%6d  %14.6E  %14.6E\n',[nod(:,1) IGr_nod./mu_sc IGr_nod.*(2*nu/(1-2*nu)/lam_sc)]');
    fclose(fid);
    fid=fopen([name '_mu' int2str(mu_sc) '.IGi'],'w')
      fprintf(fid,'%6d  %14.6E',[nod(:,1) IGi_nod]');
    fclose(fid);
else
    fprintf(1,'Shear_sinkus.mat not present: Run Direct inversion first to generate an initial guess file \n');
end
    
    

% real displacements
disp(:,1)=nod(:,1);
disp(:,2)=Dxr;
disp(:,3)=Dyr;
disp(:,4)=Dzr;

% imaginary displacements
cdisp(:,1)=nod(:,1);
cdisp(:,2)=Dxi;
cdisp(:,3)=Dyi;
cdisp(:,4)=Dzi;

%% Transform Nodal Coordinates according with DirIndex directions
nod_transform(nod,disp,cdisp,DirIndex,name);
end

function nod_transform(nod,disp,cdisp,DirIndex,name)
%ooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooo
%   NOD_TRANSFORM.M
%     This routine is used to transform the nodal co-ordinate and
%     displacement data using DirIndex
%
%   INPUTS
%     nod         [nn,4] array containing nodal coordiantes
%     disp        [nn,3] array containing real-valued displacements
%     cdisp       [nn,6] array containing complex-valued displacements
%     DirIndex    4 X 6 matrix containing information regarding how to
%                   transform displacement data and nodal co-ordinates
%     name        series name
%ooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooo

% Determine Column Locations for Transformatrion Matricies
Index_1=find(abs(DirIndex(1,:))==1);
Index_2=find(abs(DirIndex(2,:))==1);
Index_3=find(abs(DirIndex(3,:))==1);

% Get Voxel Spacing
dx = DirIndex(4,1)*1e-3;
dy = DirIndex(4,2)*1e-3;
dz = DirIndex(4,3)*1e-3;

% Nodal Co-ordinate transformation
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
new_disp(:,1) = disp(:,1);
new_disp(:,2) = disp(:,Index_1(2)-3+1)*DirIndex(1,Index_1(2));
new_disp(:,3) = disp(:,Index_2(2)-3+1)*DirIndex(2,Index_2(2));
new_disp(:,4) = disp(:,Index_3(2)-3+1)*DirIndex(3,Index_3(2));

% Imaginary Component
new_cdisp(:,1) = cdisp(:,1);
new_cdisp(:,2) = cdisp(:,Index_1(2)-3+1)*DirIndex(1,Index_1(2));
new_cdisp(:,3) = cdisp(:,Index_2(2)-3+1)*DirIndex(2,Index_2(2));
new_cdisp(:,4) = cdisp(:,Index_3(2)-3+1)*DirIndex(3,Index_3(2));

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

function gen_runfiles(add,nprocs,name)
%ooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooo
%   GEN_RUNFILES.M
%
%   DESCRIPTION:
%     GEN_RUNFILES(add,nprocs,name) generates a PBS jobfile for the 
%     parallelMRE reconstructions for use on the DISCOVERY cluster
%
%   INPUT:
%     add         a string containing the email address of the operator
%     nprocs      an integer representing the total number of cpus 
%                 requested
%     name        series name
%
%   DEFAULTS:
%     walltime - 48 hrs
%     output/error files - combined
%
%   Phillip R. Perriñez, Ph.D.
%   Thayer School of Engineering
%   Dartmouth College
%   July 2009
%ooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooo

%compute number of nodes required
nnodes=(str2num(nprocs)/8);
job=['#PBS -N ',name];
CPU =['#PBS -l nodes=',num2str(nnodes),':','ppn=8'];

%% ELASTIC (original)
fid=fopen('MRE-submit', 'wt');
fprintf(fid,'%s \n','# declare a name for this job to be sample_job');
fprintf(fid,'%s \n',job);
fprintf(fid,'%s \n','# request the queue (enter the possible names, if omitted, serial is the default)');
fprintf(fid,'%s \n','#PBS -q default');
fprintf(fid,'%s \n','# request  node');
fprintf(fid,'%s \n',CPU);
fprintf(fid,'%s \n','# request 16 hours of wall time');
fprintf(fid,'%s \n','#PBS -l walltime=16:00:00');
fprintf(fid,'%s \n','#combine PBS standard output and error files');
fprintf(fid,'%s \n','#PBS -j oe');
fprintf(fid,'%s \n','# mail is sent to you when the job starts and when it terminates or aborts');
fprintf(fid,'%s \n','#PBS -m bea');
fprintf(fid,'%s \n','# specify your email address');
fprintf(fid,'%s \n', ['#PBS -M',' ',add]);
fprintf(fid,'%s \n','# By default, PBS scripts execute in your home directory, not the');
fprintf(fid,'%s \n','# directory from which they were submitted. The following line');
fprintf(fid,'%s \n','# places you in the directory from which the job was submitted.');
fprintf(fid,'%s \n','cd $PBS_O_WORKDIR');
fprintf(fid,'%s \n','# run the program');
fprintf(fid,'%s \n','module load mpich/1.2.7p1-pgi7.1.3');
fprintf(fid,'%s \n','mpiexec -comm p4 ~/bin/3D_recon_MP_newton.x');
fprintf(fid,'%s \n','mpiexec -pernode -comm none /opt/mpich/mpich1.2.7p1-pgi7.1.3/sbin/cleanipcs');
fprintf(fid,'%s \n','exit 0');
fclose(fid);

%% ELASTIC (PARDISO)
CPU =['#PBS -l nodes=',num2str(nnodes),':','ppn=8:amd'];
fid=fopen('MREpar-submit', 'wt');
fprintf(fid,'%s \n','# declare a name for this job to be sample_job');
fprintf(fid,'%s \n',job);
fprintf(fid,'%s \n','# request the queue (enter the possible names, if omitted, serial is the default)');
fprintf(fid,'%s \n','#PBS -q default');
fprintf(fid,'%s \n','# request  node');
fprintf(fid,'%s \n',CPU);
fprintf(fid,'%s \n','# request 30 hours of wall time');
fprintf(fid,'%s \n','#PBS -l walltime=30:00:00');
fprintf(fid,'%s \n','#combine PBS standard output and error files');
fprintf(fid,'%s \n','#PBS -j oe');
fprintf(fid,'%s \n','# mail is sent to you when the job starts and when it terminates or aborts');
fprintf(fid,'%s \n','#PBS -m bea');
fprintf(fid,'%s \n','# specify your email address');
fprintf(fid,'%s \n', ['#PBS -M',' ',add]);
fprintf(fid,'%s \n','# By default, PBS scripts execute in your home directory, not the');
fprintf(fid,'%s \n','# directory from which they were submitted. The following line');
fprintf(fid,'%s \n','# places you in the directory from which the job was submitted.');
fprintf(fid,'%s \n','cd $PBS_O_WORKDIR');
fprintf(fid,'%s \n','# run the program');
fprintf(fid,'%s \n','mpirun -machinefile $PBS_NODEFILE -np 8 /home/mmcgarry/code/linearcode/pardiso/stabilized/MRE-pardiso.x');
fprintf(fid,'%s \n','mpiexec -pernode -comm=none /opt/mpich/mpich1.2.7p1-intel11.0/sbin/cleanipcs');
fprintf(fid,'%s \n','exit 0');
fclose(fid);

fid=fopen('MPICH2par-submit', 'wt');
fprintf(fid,'%s \n','#!/bin/bash -l ');
fprintf(fid,'%s \n','# declare a name for this job ');
fprintf(fid,'%s \n',job);
fprintf(fid,'%s \n','# request the queue (enter the possible names, if omitted, serial is the default)  ');
fprintf(fid,'%s \n','#PBS -q default  ');
fprintf(fid,'%s \n','# request  node  ');
fprintf(fid,'%s \n','#PBS -l nodes=1:ppn=8 ');
fprintf(fid,'%s \n','#PBS -l feature=amd ');
fprintf(fid,'%s \n','# request some hours of wall time  ');
fprintf(fid,'%s \n','#PBS -l walltime=36:00:00  ');
fprintf(fid,'%s \n','#combine PBS standard output and error files  ');
fprintf(fid,'%s \n','#PBS -j oe  ');
fprintf(fid,'%s \n','# mail is sent to you when the job starts and when it terminates or aborts  ');
fprintf(fid,'%s \n','##PBS -m bea  ');
fprintf(fid,'%s \n','# specify your email address  ');
fprintf(fid,'%s \n','##PBS -M matthew.d.mcgarry@dartmouth.edu  ');
fprintf(fid,'%s \n','# By default, PBS scripts execute in your home directory, not the  ');
fprintf(fid,'%s \n','# directory from which they were submitted. The following line  ');
fprintf(fid,'%s \n','# places you in the directory from which the job was submitted.  ');
fprintf(fid,'%s \n','cd $PBS_O_WORKDIR  ');
fprintf(fid,'%s \n','# run the program  ');
fprintf(fid,'%s \n','cat $PBS_NODEFILE | uniq > node_file  ');
fprintf(fid,'%s \n','nnodes=$(cat node_file | wc -l)  ');
fprintf(fid,'%s \n','nprocs=$(cat $PBS_NODEFILE | wc -l)  ');
fprintf(fid,'%s \n','mpdboot -n $nnodes -f node_file -r rsh  ');
fprintf(fid,'%s \n','export MKL_NUM_THREADS=1  ');
fprintf(fid,'%s \n','echo Nodes $nnodes  ');
fprintf(fid,'%s \n','echo Procs $nprocs  ');
fprintf(fid,'%s \n','cat node_file  ');
fprintf(fid,'%s \n','mpiexec -launcher rsh -n $nprocs /home/mmcgarry/code/linearcode/pardiso/stabilized/MRE-pardiso.x');
fprintf(fid,'%s \n','mpdallexit ');
fprintf(fid,'%s \n','exit 0 ');


%% POROELASTIC (PARDISO)
fid=fopen('MRPE-submit', 'wt');
fprintf(fid,'%s \n','# declare a name for this job to be sample_job');
fprintf(fid,'%s \n',job);
fprintf(fid,'%s \n','# request the queue (enter the possible names, if omitted, serial is the default)');
fprintf(fid,'%s \n','#PBS -q default');
fprintf(fid,'%s \n','# request  node');
fprintf(fid,'%s \n',CPU);
fprintf(fid,'%s \n','# request 47 hours of wall time');
fprintf(fid,'%s \n','#PBS -l walltime=47:00:00');
fprintf(fid,'%s \n','#combine PBS standard output and error files');
fprintf(fid,'%s \n','#PBS -j oe');
fprintf(fid,'%s \n','# mail is sent to you when the job starts and when it terminates or aborts');
fprintf(fid,'%s \n','#PBS -m bea');
fprintf(fid,'%s \n','# specify your email address');
fprintf(fid,'%s \n', ['#PBS -M',' ',add]);
fprintf(fid,'%s \n','# By default, PBS scripts execute in your home directory, not the');
fprintf(fid,'%s \n','# directory from which they were submitted. The following line');
fprintf(fid,'%s \n','# places you in the directory from which the job was submitted.');
fprintf(fid,'%s \n','cd $PBS_O_WORKDIR');
fprintf(fid,'%s \n','# run the program');
fprintf(fid,'%s \n','mpirun -machinefile $PBS_NODEFILE -np 8 /home/pperrine/MRPE-PARDISO/MRPE-pardiso-v2.2-ethernet.x');
fprintf(fid,'%s \n','mpiexec -pernode -comm=none /opt/mpich/mpich1.2.7p1-intel11.0/sbin/cleanipcs');
fprintf(fid,'%s \n','exit 0');
fclose(fid);

fid=fopen('MRPE-submit', 'wt');
fprintf(fid,'%s \n','#!/bin/bash -l');
fprintf(fid,'%s \n','# declare a name for this job to be sample_job ');
fprintf(fid,'%s \n',job);
fprintf(fid,'%s \n','# request the queue (enter the possible names, if omitted, serial is the default) ');
fprintf(fid,'%s \n','#PBS -q default ');
fprintf(fid,'%s \n','# request  node ');
fprintf(fid,'%s \n','#PBS -l nodes=1:ppn=8');
fprintf(fid,'%s \n','#PBS -l feature=amd');
fprintf(fid,'%s \n','# request 24 hours of wall time ');
fprintf(fid,'%s \n','#PBS -l walltime=24:00:00 ');
fprintf(fid,'%s \n','#combine PBS standard output and error files ');
fprintf(fid,'%s \n','#PBS -j oe ');
fprintf(fid,'%s \n','# mail is sent to you when the job starts and when it terminates or aborts ');
fprintf(fid,'%s \n','#PBS -m bea ');
fprintf(fid,'%s \n','# specify your email address ');
fprintf(fid,'%s \n','#PBS -M matthew.d.mcgarry@dartmouth.edu ');
fprintf(fid,'%s \n','# By default, PBS scripts execute in your home directory, not the ');
fprintf(fid,'%s \n','# directory from which they were submitted. The following line ');
fprintf(fid,'%s \n','# places you in the directory from which the job was submitted. ');
fprintf(fid,'%s \n','cd $PBS_O_WORKDIR ');
fprintf(fid,'%s \n','# run the program ');
fprintf(fid,'%s \n','cat $PBS_NODEFILE | uniq > node_file');
fprintf(fid,'%s \n','nnodes=$(cat node_file | wc -l)');
fprintf(fid,'%s \n','nprocs=$(cat $PBS_NODEFILE | wc -l)');
fprintf(fid,'%s \n','mpdboot -n $nnodes -f node_file -r rsh');
fprintf(fid,'%s \n','mpiexec -n $nprocs /home/apattiso/src/MRPE_INV/original/newversion/MRPE-pardiso.x');
fprintf(fid,'%s \n','mpdallexit');
fprintf(fid,'%s \n','exit 0 ');
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
