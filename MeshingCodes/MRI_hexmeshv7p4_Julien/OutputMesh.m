function [namefile]=OutputMesh(outstm,mask,node,deplacement,mridim,coord,regs)
% Output everything to files

%mask file
namefile.maskoutf=[outstm '.mask.mat'];
save(namefile.maskoutf,'mask');

%nod files
ind=(1:node.nn)';
nodoutf=[outstm '.nod'];
namefile.nodhmgoutf=[outstm '.hom.nod'];
nodregoutf=[outstm '.reg.nod'];


fid=fopen(nodoutf,'w');
fprintf(fid,'%7i %15.8e %15.8e %15.8e %12.5e \n',[ind node.nod ones(node.nn,1)]');
fclose(fid);

fid=fopen(namefile.nodhmgoutf,'w');
fprintf(fid,'%7i %15.8e %15.8e %15.8e  1 \n',[ind node.nod]');
fclose(fid);

fid=fopen(nodregoutf,'w');
fprintf(fid,'%7i %15.8e %15.8e %15.8e %3i \n',[ind node.nod node.regnod]');
fclose(fid);


%elm files
elmind=(1:node.nel)';
elmf=[outstm '.elm'];
fid=fopen(elmf,'w');
fprintf(fid,'%7i %8i %8i %8i %8i %8i %8i %8i %8i %8i %8i %8i %8i %8i %8i %8i %8i %8i %8i %8i %8i %8i %8i %8i %8i %8i %8i %8i \n',[elmind node.in]');
fclose(fid);

% Element file for MRE-Zone-v6 (extra column of ones to indicate
% homogeneous elementally defined properties
namefile.elmhmgf=[outstm '.hom.elm'];
fid=fopen(namefile.elmhmgf,'w');
fprintf(fid,'%7i %8i %8i %8i %8i %8i %8i %8i %8i %8i %8i %8i %8i %8i %8i %8i %8i %8i %8i %8i %8i %8i %8i %8i %8i %8i %8i %8i  1 \n',[elmind node.in]');
fclose(fid);


%index file -> Note this is interpolated index. Need to back-interpolate
%using xin and xout to get values at original MR voxels.
idxf=[outstm '.idx'];
fid=fopen(idxf,'w');
fprintf(fid,'%7i %6i %6i %6i\n',[ind node.idx]');
fclose(fid);

%disp file
namefile.dspoutf=[outstm '.dsp'];
fid=fopen(namefile.dspoutf,'w');
fprintf(fid,'%7i %15.8e %15.8e %15.8e %15.8e %15.8e %15.8e \n',...
    [ind deplacement.uvwr(:,1) deplacement.uvwi(:,1) deplacement.uvwr(:,2) deplacement.uvwi(:,2) deplacement.uvwr(:,3) deplacement.uvwi(:,3)]');
fclose(fid);

%mri file (old magim.mtrin)
magoutf=[outstm '.mri'];
fid=fopen(magoutf,'w');
fprintf(fid,'%7i %15.8e\n',[ind node.magimnod]');
fclose(fid);


%bnod file
namefile.bcoutf=[outstm '.bnd'];
fid=fopen(namefile.bcoutf,'w');
fprintf(fid,'%7i  %7i\n',[(1:node.nbnod)' node.bnod']');
fclose(fid);

% region file for soft prior -> Output the region file at the original data
% resolution so that the regions can be interpolated to the material
% property meshes inside the recon code.
namefile.regoutf=[outstm '.regstack'];
% Format:
% StackSize
% x coords
% y coords
% z coords
% slice1
% slice2
% ..
% Slice_end
fid=fopen(namefile.regoutf,'w');
fprintf(fid,'Array Dimensions \n');
fprintf(fid,'%7i %7i %7i \n',mridim);
fprintf(fid,'x coordinates \n');
fclose(fid);
dlmwrite(namefile.regoutf,coord.xin,'-append','delimiter',' ');
fid=fopen(namefile.regoutf,'a'); fprintf(fid,'y coordniates \n'); fclose(fid);
dlmwrite(namefile.regoutf,coord.yin,'-append','delimiter',' ');
fid=fopen(namefile.regoutf,'a'); fprintf(fid,'z coordniates \n'); fclose(fid);
dlmwrite(namefile.regoutf,coord.zin,'-append','delimiter',' ');
for ii=1:mridim(3)
    dlmwrite(namefile.regoutf,regs(:,:,ii),'-append','delimiter',' ');
end


end