function [eb pb] = CreateSurfaceFrom3DImage(stackfn,stackinfofn)
% Read in the 3D images defined in 'stackfn' and create a surface mesh
% 'stackfn' is a .mat file that contains a matlab structure variable in
% which a 3D matrix holds the segemented MR images.
% 'stackinfofn is a .mat file that contains a matlab structure variable in
% which a 'stackInfo' field extists representing dicom header information.
% The output is the list of elements and nodes of the surface
% Currently only works on segmented images with one label.
%
% Written by:
%   Hamid Ghadyani, Nov 2010


fprintf(' Creating a surface mesh from segmented MR images.\n');

% Conver the 3D image to .inr file format
inrsavefn = [getuserdir filesep 'iso2surf.inr'];
offsavefn = [getuserdir filesep 'iso2surf.off'];
stack_param = stack2inr(stackfn,stackinfofn,inrsavefn);

% Compute parameters for surface mesher
% We assumed that segmented data has only one label
isovalues = unique(stack_param.labels);
if ~ismember(0,isovalues)
    error(' Your segmented images should use 0 value to specify exterior pixels!');
else
    isovalue = setdiff(isovalues,0);
    if length(isovalue)~=1
        error(' Multiple material segmentation is not supported yet!');
    end
end

if stack_param.nRow>=200 && stack_param.nCol >=200
    scx = stack_param.nRow / 2.;
    scy = stack_param.nCol / 2.;
    scz = stack_param.nPln / 2.;
else
    scx = stack_param.nRow ;
    scy = stack_param.nCol ;
    scz = stack_param.nPln ;
end


bssr = scx*scy*3.;
prec = 1e-5;

% The minimum angle of each surface triangle
facet_angle = 23;
%facet_angle = 23.;
% To change size of triangles increase/decrease this value
facet_radius = 4; %4
%facet_radius = 7;
% To increase accuracy of smoothness of surface lower this value
% Beware that it can produce too many nodes on the surface
%facet_distance = 7;
facet_distance =4; %4

% Write out the parameter file:
arglist_fn = [getuserdir filesep 's2s.txt'];
fid = OpenFile(arglist_fn,'wt');
fprintf(fid,'%f\n',isovalue);
fprintf(fid,'%f\n%f\n%f\n',scx,scy,scz);
fprintf(fid,'%f\n%.18f\n',bssr,prec);
fprintf(fid,'%.8f\n',facet_angle);
fprintf(fid,'%.8f\n',facet_radius);
fprintf(fid,'%.8f\n',facet_distance);
fclose(fid);

% Call stack2surface_cgal
delete(offsavefn)
command='stack2surface_cgal';
syscommand = GetSystemCommand(command);
if ~ispc
    eval(['! chmod u+x "' syscommand '"']);
    % I had to link 32bit Mac version of CGAL with libgmp.dylib because of
    % a bug in linking with static library. So here we setup the
    % DYLD_LIBRARY_PATH variable to make sure we can find libgmp.dylib
    dyldpath=[];
    if strcmp(computer,'MACI')
        dircommand = fileparts(syscommand);
        dyldpath = ['DYLD_LIBRARY_PATH=' dircommand ' '];
    end
else
    dyldpath = [];
end
% Unset the LD path so that our standalone executables load the system
% libstdc++ rather than Matlab's.
if isunix && ~ismac
   dyldpath = 'unset LD_LIBRARY_PATH; ';
end

disp('Before Eval'); 
eval(['! ' dyldpath '"' syscommand '" ' inrsavefn ' ' offsavefn ' ' arglist_fn])
disp('After Eval');

% Read the surface
[eb pb] = readOFF(offsavefn);
% pb = floor(pb/prec)*prec;
% Postprocess the surface
nodes=unique(eb(:));
[tf eb] = ismember(eb,nodes);
pb = pb(nodes,:);
eb = FixPatchOrientation(pb,eb,[],1);

input_args.verbose=0;
input_args.type=1;
[junk,junk,junk,myst] = CheckMesh3D(eb,pb,[],input_args);
if isfield(myst,'b') && myst.b~=0 && myst.b~=4
    error(' The brain surface is not closed, single material or manifold!\n');
end
% writenodelm_surface_medit('iso2surf.mesh',eb,pb,[],0)
delete(inrsavefn,offsavefn,arglist_fn);

