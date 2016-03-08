function [e p] = ResectBrainWithStereovision(brainmeshfn, stereovisionfn, outputfn)
% This routine is used to remove part of a tetrahedral mesh (brainmesh)
% using a stereovision surface.
% It will first compute the surface of the tetrahedral mesh and then
% perform a boolean operation between that and the stereovision surface
% using CGAL Nef polyhedrons. The resulting surface will be used to remesh
% and produce new set of tetrahedral brain mesh.
% 
% Input:
%  brainmeshfn: filename containing the tetrahedral brain mesh
%  stereovisionfn: filename containing the stereovision surface mesh
%  outputfn: filename to write the new mesh to ('_resected' will be appened)
% 
% Output:
%  e: list of new tetrahedra in brain mesh
%  p: coordinates of nodes in brain mesh
% 
% Written by:
%  Hamid Ghadyani, Nov 2010

% Read in the mesh files
[es ps] = ReadSurfaceMesh(stereovisionfn);
[eb pb] = ReadTetraMesh(brainmeshfn);

% Extract the surface of the brain mesh:
if size(eb,2) >= 4
    [eb pb] = boundfaces(pb,eb);
    % Check the surface's integrity
    input_args.verbose=0;
    input_args.type=1;
    [junk,junk,junk,myst] = CheckMesh3D(eb,pb,[],input_args);
    if isfield(myst,'b') && myst.b~=0 && myst.b~=4
        error(' The brain surface is not closed, single material or manifold!\n');
    end
end
eb = FixPatchOrientation(pb,eb,[],1);

axis = textread([fileparts(stereovisionfn) filesep 'microscope_axis.txt']);
% Make a closed surface out of the stereovision surface. This is necessary
% in order for CGAL routines to work.
[es ps] = CreateClosedSurface(es,ps,axis);

% Perform the resection
[eb pb] = ResectBrain_SurfaceMesh(eb,pb,es,ps);

% These parameters control the shape and size of the new tetrahedral brain
% mesh.
% To increase number of tetrahedrons, use a smaller cell_size.
% To increase or decrease the number of nodes on the surface of the brain,
% adjust facet_size and use values of less thatn 0.4 for facet_distance

% Controls the minimum angle of triangles on the surface of the brain mesh
% If facet_angle>30, there would be no guarantee that mesh generator will
% terminate
criteria.facet_angle = 23.0; % degrees

% Controls the size of the triangles on the surface of the brain mesh.
criteria.facet_size = 2.5;

% Controls the accuracy of the representing the surface mesh. The smaller
% this number, the more nodes will appear on the surface mesh
criteria.facet_distance = 0.7;

% Controls the shape quality of tetrahedral elements, the defualt value
% usually creates good elements
criteria.cell_radisu_edge = 3;

% Controls the size of tetrahedral elements, This value should not be
% _very_ different than 'facet_size'.
% It also affects how many node numbers will be created within the
% tetrahedral mesh
criteria.cell_size = 3.5;

% Call the mesher
[e p] = MakeTetraMesh(eb, pb, criteria);
% Discard material information for now
e = e(:,1:4);

% Write out the results
bel = getBdyFromMesh (e);
if nargin == 2
    foo = remove_extension(brainmeshfn);
    fn = [foo '_resected'];
elseif nargin >= 3
    foo = remove_extension(outputfn);
    fn = [foo '_resected'];
end
WriteOutputFiles(fn,e,p,bel)

fprintf(' Done.\n');


