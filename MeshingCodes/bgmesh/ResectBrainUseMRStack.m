function [e p] = ResectBrainUseMRStack(stackfn, stackinfofn, ...
                         stereovisionfn, resection_normal, ...
                         outputfn, meshing_criteria)
% [e p] = ResectBrainUseMRStack(stackfn, stackinfofn, ...
%                         stereovisionfn, resection_normal, ...
%                         outputfn, meshing_criteria)
% 
% This routine reads the segmented images of a brain (stored in stackfn and
% stackinfofn), creates a surface mesh and then removes part of the
% surface mesh using the 'stereovisionfn'. The resulting surface mesh is
% then used to create a tetrahedral mesh (saved in 'outputfn')
% 
% It will first compute the surface of the tetrahedral mesh and then
% perform a boolean operation between that and the stereovision surface
% using CGAL Nef polyhedrons. 
% 
% Input:
% 
%  stackfn: Name of a .mat file containing segmented images of the brain
%  stackinfofn: Name of a .mat file containing 'stackInfo' (dicom header
%               information)
%  stereovisionfn: filename containing the stereovision surface mesh
%  resection_normal: contains information about the normal (1,:) to brain
%                    surface before surgery and center of craniotomy (2,:)
%  outputfn: filename to write the new mesh to ('_resected' will be appened)
%  meshing_criteria: A structure holding various tetrahedral meshing
%                    See documentation in SetupMeshingCriteria function.
% 
% 
% Output:
%  e: list of new tetrahedra in brain mesh
%  p: coordinates of nodes in brain mesh
% 
% 
% Written by:
%  Hamid Ghadyani, Nov 2010
% Modified by:
%  Hamid Ghadyani, Oct 2011
%       Added mesh refinement options around the tumor region.
%       Updated to work with new pieces of input information from operating
%       room.
% 
% Ex: 
% [e p] = ResectBrainUseMRStack('OR2011_03_28_003/Mesh/brain_bw.mat','OR2011_03_28_003/Mesh/stackInfo.mat','OR2011_03_28_003/Cameras/frame003/frame003',[-0.447117051926064	-0.0988338307372974	-0.0307691789653279;55.4879908256881	110.117376146789	106.746926605505],'foo');
% 
if nargin < 6
    meshing_criteria = [];
end
    
warning('off','MATLAB:DELETE:FileNotFound');
tic
% Read in the stereovision mesh files
[es ps] = ReadSurfaceMesh(stereovisionfn);
% microscope_axis = textread([fileparts(stereovisionfn) filesep 'microscope_axis.txt']);

% Check if the stereovison surface is self-intersecting or not:
r = surface_self_intersect_check(es,ps);
if r ~= 0
    writenodelm_surface_medit('bad_stereo.mesh',es,ps);
    fprintf(' stereovision surface seems to be intersecting itself.')
    fprintf(' check your surface to ensure valid final mesh.')
end

% Make a closed surface out of the stereovision surface. This is necessary
% in order for CGAL routines to work.
[es ps] = CreateClosedSurface(es,ps,resection_normal);

% Create a surface mesh out of the segmented brain images:
[eb pb] = CreateSurfaceFrom3DImage(stackfn,stackinfofn);

% % Plot the brain and stereovision surfaces
% trisurf(eb,pb(:,1),pb(:,2),pb(:,3),ones(length(pb),1),'FaceAlpha',1);
% hold on
% trisurf(es,ps(:,1),ps(:,2),ps(:,3),2*ones(size(ps,1),1),'FaceAlpha',0.6);
% axis equal

% Perform the resection
[ebn pbn] = ResectBrain_SurfaceMesh(eb,pb,es,ps);

% re-mesh the resected brain surface
% subdivide_depth specifies how refine the resulting mesh will be. The
% higer it is the more refinement (nodes) will be generated.
subdivide_depth = 1; 
[ebn pbn] = SubdivideSurfaceMesh(ebn,pbn,subdivide_depth);

criteria = SetupMeshingCriteria(meshing_criteria);

% Call the mesher
[e p] = MakeTetraMesh(ebn, pbn, criteria);
% Discard material information for now
e = e(:,1:4);

% Renumber nodes to minimize matrix bandwidth
[e p]=bandwidth_reduction(e,p);

% Write out the results
bel = getBdyFromMesh (e);
if nargin == 3
    foo = remove_extension(stackfn);
    fn = [foo '_resected'];
elseif nargin >= 4
    foo = remove_extension(outputfn);
    fn = [foo '_resected'];
end
WriteOutputFiles(fn,e,p,bel)

fprintf(' Done.\n');

warning('on','MATLAB:DELETE:FileNotFound');
toc
fprintf('\n');







function criteria = SetupMeshingCriteria(meshing_criteria)

fieldnames = {'tumor_cen','surf_cen','ref_ratio','facet_angle',...
    'facet_size','facet_distance','cell_radius_edge','cell_size',...
    'R1','R2','S1','S2'};

if isempty(meshing_criteria)
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
    criteria.facet_size = 3.0;

    % Controls the accuracy of the representing the surface mesh. The smaller
    % this number, the more nodes will appear on the surface mesh
    criteria.facet_distance = 0.6;

    % Controls the shape quality of tetrahedral elements, the defualt value
    % usually creates good elements
    criteria.cell_radius_edge = 3;

    % Controls the size of tetrahedral elements, This value should not be
    % _very_ different than 'facet_size'.
    % It also affects how many node numbers will be created within the
    % tetrahedral mesh
    criteria.cell_size = 3;

    % Specify the center of tumor
    criteria.tumor_cen = [93.7500  117.1875   97.5000];
    
    % Specify the center of craniotomy
    criteria.surf_cen = [55.4880  110.1174  106.7469];
    
    % Following options specify the mesh density around the surf_cen point
    criteria.R1 = 30;  % how far around the surf_cen we need to refine
    criteria.S1 = 1.6; % Size of cell around the surf_cen
    
    % Following options specify the mesh density around the tumor_cen point
    criteria.R2 = 5; % how far around the surf_cen we need to refine
    criteria.S2 = 3; % Size of cell around the tumor_cen

    % Specify the scale of refinement within the above sphere.
    % This number will be multiplied by criteria.cell_size to determine the
    % size of the tetrahedra within sphere. A value of 1.0 will create a
    % unifirm mesh.
    criteria.ref_ratio = 0.4;
else
    for i=1:length(fieldnames)
        if ~isfield(meshing_criteria,fieldnames{i})
            error([' Provided parameters for meshing are not enough. ' ...
            'Check out the documentation in SetupMeshingCriteria ' ...
            'function within ResectBrainUserMRStack.m file']);
        end
    end
    criteria = meshing_criteria;
end
