function [e p] = MakeTetraMesh(eb, pb, criteria)
% Using the polyhedron defined in 'eb' and 'pb' and size/shape settings in
% 'criteria', this routine uses CGAL library to create a tetrahedral mesh.
% 
% criteria is supposed to have following fields:
% facet_angle, facet_size, facet_distance, cell_radius_edge, cell_size
% 
% These parameters control the shape and size of the new tetrahedral brain
% mesh.
% To increase number of tetrahedrons, use a smaller cell_size.
% To increase or decrease the number of nodes on the surface of the brain,
% adjust facet_size and use values of less thatn 0.5 for facet_distance

% facet_angle controls the minimum angle of triangles on the surface of the brain mesh
% If facet_angle>30, there would be no guarantee that mesh generator will
% terminate
% criteria.facet_angle = 23.0; % degrees

% facet_size controls the size of the triangles on the surface of the brain mesh.
% criteria.facet_size = 2.5;

% facet_distance controls the accuracy of the representing the surface mesh. The smaller
% this number, the more nodes will appear on the surface mesh
% criteria.facet_distance = 0.7;

% cell_radisu_edge controls the shape quality of tetrahedral elements, the defualt value
% usually creates good elements
% criteria.cell_radisu_edge = 3;

% cell_size controls the size of tetrahedral elements, This value should not be
% _very_ different than 'facet_size'.
% It also affects how many node numbers will be created within the
% tetrahedral mesh
% criteria.cell_size = 3.5;
% 
% Written by:
%  Hamid Ghadyani, Nov 2010

fprintf(' Making tetrahedral mesh from resected brain surface mesh.\n');

% Default values
if nargin <3
    facet_angle=22.0;
    facet_size = 5;
    facet_distance = 0.8;
    cell_radius_edge = 4.0;
    cell_size = 2.0;
    tumor_cen = [100.0000 100.000 20.00];
    surf_cen = [100.0000 100.0000 25.000];
    R1 = 30; S1= 3; R2 = 5; S2 = 3;
    ref_ratio = 1;
else
    facet_angle = criteria.facet_angle;
    facet_size = criteria.facet_size;
    facet_distance = criteria.facet_distance;
    cell_radius_edge = criteria.cell_radius_edge;
    cell_size = criteria.cell_size;
    tumor_cen = criteria.tumor_cen;
    surf_cen = criteria.surf_cen;
    R1 = criteria.R1; 
    S1 = criteria.S1;
    R2 = criteria.R2; S2 = criteria.S2;
    ref_ratio = criteria.ref_ratio;
end

inputfn = [getuserdir filesep 'inputpolyhedron.off'];
criteriafn = [getuserdir filesep 'meshcriteria.txt'];
outmeshfn = [getuserdir filesep 'outmesh.mesh'];

%delete(inputfn,criteriafn,outmeshfn);
writeOFF(inputfn, eb, pb);

fid = OpenFile(criteriafn,'wt');
fprintf(fid,'%f\n%f\n%f\n%f\n%f\n', ...
        facet_angle, facet_size, facet_distance,...
        cell_radius_edge, cell_size);
fprintf(fid,'%f %f %f\n',tumor_cen);
fprintf(fid,'%f %f %f\n',surf_cen);
fprintf(fid,'%f\n%f\n%f\n%f\n%f\n',R1,R2,S1,S2,ref_ratio);
fprintf(fid,...
   '\n\n\n%s\n%s\n%s\n%s\n%s\n%s\n%s\n%s\n%s\n%s\n%s\n%s\n%s\n%s\n%s\n',...
        'cfs >> facet_angle;',...
        'cfs >> facet_size;',...
        'cfs >> facet_distance;',...
        'cfs >> cell_radius_edge;',...
        'cfs >> cell_size;',...
        'cfs >> xa; // Point A (tumor side)',...
        'cfs >> ya;',...
        'cfs >> za;',...
        'cfs >> xb; // Point B (brain side)',...
        'cfs >> yb;',...
        'cfs >> zb;',...
        'cfs >> r1; // radius on brain side',...
        'cfs >> r2; // radius on tumor side',...
        'cfs >> s1; // tet size on brain side',...
        'cfs >> s2; // tet size on tumor side',...
        'cfs >> ref_ratio;');
fclose(fid);

syscommand = GetSystemCommand('make_tetmesh_refine_sphere');
if ~ispc
    eval(['! chmod u+x "' syscommand '"']);
end
% Unset the LD path so that our standalone executables load the system
% libstdc++ rather than Matlab's.
dyldpath = '';
if isunix && ~ismac
    dyldpath = 'unset LD_LIBRARY_PATH; ';
end

makemeshcommand = ['! "' syscommand '" ' inputfn ' ' outmeshfn ' ' criteriafn];
%makemeshcommand = ['! ' dyldpath '"' syscommand '" ' inputfn ' ' outmeshfn ' ' criteriafn];

eval(makemeshcommand);
%delete(criteriafn,inputfn);

[e p] = readMEDIT(outmeshfn);

end




% shfn);



