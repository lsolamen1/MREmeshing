function [e p] = MakeTetraMeshFromTaggedStack(stack, criteria)
% 
% Written by:
%  Hamid Ghadyani, Nov 2010

fprintf(' Making tetrahedral mesh from image stack.\n');

if nargin < 2
    basedir = '/Users/hamid_r_ghadyani/Dropbox/Research/wip/xiao/release/brain_tagged';
    cd(basedir)
    load('stack')
    load('stackInfo')
    
    
    criteria = stackInfo;
    
    criteria.facet_angle = 22.0;
    criteria.facet_size = 3.0;
    criteria.facet_distance = 0.8;
    criteria.cell_radius_edge = 4.0;
    criteria.cell_size = 5.0;
    criteria.delmedit=0;
    
    criteria.special_subdomain_size = 1.;
    criteria.special_subdomain_label = 2;
end

stack = uint8(stack);
[e p]=RunCGALMeshGenerator(stack,criteria);



