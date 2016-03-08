function bel = getBdyFromMesh (t)
%
% function bel = getBdyFromMesh (t)
%
% This program extracts the boundary elements from a given connectivity of
% tetrahedral elements (N-by-4), which is the output of spmesh, or 2D
% triangular elements (N-by-3).
% Returns bel (M-by-3) (given 3D mesh) or bel (M-by-2) (given 2D mesh).
% The normal of each bel should be pointing outward.
%
% Based on function [bdye bdyn] = extract_boundary
% by Hamid Ghadyani (hamid.r.ghadyani@dartmouth.edu), Apr 1st, 2009
%
% Removed UI's, no need for input/output files, removed unnecessary
% variables.
% Also compared the results with 3d_util.exe. The boundary elements are the
% same, the circular order (and thus the normal of the surface) of each bel
% is the same as 3d_util.exe.
% 
% Xiaoyao Fan (04/03/2009).


if size(t,2)==3
    edges=[t(:,[1,2]);
           t(:,[1,3]);
           t(:,[2,3])];
elseif size(t,2)==4
    edges = [t(:,[1 2 3]);
             t(:,[2 4 3]);
             t(:,[1 3 4]);
             t(:,[4 2 1])];
else
    error('Cannot handle this type of mesh!');
end
orig=edges;
[edges]=sort(edges,2);
[foo,ix,jx]=unique(edges,'rows');
vec=histc(jx,1:max(jx));
bel=orig(ix(vec==1),:);