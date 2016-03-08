function mesh2vtk(elm,node, tfn);
% function mesh2vtk(elm,node, tfn);
% converts element and node listings into vtk format and save to file
% (tfn). Only for 3D triangular mesh for now.
%
% Songbai Ji 6/24/2006

if nargin ~= 3
    error('Input argument must be 3');
end

fid = fopen (tfn, 'w');
if fid == -1
    error([tfn, ' cannot be opened for write!']);
end
nPnts = size(node,1);
fprintf(fid, ['# vtk DataFile Version 3.0', char(10), ...
        'vtk output', char(10), ...
        'ASCII', char(10), ...
        'DATASET POLYDATA', char(10), ...
        'POINTS ', num2str(nPnts), ' float', char(10)]);
nLines = floor(nPnts/3);  %% vtk node point each row has 3 (or less, for the last row) pnts
% that makes up 9 coordinates (for 3D).
i=1;
for j = 1 : nLines
    % 3(points)*3(coordinates:x,y,z)=9
    string = [sprintf('%g %g %g %g %g %g %g %g %g\n', node(i:i+2,:)')];
    % commented is slow
    %         fprintf(fid, [num2str( [node(i,:), node(i+1, :), node(i+2,:)] ), char(10)]);
    fprintf(fid,string);
    i = i+3;
end
remainder = rem(nPnts,3);
for j = nPnts-remainder+1 : nPnts
    fprintf(fid, num2str( [node(j,:)])); 
end
fprintf(fid, char(10));

% now the elements
nElms = length(elm(:,1));
fprintf(fid, ['POLYGONS ', num2str(nElms), ' ', num2str(nElms*4), char(10)]);
elm = elm-1; %vtk starts index from 0
for i = 1 : nElms
    string = [sprintf('%g %g %g %g\n', 3, elm(i,:)')];
    % commented is slow
    %         fprintf(fid, ['3 ', num2str(elm(i,:)), char(10)]);
    fprintf(fid,string);
end

fprintf(fid, ['CELL_DATA ', sprintf('%g\n', nElms)]);
fprintf(fid, ['POINT_DATA ',sprintf('%g\n', nPnts)]);
fclose(fid);