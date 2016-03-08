function WriteOutputFiles(fn,e,p,bel)
% Writes the tetrahedral mesh in 'e' and 'p' to files:
% fn.elm, fn.nod and fn.bel
% 
% Written by:
%  Hamid Ghadyani, Nov 2010

fprintf(' Writing output files.\n');

if nargin~=4
    bel = getBdyFromMesh(e);
end

bel = [(1:length(bel))' bel];
e = [(1:length(e))' e];
p = p/1000;
p = [(1:length(p))' p];

fid = OpenFile([fn '.elm'], 'wt');

fprintf (fid, '# number of elements: %d\n', size(e,1));
fprintf (fid, '%d\t%d\t%d\t%d\t%d\t1\n', e');
fclose(fid);

fid = OpenFile([fn '.nod'], 'wt');

fprintf (fid, '# number of nodes: %d\n', size(p,1));
fprintf (fid, '%d\t%.18f\t%.18f\t%.18f\n', p');
fclose(fid);

fid = OpenFile([fn '.bel'], 'wt');

fprintf (fid, '# number of bels: %d\n', size(bel,1));
fprintf (fid, '%d\t%d\t%d\t%d\t1\t0\n', bel');
fclose(fid);
