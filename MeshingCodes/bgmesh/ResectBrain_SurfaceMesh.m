function [ee pp] = ResectBrain_SurfaceMesh(eb,pb,es,ps)
% Assuming two polyhedrons are defined in A = ['eb' 'pb'] and B = ['es' 'ps']
% It uses CGAL Nef polyhedron concepts to perform a boolean operation which
% subtracts B from A and returns it in C = ['ee' 'pp']: C = A - B;
% 
% Written by:
%  Hamid Ghadyani, Nov 2010

fprintf(' Performing mesh resection.\n');

afilename = [getuserdir filesep 'A.off'];
bfilename = [getuserdir filesep 'B.off'];
dfilename = [getuserdir filesep 'diff.off'];
delete(afilename,bfilename,dfilename);
% Write polyhedrons to OFF files.
writeOFF(afilename,eb,pb);
writeOFF(bfilename,es,ps);

% Write polyhedrons in .mesh format
writenodelm_surface_medit([getuserdir filesep 'a.mesh'],eb,pb,[],0);
writenodelm_surface_medit([getuserdir filesep 'b.mesh'],es,ps,[],0);

syscommand = GetSystemCommand('subtract_polyhedron_cgal');
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
end
% Unset the LD path so that our standalone executables load the system
% libstdc++ rather than Matlab's.
if isunix && ~ismac
    dyldpath = 'unset LD_LIBRARY_PATH; ';
end
subtractcommand = ['! ' dyldpath '"' syscommand '" ' afilename ' ' bfilename ' ' dfilename];
eval(subtractcommand);

% Read the result
[ee pp] = readOFF(dfilename);

if isempty(ee) || isempty(pp)
    error('Resection didn''t produce a proper brain mesh!')
end

% Subdivide the surface to avoid problems related to conversion from Nef
% polyhedron

% syscommand = GetSystemCommand('subdivide_surface_cgal');
% if ~ispc
%     eval(['! chmod u+x "' syscommand '"']);
% end
% subcommand = ['! "' syscommand '" 1 ' dfilename dfilename];
% % eval(subcommand);
% 
% [ee pp] = readOFF(dfilename);
% if isempty(ee) || isempty(pp)
%     error('Resection didn''t produce a proper brain mesh!')
% end

delete(afilename,bfilename,dfilename);
delete('a.mesh','b.mesh');

