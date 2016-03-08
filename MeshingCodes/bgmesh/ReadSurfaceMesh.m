function [es ps] = ReadSurfaceMesh(fn)

if ~isempty(strfind(fn,'.off'))
    [es ps] = readOFF(fn);
else
    [h es] = hdrload([fn '.elm']);
    [h ps] = hdrload([fn '.nod']);
    ps = ps(:,2:4);
    es = es(:,2:4);
end

    