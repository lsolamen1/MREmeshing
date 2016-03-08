function [eb pb] = ReadTetraMesh(fn)

if ~isempty(strfind(fn,'.node'))
    [eb pb] = read_nod_elm(fn,1);
elseif ~isempty(strfind(fn,'.mesh'))
    [eb pb] = readMEDIT(fn);
elseif ~isempty(strfind(fn,'.off'))
    [eb pb] = readOFF(fn);
else
    [h pb] = hdrload([fn '.nod']);
    [h eb]=hdrload([fn '.elm']);

    pb = pb(:,2:4)*1000; % in mm
    eb = eb(:,2:5);
end
