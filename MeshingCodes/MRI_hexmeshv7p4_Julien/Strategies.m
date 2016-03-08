function [strat]=Strategies(usedef,default)

% Meshing Approach: one node per voxel, or interpolate to give nodes per wavelength
if(~usedef) 
    disp(['Meshing Strategies: 1=one node per MR voxel, 2=Target #nodes per wavelength, 3=Target Resolution'])
    strat.mesh=input(['Meshing Strategy, (Default ' int2str(default.meshstratdef) ') >>']);
end
if ~exist('meshstrat','var')||isempty(start.mesh)
    strat.mesh=default.meshstratdef;
end

% Zone Sizing Approach: nodes per zone or wavelengths per zone

if(~usedef) 
    disp(['Zone sizing Strategies: 1=fit nodes per zone, 2=fit wavelengths per zone'])
    strat.zone=input(['Zone sizing Strategy, (Default ' int2str(default.zonestratdef) ') >>']);
end
if ~exist('zonestrat','var')||isempty(strat.zone)
    strat.zone=default.zonestratdef;
end
end