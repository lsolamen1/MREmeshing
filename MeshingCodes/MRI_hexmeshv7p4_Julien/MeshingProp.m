function [meshprop]=MeshingProp(usedef,default,freqHz,strat)

if(strat.mesh==1)
    meshprop.wl=0;
    meshprop.npw=0;
    meshprop.targetres=0;
    meshprop.muest=default.mudef;
end
if(strat.mesh==2)||(strat.zone==2)
    % Estimated properties.    
    if(~usedef) 
        meshprop.muest=input(['Estimate of Shear Modulus (Default ' num2str(default.mudef) ') >>']);
    end
    if ~exist('meshprop.muest','var')||(isempty(meshprop.muest))
        meshprop.muest=default.mudef;
    end    

    meshprop.wl=1/freqHz*sqrt(meshprop.muest/default.rhoest); % Wavelength
end

if(strat.mesh==2)
    % Desired nodes per wavelength.
    if(~usedef) 
        meshprop.npw=input(['Desired Nodes per wavelength: (Default ' num2str(default.npwdef) ') >>']);
    end
    if ~exist('meshprop.npw','var')||(isempty(meshprop.npw))
        meshprop.npw=default.npwdef;
    end
    meshprop.targetres=0;
end

if(strat.zone==1)
    % Nodes per zone
    if(~usedef) 
        meshprop.npz=input(['Target Number of Nodes per zone (Default ' num2str(default.npzdef) ') >>']);
    end
    if ~exist('meshprop.npz','var')||(isempty(meshprop.npz))
        meshprop.npz=default.npzdef; 
    end 
    meshprop.wlperzone=0;
elseif(strat.zone==2)
    %Wavelengths per zone
    if(~usedef) 
        meshprop.wlperzone=input(['Target Number of wavelengths per zone (Default ' num2str(default.wlperzonedef) ') >>']);
    end
    if ~exist('meshprop.wlperzone','var')||(isempty(meshprop.wlperzone))
        meshprop.wlperzone=default.wlperzonedef;  
    end 
    meshprop.npz=0;
end

if(strat.mesh==3)
    if(~usedef) 
        meshprop.targetres=input(['Target Resolution (Default ' num2str(default.resdef) ') >>']);
    end
    if ~exist('meshprop.targetres','var')||(isempty(meshprop.targetres))
        meshprop.targetres=default.resdef;
    end
    meshprop.wl=0;
    meshprop.npw=0;
    meshprop.muest=default.mudef;
end
    
end
