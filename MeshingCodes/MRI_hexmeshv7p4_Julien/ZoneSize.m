function [zoneprop]=ZoneSize(node,strat,meshprop)
% Zonestrat=1: target number of nodes per zone
% Zonestrat=2: target number of wavelengths per zone
zoneprop.znovlp=[0.15 0.15 0.15]; % Zone overlap
znovlp1=1+zoneprop.znovlp;
%Calculate edge lenth factors
rangedim=max(node.nod)-min(node.nod);
disp(['Mesh Size Ratio (x,y,z): ',num2str(rangedim./min(rangedim))])

if(strat.zone==1)
    
    % -> Aim for approximatly cubic zones, with a specified number of nodes per zone
    % npz = Target number of nodes per zone (Note that this strategy usually ends up with zones ~70% of target size).
    disp(['zonestrat=1, aiming for ' int2str(meshprop.npz) ' nodes per zone'])
    zfac=rangedim./min(rangedim); % Ratio of mesh dimensions
    
    % Nn/Nz=Npz/(znovlp+1)^3
    % Nz=Nn*(znovrlp+1)^3/Npz
    % zx*zy*zz = Nz
    % zfacx*F * zfacy*F * zfacz*F = Nz
    % F = (Nz/(zfacx*zfacy*zfacz))^(1/3)
    % zx =zfacx*F, zy=zfacy*F, zz=zfacz*Ffprintf(fid,'%c',direc);fprintf(fid,'/');fprintf(fid,'%c',elmf);fprintf(fid,'/n')
    
    Nz=node.nn*znovlp1(1)*znovlp1(2)*znovlp1(3)/meshprop.npz;
    F=(Nz/(zfac(1)*zfac(2)*zfac(3)))^(1/3);
    znedge=zfac.*F;
    
    znedgeint=round(znedge);
    for ii=1:3
        if znedgeint(ii)==0
            znedgeint(ii)=1;
        end
    end
    disp(['v <= 7.04 Zone edge factors calculated = ' int2str(znedgeint)])
    disp(['Estimated nodes per zone = ' num2str(node.nn*prod(znovlp1)/(prod(znedgeint)))])
    
    % v7.05 = direct specification of zone sizes, [L1 L2 L3]
    % Aim for cubic zones, NPZ =(L*(1+2*ovlp(1)))/meshres(1))*(L*(1+2*ovlp(2)))/meshres(2))*(L*(1+2*ovlp(3)))/meshres(3))
    % L^3=NPZ*mesres(1)*meshres(2)*meshres(3)/(1+2*ovlp(1))*(1+2*ovlp(2))*1+2*ovlp(3)))
    zoneprop.znedgelength(1:3)=( meshprop.npz*meshprop.meshres(1)*meshprop.meshres(2)*meshprop.meshres(3) / ((1+2*zoneprop.znovlp(1))*(1+2*zoneprop.znovlp(2))*(1+2*zoneprop.znovlp(3))) )^(1/3);
    disp(['v7.05 zone edge length ' num2str(zoneprop.znedgelength)])
    disp(['Estimated nodes per zone = ' int2str(prod(zoneprop.znedgelength.*(1+2*zoneprop.znovlp)./meshprop.meshres))])
elseif(strat.zone==2)
    disp(['zonestrat=2, aiming for ' int2str(meshprop.wlperzone) ' wavelengths per zone'])
    
    znesz=meshprop.wlperzone*meshprop.wl;
    znedgeint=round(rangedim./znesz.*(1+2*zoneprop.znovlp));
    for ii=1:3
        if(znedgeint(ii)<3) % Mostly Edge zones, only one SZ overlap
            znedgeint(ii)=round(rangedim(ii)./znesz.*(1+zoneprop.znovlp(ii)));
        end
        if znedgeint(ii)==0
            znedgeint(ii)=1;
        end
    end
    
    disp(['v <= 7.04 Zone edge factors calculated = ' int2str(znedgeint)])
    disp(['Estimated nodes per zone = ' num2str(node.nn*prod(znovlp1)/(prod(znedgeint)))])
    
    % v7.05 = direct specification of zone sizes, [L1 L2 L3]
    % Actual SZ length = (1+2*znovlp(i))*L(i)
    for ii=1:3
        zoneprop.znedgelength(ii)=znesz/(1+2*zoneprop.znovlp(ii));
    end
    disp(['v7.05 zone edge length ' num2str(zoneprop.znedgelength)])
    disp(['Estimated nodes per zone = ' num2str(prod(znesz./meshprop.meshres)) ])
end

disp(['Meshing Complete, ' int2str(node.nn) ' Nodes and ' int2str(node.nel) ' elements'])
end