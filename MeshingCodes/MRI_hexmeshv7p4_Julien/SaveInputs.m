function SaveInputs(outstm,strat,interpo,node,coord,dispscale,dispscalar,default,meshprop,deplacement,msk)
%save data required to get back to original MR voxels (including mask filename)
if isempty(msk)
    msk=namefile.maskoutf;
end

save([outstm '.InterpLocations.mat'],'-struct','interpo','maskint')
save([outstm '.InterpLocations.mat'],'-struct','node','nodin','-append')
save([outstm '.InterpLocations.mat'],'-struct','coord','xout','yout','zout','xin','yin','zin','-append')
save([outstm '.InterpLocations.mat'],'msk','-append')

%save interpolated motion data for possible later use
save([outstm '.InterpData.mat'],'-struct','interpo','maskint','MagIm_int')
save([outstm '.InterpData.mat'],'-struct','deplacement','Urint','Uiint','-append')

%save inputs to meshing code to link with mesh
save([outstm '.meshinput.mat'],'-struct','strat','mesh','zone')
save([outstm '.meshinput.mat'],'dispscale','dispscalar','-append')

if((strat.mesh==2)||(strat.zone==2))
    save([outstm '.meshinput.mat'],'-struct','meshprop','muest')
    save([outstm '.meshinput.mat'],'-struct','default','rhoest','-append')
end
if(strat.mesh==2)||(strat.mesh==3)
    save([outstm '.meshinput.mat'],'-struct','meshprop','targetres','meshres','-append')
end
if(strat.mesh==2)
    save([outstm '.meshinput.mat'],'-struct','meshprop','npw','-append')
end

if(strat.zone==2)
    save([outstm '.meshinput.mat'],'-struct','meshprop','wlperzone','-append')
end
end