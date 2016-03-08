function [outstm]=Output_File_Sem(usedef,strat,meshprop)

junk=pwd;                   % Show (print) current working directory
dirloc=find(junk=='\');     
outstmdef=junk(dirloc(end)+1:end);

if(strat.mesh==1)
    outstmdef=[outstmdef '_voxelmesh'];
elseif(strat.mesh==2)
    outstmdef=[outstmdef '_npw' int2str(meshprop.npw)];
elseif(strat.mesh==3)
    outstmdef=[outstmdef '_res' sprintf('%2.1f-%2.1f-%2.1f',meshprop.targetres.*1000)];
end
    
if(~usedef) 
    outstm=input(['output file stem (default is ' outstmdef '):  >>'],'s');
end
if ~exist('outstm','var')||isempty(outstm)
    outstm=outstmdef;
end

end