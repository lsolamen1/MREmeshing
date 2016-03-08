function gen_MRPEkappa_reconfile(name,dir,itmax,freq,reg,spfilt,udopt,znopt,diag,ovrlp,sim)
%ooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooo
%   GEN_MRPE_RECONFILE.M
%     This routine is used to generate the input file for the poroelastic 
%     reconstruction
%
%   INPUTS
%     name        series name
%     dir         directory for reconstruction files
%     itmax       maximum number of iterations
%     freq        MRE frequency
%     reg         Marquardt regularization parameter
%     spfilt      spatial filtering level
%     udopt       update option
%     znopt       zone option
%     diag        cube diagonal
%     ovrlp       percent zone overlap
%ooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooo

fid = fopen('3D_MRPE_3parm_Files.dat','wt');
fprintf(fid,'3D Subzone Inversion Run File [PARDISO]\n');
fprintf(fid,'Number of OpenMP Processes: <numprocs>\n');
fprintf(fid,'8\n');
fprintf(fid,'Number of Iterations:\n');
fprintf(fid,[itmax,'\n']);
fprintf(fid,'Scaling Factors: <G> <Poisson> <Kp>\n');
fprintf(fid,'3.d3,0.35d0,8.d0\n');
fprintf(fid,'Lower Constraints: <G> <L> <Kp>\n');
fprintf(fid,'1.d0,1.d0,12.d0\n');
fprintf(fid,'Upper Constraints: <G> <L> <Kp>\n');
fprintf(fid,'1.d5,1.d6,5.d0\n');
fprintf(fid,'Frequency (Hz):\n');
fprintf(fid,[num2str(freq,'%.4f'),'\n']);
fprintf(fid,'porosity,density of solid,fluid,apparent:\n');
fprintf(fid,'0.2d0,1020.d0,1000.d0,150.d0\n');
fprintf(fid,'Joachimowicz: <JchG> <JchL> <JchK>\n');
fprintf(fid,'75.d0,75.d0,75.d0\n');
fprintf(fid,'Tikhonov: <TKG> <TKL> <TKK>\n');
fprintf(fid,'0.d0,0.d0,0.d0\n');
fprintf(fid,'Total Variation Regularization Parameter: <TVG> <TVL> <TVK>\n');
fprintf(fid,'0.d0,0.d0,0.d0\n');
fprintf(fid,'Soft Prior Regularization: <SP>\n');
fprintf(fid,'0.d0,0.d0,0.d0\n');
fprintf(fid,'Marquadt Regularization: <MQ>\n');
fprintf(fid,[reg,'\n']);
fprintf(fid,'Spatial Filtering Weight: <theta>\n');
fprintf(fid,[spfilt,'\n']);
fprintf(fid,'Parameter Update Option: <1> global <2> continuous\n');
fprintf(fid,'%d\n',udopt);
fprintf(fid,'Zone Option: <1> cubes <2> rectangular prisms\n');
fprintf(fid,'%d\n',znopt);
fprintf(fid,'Cube Diagonal:\n');
fprintf(fid,'%.4f\n',diag);
fprintf(fid,'Zone Edge Factor: <znedgfac(1:3)>\n');
fprintf(fid,'1 1 1\n');
fprintf(fid,'Zone Overlap Factor: <znovrlp(1:3)>\n');
fprintf(fid,'%.2f %.2f %.2f\n',ovrlp,ovrlp,ovrlp);
fprintf(fid,'Update File Name:\n');
fprintf(fid,'MRPE.upd\n');
fprintf(fid,'File List:\n');
if(sim==1)
    fprintf(fid,['~/meshes/',name,'.mesh/',name,'.nod\n']);
    fprintf(fid,['~/meshes/',name,'.mesh/',name,'.elm\n']);
    fprintf(fid,['~/meshes/',name,'.mesh/',name,'.bel\n']);
else
    fprintf(fid,['../',name,'.nod\n']);
    fprintf(fid,['../',name,'.elm\n']);
    fprintf(fid,['../',name,'.bel\n']);
end
fprintf(fid,[name,'.inv.bcs\n']);
fprintf(fid,[name,'.v3c.ci\n']);
fprintf(fid,['../',name,'.inv.mtr\n']);
fprintf(fid,'pressure.s1c.ci\n');
fprintf(fid,'../regions.dat\n');
fprintf(fid,[dir,'/MRPE_3param-',name,'-inv\n']);
fprintf(fid,[dir,'/3Paramdisplace.v3c\n']);
fprintf(fid,[dir,'/3Parampressure.s1c\n']);
fprintf(fid,[dir,'/3Paramrelerr.dat\n']);
fclose(fid);

%% v11 runfile
fid = fopen('MRPE_v11_runfile.dat','wt');
fprintf(fid,'3D Subzone Inversion Run File [PARDISO]\n');
fprintf(fid,'Number of OpenMP Processes: <numprocs>\n');
fprintf(fid,'8\n');
fprintf(fid,'Number of Iterations:\n');
fprintf(fid,[itmax,'\n']);
fprintf(fid,'Scaling Factors: <G> <Poisson> <Kp>\n');
fprintf(fid,'3.d3,0.35d0,8.d0\n');
fprintf(fid,'Lower Constraints: <G> <L> <Kp>\n');
fprintf(fid,'1.d0,1.d0,12.d0\n');
fprintf(fid,'Upper Constraints: <G> <L> <Kp>\n');
fprintf(fid,'1.d5,1.d6,5.d0\n');
fprintf(fid,'Frequency (Hz):\n');
fprintf(fid,[num2str(freq,'%.4f'),'\n']);
fprintf(fid,'porosity,density of solid,fluid,apparent:\n');
fprintf(fid,'0.2d0,1020.d0,1000.d0,150.d0\n');
fprintf(fid,'Joachimowicz: <JchG> <JchL> <JchK>\n');
fprintf(fid,'75.d0,75.d0,75.d0\n');
fprintf(fid,'Tikhonov: <TKG> <TKL> <TKK>\n');
fprintf(fid,'0.d0,0.d0,0.d0\n');
fprintf(fid,'Total Variation Regularization Parameter: <TVG> <TVL> <TVK>\n');
fprintf(fid,'0.d0,0.d0,0.d0\n');
fprintf(fid,'Soft Prior Regularization: <SP>\n');
fprintf(fid,'0.d0,0.d0,0.d0\n');
fprintf(fid,'Marquadt Regularization: <MQ>\n');
fprintf(fid,[reg,'\n']);
fprintf(fid,'Spatial Filtering Weight: <theta>\n');
fprintf(fid,[spfilt,'\n']);
fprintf(fid,'Parameter Update Option: <1> global <2> continuous\n');
fprintf(fid,'%d\n',udopt);
fprintf(fid,'Zone Option: <1> cubes <2> rectangular prisms\n');
fprintf(fid,'%d\n',znopt);
fprintf(fid,'Cube Diagonal:\n');
fprintf(fid,'%.4f\n',diag);
fprintf(fid,'Zone Edge Factor: <znedgfac(1:3)>\n');
fprintf(fid,'1 1 1\n');
fprintf(fid,'Zone Overlap Factor: <znovrlp(1:3)>\n');
fprintf(fid,'%.2f %.2f %.2f\n',ovrlp,ovrlp,ovrlp);
fprintf(fid,'Update File Name:\n');
fprintf(fid,'MRPE.upd\n');
fprintf(fid,'File List:\n');
if(sim==1)
    fprintf(fid,'Node File:\n');
    fprintf(fid,['~/meshes/',name,'.mesh/',name,'.nod\n']);
    fprintf(fid,'Element File:\n');
    fprintf(fid,['~/meshes/',name,'.mesh/',name,'.elm\n']);
    fprintf(fid,'Boundary Element File:\n');
    fprintf(fid,['~/meshes/',name,'.mesh/',name,'.bel\n']);
else
    fprintf(fid,'Node File:\n');
    fprintf(fid,['../',name,'.nod\n']);
    fprintf(fid,'Element File:\n');
    fprintf(fid,['../',name,'.elm\n']);
    fprintf(fid,'Boundary Element File:\n');
    fprintf(fid,['../',name,'.bel\n']);
end
fprintf(fid,'Boundary Conditon file (pressure only):\n');
fprintf(fid,[name,'.inv.pressure.bcs\n']);
fprintf(fid,'Displacement file (not complex input):\n');
fprintf(fid,['../',name,'.v3c\n']);
fprintf(fid,'inverse material file:\n');
fprintf(fid,['../',name,'.inv.mtr\n']);
fprintf(fid,'Soft Prior regions file:\n');
fprintf(fid,'../regions.dat\n');
fprintf(fid,'output material file:\n');
fprintf(fid,[dir,'/MRPE_3param-',name,'-inv\n']);
fprintf(fid,'output displacement file:\n');
fprintf(fid,[dir,'/3Paramdisplace.v3c\n']);
fprintf(fid,'output pressure file:\n');
fprintf(fid,[dir,'/3Parampressure.s1c\n']);
fprintf(fid,'output error file:\n');
fprintf(fid,[dir,'/3Paramrelerr.dat\n']);
fclose(fid);

end