function genrunfilev8poro(base,iter,spregions);

nparaminput=input('How many properties do you want to reconstruct? (Default =2)>>');
if nparaminput~=2 || nparaminput~=3;
    nparaminput=2;
end



fid=fopen('runfile_v8poro.dat','wt');


fprintf(fid,'Subzone Reconstruction Data File\n');
fprintf(fid,'Problem Type <0=Forward Problem Only> <1+=Inverse Problem>\n');
fprintf(fid,'1\n');
fprintf(fid,'Node File:\n');
fprintf(fid,['../',base,'.nod\n']);
fprintf(fid,'Element File:\n');
fprintf(fid,['../',base,'.elm\n']);
fprintf(fid,'Boundary Node File: \n');
fprintf(fid,['../',base,'.bnod\n']);
fprintf(fid,'Pressure BC File: (format: <#nbc, bnod, typ, realval, imagval>)ed \n');
fprintf(fid,['../',base,'.pbcs\n']);
fprintf(fid,'Region Stack File:\n');
fprintf(fid,['../',base,'.regstack\n']);
fprintf(fid,'Measured Displacement File:\n');
fprintf(fid,'1\n');

load HeaderData.mat
fprintf(fid,'%10.4fd0',freqHz);fprintf(fid,'\n');
fprintf(fid,['../',base,'.dsp\n']);
fprintf(fid,'Initial Solution File:\n');
fprintf(fid,'0\n');
fprintf(fid,'Output File Stem:\n');
fprintf(fid,['INV/',base,'\n']);
fprintf(fid,'Print Detailed Runtime Debugging and Execution Comments <verb> <file>:\n');
fprintf(fid,'.false.,.true.\n');
fprintf(fid,'Material Model <1=isotropic> <2=orthotropic> <3=iso compress> <4=poroelastic>:\n');
fprintf(fid,'4\n'); 
fprintf(fid,'Number of Material Properties:\n');


fprintf(fid,'4 \n');
fprintf(fid,'Material Description Style <1=nodal> <2=element> [<shear modulus> <density> <bulk modulus>]: \n');
fprintf(fid,'1,1,1,1  \n');
fprintf(fid,'Number of Parameters per Material Point: \n');
fprintf(fid,'1,1,1,1 \n');
fprintf(fid,'Property Scalars <Shear Mod, Lamda Mod, Hydraulic conductivity, porosity>: \n');
fprintf(fid,' 3300.d0, 100.d0, 10000.d0, 1.d-4, -8.d0, -24.d0, 0.2d0, 1.d-4\n');
fprintf(fid,'Other Poroealstic constants assumed constant <rho, rhof, rhoa>\n');
fprintf(fid,' 1020.d0, 1000.d0, 150.d0\n');
fprintf(fid,'Multifrequency Indicator\n');
fprintf(fid,'.false.,.false.,.false.,.false.\n');
fprintf(fid,'Multifrequency Alpha\n');
fprintf(fid,'0.d0,0.d0,0.d0,0.d0\n');
fprintf(fid,'0.d0,0.d0,0.d0,0.d0\n');
fprintf(fid,'Number of Material Property Mesh Resolutions: \n');
fprintf(fid,'3\n');
fprintf(fid,'Number of material mesh structures\n');
fprintf(fid,'3\n');
fprintf(fid,'Matrial mesh iteration stucture (one less than number of structures)\n');
fprintf(fid,'3,6\n');
fprintf(fid,'Material Property Mesh Resolutions (x,y,z): \n');
fprintf(fid,'0.0025000 , 0.0025000 , 0.0025000 \n');
fprintf(fid,'0.0125000 , 0.0125000 , 0.0125000 \n');
fprintf(fid,'0.0250000 , 0.0250000 , 0.0250000 \n');
fprintf(fid,' Material Mesh Index (1st line real part, 2nd line imag)  \n');
fprintf(fid,' 1 1 1 1\n');
fprintf(fid,' 3 3 3 3\n');
fprintf(fid,' Material Mesh Index (1st line real part, 2nd line imag) \n');
fprintf(fid,' 1 1 1 1\n');
fprintf(fid,' 3 3 3 3\n');
fprintf(fid,' Material Mesh Index (1st line real part, 2nd line imag) \n');
fprintf(fid,' 1 1 1 1\n');
fprintf(fid,' 3 3 3 3\n');
fprintf(fid,'Reconstruction Indicators: \n');


if nparaminput==2; 
fprintf(fid,'.true.,.false. \n');
fprintf(fid,'.true.,.false. \n');
fprintf(fid,'.false.,.false. \n');
fprintf(fid,'.false.,.false. \n');

elseif nparaminput==3;
fprintf(fid,'.true.,.false. \n');
fprintf(fid,'.true.,.false. \n');
fprintf(fid,'.true.,.false. \n');
fprintf(fid,'.false.,.false. \n');
    
else
    error('# of parameters to reconstruct is not clear');
end

fprintf(fid,'Property Estimate Variance Calculation: \n');
fprintf(fid,'.false. \n');
fprintf(fid,'Zone Sizes (not including overlap factor, actual size = (1+2*ovlp)*siz) [x y z]: \n');
fprintf(fid,' 1.95633e-02, 1.95633e-02, 1.95633e-02\n');
fprintf(fid,'Zone Overlap Percentage [x y z]: \n');
fprintf(fid,'0.150d0,0.150d0,0.150d0\n');
fprintf(fid,'Iteration Limits [global zone]: \n');
fprintf(fid,[iter,',1\n']);
fprintf(fid,'Minimum Parameter Update Size [global zone line]: \n');
fprintf(fid,'0.1d0, 0.1d0, 0.1d0\n');
fprintf(fid,'Number of Zone Iteration Structures: \n');
fprintf(fid,'3\n'); 
fprintf(fid,'Iteration Limits for Zone Iteration Structures (NOTE:  provide one less limit than # of structures!!!): \n');


if isequal(iter,'100');
    fprintf(fid,'50,80 \n');
elseif isequal(iter,'200');
    fprintf(fid,'20,150 \n');
else 
    fprintf(fid,'50,80 \n');
end

fprintf(fid,'Zone Iteration Structures [<# of CG iters> <# of GN iters> <!!! QN CURRENTLY UNAVAILABLE !!!> <# of line search iters>]: \n');
fprintf(fid,'1,0,0,1 \n');
fprintf(fid,'2,0,0,2 \n');
fprintf(fid,'3,0,0,3 \n');
fprintf(fid,'Number of Processors per MUMPS Communicator \n');
fprintf(fid,'1 \n');
fprintf(fid,'Maximum Amount of RAM per Processor [MB] \n');
fprintf(fid,'999999999 \n');
fprintf(fid,'ooooooooooo  REGULARIZATION DESCRIPTORS  ooooooooooooooooooooooooooooooooooooooooooooooooooooooooooo \n');
fprintf(fid,'Regularization Indicators [<TV> <SF> <TK> <MQ> <JW> <constraints> <CG residual> <VH> <McG> <soft prior>]: \n');


if spregions==1
  fprintf(fid,'.false.,.true.,.false.,.false.,.false.,.false.,.true.,.true.,.true.,.true. \n');
else 
fprintf(fid,'.false.,.true.,.false.,.false.,.false.,.false.,.true.,.true.,.true.,.false. \n');
end

fprintf(fid,'Number of constant regularization iterations (final N iterations use final regularization weights) \n');
fprintf(fid,'30 \n');
fprintf(fid,'Number of Parameters to Treat with Total Variation: \n');
fprintf(fid,'4 \n');
fprintf(fid,'TV Parameter List: \n');
fprintf(fid,'1,2,3,4 \n');
fprintf(fid,'TV Delta Values: <1st line real, 2nd line imag> \n');
fprintf(fid,'1.d-19,1.d-19,1.d-19 ,1.d-19\n');
fprintf(fid,'1.d-19,1.d-19,1.d-19 ,1.d-19\n');
fprintf(fid,'TV Initial Weights: <1st line real, 2nd line imag> \n');
fprintf(fid,'5.d-16,5.d-16,5.d-16 ,5.d-16\n');
fprintf(fid,'5.d-16,5.d-16,5.d-16 ,5.d-16\n');
fprintf(fid,'TV Final Weights: <1st line real, 2nd line imag> \n');
fprintf(fid,'5.d-16,5.d-16,5.d-16 ,5.d-16\n');
fprintf(fid,'5.d-16,5.d-16,5.d-16 ,5.d-16\n');
fprintf(fid,'Number of Parameters to Treat with Spatial Filtering: \n');
fprintf(fid,'4 \n');
fprintf(fid,'SF Parameter List: \n');
fprintf(fid,'1,2,3,4 \n');
fprintf(fid,'SF Gaussian filter initial Widths: <1st line real, 2nd line imag> \n');
fprintf(fid,'0.003d0,0.003d0,0.003d0,0.003d0 \n');
fprintf(fid,'0.003d0,0.003d0,0.003d0,0.003d0 \n');
fprintf(fid,'SF Final Gaussian width: <1st line real, 2nd line imag> \n');
fprintf(fid,'0.0015d0,0.0015d0,0.0015d0,0.0015d0 \n');
fprintf(fid,'0.0015d0,0.0015d0,0.0015d0,0.0015d0 \n');
fprintf(fid,'Number of Parameters to Treat with Tikhonov Regularization: \n');
fprintf(fid,'4 \n');
fprintf(fid,'TK Parameter List: \n');
fprintf(fid,'1,2,3,4 \n');
fprintf(fid,'TK Initial Weights: <1st line real, 2nd line imag> \n');
fprintf(fid,'1.d-18,1.d-18,1.d-18,1.d-18 \n');
fprintf(fid,'1.d-18,1.d-18,1.d-18,1.d-18 \n');
fprintf(fid,'TK Final Weights: <1st line real, 2nd line imag> \n');
fprintf(fid,'1.d-18,1.d-18,1.d-18,1.d-18\n');
fprintf(fid,'1.d-18,1.d-18,1.d-18,1.d-18\n');
fprintf(fid,'Number of Parameters to Treat with Marquardt Regularization: \n');
fprintf(fid,'4 \n');
fprintf(fid,'Distance of alpha from 1 at which MQ weights are adjusted: \n');
fprintf(fid,'0.25d0 \n');
fprintf(fid,'MQ Parameter List: \n');
fprintf(fid,'1,2,3,4 \n');
fprintf(fid,'MQ Weight Delta: <1st line real, 2nd line imag> \n');
fprintf(fid,'0.5d0,0.5d0,0.5d0,0.5d0 \n');
fprintf(fid,'0.5d0,0.5d0,0.5d0,0.5d0 \n');
fprintf(fid,'MQ Initial Weights: <1st line real, 2nd line imag> \n');
fprintf(fid,'1.d2,1.d2,1.d2,1.d2 \n');
fprintf(fid,'1.d2,1.d2,1.d2,1.d2 \n');
fprintf(fid,'MQ Minimum Weights: <1st line real, 2nd line imag> \n');
fprintf(fid,'1.d-11,1.d-11,1.d-11,1.d-11\n');
fprintf(fid,'1.d-11,1.d-11,1.d-11,1.d-11 \n');
fprintf(fid,'Number of Parameters to Treat with Joachimowicz Regularization: \n');
fprintf(fid,'4\n');
fprintf(fid,'JW Parameter List: \n');
fprintf(fid,'1,2,3,4 \n');
fprintf(fid,'JW Initial Weights: <1st line real, 2nd line imag> \n');
fprintf(fid,'5.d0,5.d0,5.d0,5.d0\n');
fprintf(fid,'5.d0,5.d0,5.d0,5.d0\n');
fprintf(fid,'JW Final Weights: <1st line real, 2nd line imag> \n');
fprintf(fid,'5.d0,5.d0,5.d0,5.d0 \n');
fprintf(fid,'5.d0,5.d0,5.d0,5.d0\n');
fprintf(fid,'Number of Parameters to Treat with Constraints: \n');
fprintf(fid,'4\n');
fprintf(fid,'Constraint Parameter List: \n');
fprintf(fid,'1,2,3,4 \n');
fprintf(fid,'Constraint Weights <converges to exact solution as weight --> inf>: \n');
fprintf(fid,'1.d-14,1.d-14,1.d-14,1.d-14 \n');
fprintf(fid,'Constraint Minimums: \n');
fprintf(fid,'200.d0, 0.d0, 200.d0, 0.d0, -11.d0, -24.d0, 0.01d0, 0.d0 \n');
fprintf(fid,'Constraint Maximums: \n');
fprintf(fid,'5.d5, 5.d5, 5.d5, 5.d5, -5.d0, -24.d0, 0.8d0, 0.d0 \n');
fprintf(fid,'CG Residual Scaling Initial Weights: <1st line real, 2nd line imag> \n');
fprintf(fid,'0.04d0,0.04d0,0.04d0,0.04d0 \n');
fprintf(fid,'0.04d0,0.04d0,0.04d0,0.04d0 \n');
fprintf(fid,'CG Residual Scaling Final Weights: <1st line real, 2nd line imag> \n');
fprintf(fid,'0.01d0,0.01d0,0.01d0,0.01d0 \n');
fprintf(fid,'0.01d0,0.01d0,0.01d0,0.01d0 \n');
fprintf(fid,'Van Houten Regularization Level: \n');
fprintf(fid,'1.2d0 \n');
fprintf(fid,'Number of parameters to treat with soft prior: \n');
fprintf(fid,'4 \n');
fprintf(fid,'SP Parameter List: \n');
fprintf(fid,'1,2,3,4 \n');
fprintf(fid,'SP Initial Weights: <1st line real, 2nd line imag> \n');
fprintf(fid,'1.d-22,1.d-22,1.d-22,1.d-22 \n');
fprintf(fid,'1.d-22,1.d-22,1.d-22,1.d-22 \n');
fprintf(fid,'SP Final Weights: <1st line real, 2nd line imag> \n');
fprintf(fid,'1.d-22,1.d-22,1.d-22,1.d-22 \n');
fprintf(fid,'1.d-22,1.d-22,1.d-22,1.d-22 \n');
fprintf(fid,'SP start iteration delay (first N iterations do not use soft prior) \n');
fprintf(fid,'0 \n');

fclose(fid);

fid2=fopen('MPICH2-Submitv8poro','wt');
fprintf(fid2,'#!/bin/bash -l \n');
fprintf(fid2,'# declare a name for this job \n');
fprintf(fid2,'#PBS -N MRE3D_150415_continuous_voxelmesh_G3300.v7.3.inv.iso.incomp.visc_SPoff_SF0p0015_CG2p2\n');
fprintf(fid2,'# request the queue (enter the possible names, if omitted, serial is the default) \n');
fprintf(fid2,'#PBS -q default  \n');
fprintf(fid2,'# request  node  \n');
fprintf(fid2,'#PBS -l nodes=4:ppn=8 \n');
fprintf(fid2,'#PBS -l feature=amd \n');
fprintf(fid2,'# request some hours of wall time  \n');
fprintf(fid2,'#PBS -l walltime=36:00:00  \n');
fprintf(fid2,'#combine PBS standard output and error files  \n');
fprintf(fid2,'#PBS -j oe  \n');
fprintf(fid2,'# mail is sent to you when the job starts and when it terminates or aborts  \n');
fprintf(fid2,'##PBS -m bea  \n');
fprintf(fid2,'# specify your email address  \n');
fprintf(fid2,'##PBS -M matthew.d.mcgarry@dartmouth.edu  \n');
fprintf(fid2,'# Change to Submission Directory  \n');
fprintf(fid2,'cd $PBS_O_WORKDIR  \n');
fprintf(fid2,'# run the program \n');
fprintf(fid2,'cat $PBS_NODEFILE | uniq > node_file \n');
fprintf(fid2,'nodes=$(cat node_file | wc -l) \n');
fprintf(fid2,'nprocs=$(cat $PBS_NODEFILE | wc -l) \n');
fprintf(fid2,'export MKL_NUM_THREADS=1 \n');
fprintf(fid2,'echo Nodes $nnodes \n');
fprintf(fid2,'echo Procs $nprocs \n');
fprintf(fid2,'mpiexec -n $nprocs -hostfile $PBS_NODEFILE /ihome/lsolamen/code/MREv8/MRE-Zone.mp2 runfile_v8poro.dat\n');
fprintf(fid2,'exit 0  \n');
fclose(fid2);



