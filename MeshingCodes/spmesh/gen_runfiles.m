function gen_runfiles(add,nprocs,name)
%ooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooo
%   GEN_RUNFILES.M
%
%   DESCRIPTION:
%     GEN_RUNFILES(add,nprocs,name) generates a PBS jobfile for the 
%     parallelMRE reconstructions for use on the DISCOVERY cluster
%
%   INPUT:
%     add         a string containing the email address of the operator
%     nprocs      an integer representing the total number of cpus 
%                 requested
%     name        series name
%
%   DEFAULTS:
%     walltime - 48 hrs
%     output/error files - combined
%
%   Phillip R. Perriñez, Ph.D.
%   Thayer School of Engineering
%   Dartmouth College
%   July 2009
%ooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooo

%compute number of nodes required
nnodes=(str2num(nprocs)/8);
job=['#PBS -N ',name];
CPU =['#PBS -l nodes=',num2str(nnodes),':','ppn=8'];

%% ELASTIC (PARDISO)
CPU =['#PBS -l nodes=',num2str(nnodes),':','ppn=8:amd'];
fid=fopen('MREpar-submit', 'wt');
fprintf(fid,'%s \n','# declare a name for this job to be sample_job');
fprintf(fid,'%s \n',job);
fprintf(fid,'%s \n','# request the queue (enter the possible names, if omitted, serial is the default)');
fprintf(fid,'%s \n','#PBS -q default');
fprintf(fid,'%s \n','# request  node');
fprintf(fid,'%s \n',CPU);
fprintf(fid,'%s \n','# request 48 hours of wall time');
fprintf(fid,'%s \n','#PBS -l walltime=48:00:00');
fprintf(fid,'%s \n','#combine PBS standard output and error files');
fprintf(fid,'%s \n','#PBS -j oe');
fprintf(fid,'%s \n','# mail is sent to you when the job starts and when it terminates or aborts');
fprintf(fid,'%s \n','#PBS -m bea');
fprintf(fid,'%s \n','# specify your email address');
fprintf(fid,'%s \n', ['#PBS -M',' ',add]);
fprintf(fid,'%s \n','# By default, PBS scripts execute in your home directory, not the');
fprintf(fid,'%s \n','# directory from which they were submitted. The following line');
fprintf(fid,'%s \n','# places you in the directory from which the job was submitted.');
fprintf(fid,'%s \n','cd $PBS_O_WORKDIR');
fprintf(fid,'%s \n','# run the program');
fprintf(fid,'%s \n','mpirun -machinefile $PBS_NODEFILE -np 16 /home/apattiso/src/MRE_INV/MRE-pardiso.x');
fprintf(fid,'%s \n','mpiexec -pernode -comm=none /opt/mpich/mpich1.2.7p1-intel11.0/sbin/cleanipcs');
fprintf(fid,'%s \n','exit 0');
fclose(fid);

%% POROELASTIC 2parm (PARDISO)
fid=fopen('MRPE-2prm-submit', 'wt');
fprintf(fid,'%s \n','#!/bin/bash -l');
fprintf(fid,'%s \n','# declare a name for this job to be sample_job');
fprintf(fid,'%s \n',['#PBS -N 2P_',name]);
fprintf(fid,'%s \n','# request the queue (enter the possible names, if omitted, serial is the default)');
fprintf(fid,'%s \n','#PBS -q default');
fprintf(fid,'%s \n','# request  node');
fprintf(fid,'%s \n','#PBS -l nodes=2:ppn=8');
fprintf(fid,'%s \n','#PBS -l feature=amd');
fprintf(fid,'%s \n','# request 24 hours of wall time');
fprintf(fid,'%s \n','#PBS -l walltime=24:00:00');
fprintf(fid,'%s \n','#combine PBS standard output and error files');
fprintf(fid,'%s \n','#PBS -j oe');
fprintf(fid,'%s \n','# mail is sent to you when the job starts and when it terminates or aborts');
fprintf(fid,'%s \n','##PBS -m bea');
fprintf(fid,'%s \n','## specify your email address');
fprintf(fid,'%s \n', ['##PBS -M ',add]);
fprintf(fid,'%s \n','# By default, PBS scripts execute in your home directory, not the');
fprintf(fid,'%s \n','# directory from which they were submitted. The following line');
fprintf(fid,'%s \n','# places you in the directory from which the job was submitted.');
fprintf(fid,'%s \n','cd $PBS_O_WORKDIR');
fprintf(fid,'%s \n','# run the program');
fprintf(fid,'%s \n','export MKL_NUM_THREADS=1');
fprintf(fid,'%s \n','cat $PBS_NODEFILE | uniq > node_file ');
fprintf(fid,'%s \n','nnodes=$(cat node_file | wc -l) ');
fprintf(fid,'%s \n','nprocs=$(cat $PBS_NODEFILE | wc -l)');
fprintf(fid,'%s \n','echo Nodes $nnodes'); 
fprintf(fid,'%s \n','echo Procs $nprocs');
fprintf(fid,'%s \n',['mpirun -np $nprocs -hostfile $PBS_NODEFILE /ihome/mmcgarry/code/poroMRE/2Precon/MRPE-2Paramrecon_v4.x']);
fprintf(fid,'%s \n','exit 0');
fclose(fid);

%% POROELASTIC 3parm (PARDISO)
fid=fopen('MRPE-3prm-submit', 'wt');
fprintf(fid,'%s \n','#!/bin/bash -l');
fprintf(fid,'%s \n','# declare a name for this job to be sample_job');
fprintf(fid,'%s \n',['#PBS -N 3P_',name]);
fprintf(fid,'%s \n','# request the queue (enter the possible names, if omitted, serial is the default)');
fprintf(fid,'%s \n','#PBS -q default');
fprintf(fid,'%s \n','# request  node');
fprintf(fid,'%s \n','#PBS -l nodes=2:ppn=8');
fprintf(fid,'%s \n','#PBS -l feature=amd');
fprintf(fid,'%s \n','# request 24 hours of wall time');
fprintf(fid,'%s \n','#PBS -l walltime=24:00:00');
fprintf(fid,'%s \n','#combine PBS standard output and error files');
fprintf(fid,'%s \n','#PBS -j oe');
fprintf(fid,'%s \n','# mail is sent to you when the job starts and when it terminates or aborts');
fprintf(fid,'%s \n','##PBS -m bea');
fprintf(fid,'%s \n','## specify your email address');
fprintf(fid,'%s \n', ['##PBS -M ',add]);
fprintf(fid,'%s \n','# By default, PBS scripts execute in your home directory, not the');
fprintf(fid,'%s \n','# directory from which they were submitted. The following line');
fprintf(fid,'%s \n','# places you in the directory from which the job was submitted.');
fprintf(fid,'%s \n','cd $PBS_O_WORKDIR');
fprintf(fid,'%s \n','# run the program');
fprintf(fid,'%s \n','# run the program'); 
fprintf(fid,'%s \n','export MKL_NUM_THREADS=1'); 
fprintf(fid,'%s \n','cat $PBS_NODEFILE | uniq > node_file');  
fprintf(fid,'%s \n','nnodes=$(cat node_file | wc -l)');  
fprintf(fid,'%s \n','nprocs=$(cat $PBS_NODEFILE | wc -l)'); 
fprintf(fid,'%s \n','echo Nodes $nnodes'); 
fprintf(fid,'%s \n','echo Procs $nprocs'); 
fprintf(fid,'%s \n','mpirun -np $nprocs -hostfile $PBS_NODEFILE /ihome/mmcgarry/code/poroMRE/AJ3Preconv10/MRPE-3Precon_v10.x'); 
fprintf(fid,'%s \n','exit 0');
fclose(fid);

%% POROELASTIC 3parmv11 (PARDISO)
fid=fopen('MRPE-3prm-v11-submit', 'wt');
fprintf(fid,'%s \n','#!/bin/bash -l');
fprintf(fid,'%s \n','# declare a name for this job to be sample_job');
fprintf(fid,'%s \n',['#PBS -N 3P_',name]);
fprintf(fid,'%s \n','# request the queue (enter the possible names, if omitted, serial is the default)');
fprintf(fid,'%s \n','#PBS -q default');
fprintf(fid,'%s \n','# request  node');
fprintf(fid,'%s \n','#PBS -l nodes=2:ppn=8');
fprintf(fid,'%s \n','#PBS -l feature=amd');
fprintf(fid,'%s \n','# request 30 hours of wall time');
fprintf(fid,'%s \n','#PBS -l walltime=30:00:00');
fprintf(fid,'%s \n','#combine PBS standard output and error files');
fprintf(fid,'%s \n','#PBS -j oe');
fprintf(fid,'%s \n','## mail is sent to you when the job starts and when it terminates or aborts');
fprintf(fid,'%s \n','##PBS -m bea');
fprintf(fid,'%s \n','## specify your email address');
fprintf(fid,'%s \n', ['##PBS -M ',add]);
fprintf(fid,'%s \n','# By default, PBS scripts execute in your home directory, not the');
fprintf(fid,'%s \n','# directory from which they were submitted. The following line');
fprintf(fid,'%s \n','# places you in the directory from which the job was submitted.');
fprintf(fid,'%s \n','cd $PBS_O_WORKDIR');
fprintf(fid,'%s \n','# run the program');
fprintf(fid,'%s \n','export MKL_NUM_THREADS=1');
fprintf(fid,'%s \n','cat $PBS_NODEFILE | uniq > node_file ');
fprintf(fid,'%s \n','nnodes=$(cat node_file | wc -l) ');
fprintf(fid,'%s \n','nprocs=$(cat $PBS_NODEFILE | wc -l)');
fprintf(fid,'%s \n','echo Nodes $nnodes'); 
fprintf(fid,'%s \n','echo Procs $nprocs');
fprintf(fid,'%s \n',['mpirun -np $nprocs -hostfile $PBS_NODEFILE /ihome/mmcgarry/code/poroMRE/HCrecon/MRPE-3P-inverse_v11.discov MRPE_v11_runfile.dat']);
fprintf(fid,'%s \n','exit 0');
fclose(fid);

end