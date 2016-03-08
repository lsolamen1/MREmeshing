function gen_MREpar_runfile(email,numprocs,series_name)
%*************************************************************************
%  GEN_RUNFILE.M
%
%  DESCRIPTION:
%   GEN_MREPAR_RUNFILE(email,numprocs) generates a PBS jobfile for the 
%   MRE elasticity-based parallel reconstruction for use on the DISCOVERY
%   cluster.
%
%  INPUT:
%   email - a string containing the email address of the operator.
%   proc - an integer representing the total number of cpus requested.
%
%  DEFAULTS:
%   walltime - 48 hrs
%   output/error files - combined
%
%  OUTPUT:
%   Generates a PBS job file named parallel_recon
%
%  PHILLIP R. PERRINEZ, PHD
%  Thayer School of Engineering
%  Dartmouth College
%  November 2008
%************************************************************************

npro=(numprocs/8);
name=['#PBS -N ',series_name];
CPU =['#PBS -l nodes=',num2str(npro),':','ppn=8:quadcore'];

fid=fopen('MREpar-submit', 'w');
fprintf(fid,'%s \n','# declare a name for this job to be sample_job');
fprintf(fid,'%s \n',name);
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
fprintf(fid,'%s \n', ['#PBS -M',' ',email]);
fprintf(fid,'%s \n','# By default, PBS scripts execute in your home directory, not the'); 
fprintf(fid,'%s \n','# directory from which they were submitted. The following line'); 
fprintf(fid,'%s \n','# places you in the directory from which the job was submitted.');
fprintf(fid,'%s \n','cd $PBS_O_WORKDIR');
fprintf(fid,'%s \n','# run the program');
fprintf(fid,'%s \n','mpirun -machinefile $PBS_NODEFILE -np 16 /home/pperrine/MRE-PARDISO/MRE-pardiso.x');
fprintf(fid,'%s \n','mpiexec -pernode -comm=none /opt/mpich/mpich1.2.7p1-intel9.1/sbin/cleanipcs');
fprintf(fid,'%s \n','exit 0');
fclose(fid);
