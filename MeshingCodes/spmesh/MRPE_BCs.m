%oooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooo
% Phillip R. Perriñez, Ph.D.
% October 2008
%
% PURPOSE:
%   Define boundary conditions file for pertaining to the
%   poroelasticity based elastic parameter reconstruction.
%
% FILES:
%   Nodal Positions     (.nod)
%   Boundary Nodes      (.bnod)
%   Boundary Elements   (.bel)
%   Displacement        (filename.v3c)
%oooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooo
function MRPE_BCs(name,SIM,type)


if (nargin==1)
    SIM=input('Is this simulated or real data? (sim=1, real=2) >> ');
    type=input('Is this phantom or brain data? (phan=1, brain=2, IAphan=3, all typ 1=4, btmx typ2=5) >> ');
end

if SIM == 1
    nod = load(['~/meshes/',name,'.mesh/',name,'.nod']);
    bnods = load(['~/meshes/',name,'.mesh/',name,'.bnod']);
else
    nod = load(['./',name,'.nod']);
    bnods = load(['./',name,'.bnod']);
end

% Load displace.v3c
if ~isempty(dir('*.v3c*'));
    UC_name = [name,'.v3c'];
elseif ~isempty(dir('*.dsp*'));
    UC_name = [name,'.dsp'];
else
    error('Cannot Find Displacement File');
end

UC = load(UC_name);
Ux = UC(:,2)+1i.*UC(:,3);
Uy = UC(:,4)+1i.*UC(:,5);
Uz = UC(:,6)+1i.*UC(:,7);

nn = size(bnods,1);

% Distinguish Types of BC's:
outlist = ones(nn,1);

% SURFACES
ind = bnods(:,2);
bnod(:,1) = ind;
bnod(:,2:4) = nod(ind,2:4);

% ALL SURFACES
s1 = bnod(:,2:4);

% % BOTTOM
ind = find(bnod(:,4)<=(min(bnod(:,4)))+0.01*range(bnod(:,4))); %For bgmesh to define the bottom slice pressure boundary condition.
%ind = find(bnod(:,4)<=1.1*min(bnod(:,4));
s2 = bnod(ind,2:4); outlist(ind) = 2;
    
% TOP
ind = find(bnod(:,4)>=0.99*max(bnod(:,4)));
s3 = bnod(ind,2:4); outlist(ind) = 3;

% min x
ind = find(bnod(:,2)<=1.01*min(bnod(:,2)));
s4 = bnod(ind,2:4); outlist(ind) = 4;

%*************************************************************************
bc_name = [name,'.inv.bcs'];
fid=fopen(bc_name,'wt');
k = 0; kk=0; jj=0;
for ii = 1:nn
    n = bnod(ii,1);
    if type==1  % phantom
        if ((outlist(ii) == 1)||(outlist(ii) == 4)) % All SURFACES
            k = k+1;
            bc=[1.0 0.0 real(Ux(n)) imag(Ux(n)) 1.0 0.0 real(Uy(n)) imag(Uy(n)) 1.0 0.0 real(Uz(n)) imag(Uz(n)) 1.0 0.0 0.0];
            fprintf(fid,'%d %d  %.1f %e (%.7e,%.7e)  %.1f %e (%.7e,%.7e)  %.1f %e (%.7e,%.7e)  %.1f %e %e\n',k,n,bc);
        elseif outlist(ii) == 2 % BOTTOM SURFACE
            k = k+1;
            bc=[1.0 0.0 real(Ux(n)) imag(Ux(n)) 1.0 0.0 real(Uy(n)) imag(Uy(n)) 1.0 0.0 real(Uz(n)) imag(Uz(n)) 2.0 0.0 0.0];
            fprintf(fid,'%d %d  %.1f %e (%.7e,%.7e)  %.1f %e (%.7e,%.7e)  %.1f %e (%.7e,%.7e)  %.1f %e %e\n',k,n,bc);      
        elseif outlist(ii) == 3 % TOP SURFACE
            k = k+1;
            bc=[1.0 0.0 real(Ux(n)) imag(Ux(n)) 1.0 0.0 real(Uy(n)) imag(Uy(n)) 1.0 0.0 real(Uz(n)) imag(Uz(n)) 1.0 0.0 0.0];
            fprintf(fid,'%d %d  %.1f %e (%.7e,%.7e)  %.1f %e (%.7e,%.7e)  %.1f %e (%.7e,%.7e)  %.1f %e %e\n',k,n,bc);     
        end
        
    elseif type==2  % brain
        if ((outlist(ii) == 1)||(outlist(ii) == 4)) % All SURFACES
            k = k+1;
            bc=[1.0 0.0 real(Ux(n)) imag(Ux(n)) 1.0 0.0 real(Uy(n)) imag(Uy(n)) 1.0 0.0 real(Uz(n)) imag(Uz(n)) 2.0 0.0 0.0];
            fprintf(fid,'%d %d  %.1f %e (%.7e,%.7e)  %.1f %e (%.7e,%.7e)  %.1f %e (%.7e,%.7e)  %.1f %e %e\n',k,n,bc);
        elseif outlist(ii) == 2 % BOTTOM SURFACE
            k = k+1;
            bc=[1.0 0.0 real(Ux(n)) imag(Ux(n)) 1.0 0.0 real(Uy(n)) imag(Uy(n)) 1.0 0.0 real(Uz(n)) imag(Uz(n)) 1.0 0.0 0.0];
            fprintf(fid,'%d %d  %.1f %e (%.7e,%.7e)  %.1f %e (%.7e,%.7e)  %.1f %e (%.7e,%.7e)  %.1f %e %e\n',k,n,bc);
        elseif outlist(ii) == 3 % TOP SURFACE
            k = k+1;
            bc=[1.0 0.0 real(Ux(n)) imag(Ux(n)) 1.0 0.0 real(Uy(n)) imag(Uy(n)) 1.0 0.0 real(Uz(n)) imag(Uz(n)) 1.0 0.0 0.0];
            fprintf(fid,'%d %d  %.1f %e (%.7e,%.7e)  %.1f %e (%.7e,%.7e)  %.1f %e (%.7e,%.7e)  %.1f %e %e\n',k,n,bc);
        end
        
    elseif type==3  % IA phantom
        if ((outlist(ii) == 1)||(outlist(ii) == 4)) % All SURFACES
            k = k+1;
            bc=[1.0 0.0 real(Ux(n)) imag(Ux(n)) 1.0 0.0 real(Uy(n)) imag(Uy(n)) 1.0 0.0 real(Uz(n)) imag(Uz(n)) 2.0 0.0 0.0];
            fprintf(fid,'%d %d  %.1f %e (%.7e,%.7e)  %.1f %e (%.7e,%.7e)  %.1f %e (%.7e,%.7e)  %.1f %e %e\n',k,n,bc);
        elseif outlist(ii) == 2 % BOTTOM SURFACE
            k = k+1;
            bc=[1.0 0.0 real(Ux(n)) imag(Ux(n)) 1.0 0.0 real(Uy(n)) imag(Uy(n)) 1.0 0.0 real(Uz(n)) imag(Uz(n)) 2.0 0.0 0.0];
            fprintf(fid,'%d %d  %.1f %e (%.7e,%.7e)  %.1f %e (%.7e,%.7e)  %.1f %e (%.7e,%.7e)  %.1f %e %e\n',k,n,bc);
        elseif outlist(ii) == 3 % TOP SURFACE
            k = k+1;
            bc=[1.0 0.0 real(Ux(n)) imag(Ux(n)) 1.0 0.0 real(Uy(n)) imag(Uy(n)) 1.0 0.0 real(Uz(n)) imag(Uz(n)) 1.0 0.0 0.0];
            fprintf(fid,'%d %d  %.1f %e (%.7e,%.7e)  %.1f %e (%.7e,%.7e)  %.1f %e (%.7e,%.7e)  %.1f %e %e\n',k,n,bc);
        end  
    elseif type==4  % All typ 1 pressure BCs
        if ((outlist(ii) == 1)||(outlist(ii) == 4)) % All SURFACES
            k = k+1;
            bc=[1.0 0.0 real(Ux(n)) imag(Ux(n)) 1.0 0.0 real(Uy(n)) imag(Uy(n)) 1.0 0.0 real(Uz(n)) imag(Uz(n)) 1.0 0.0 0.0];
            fprintf(fid,'%d %d  %.1f %e (%.7e,%.7e)  %.1f %e (%.7e,%.7e)  %.1f %e (%.7e,%.7e)  %.1f %e %e\n',k,n,bc);
        elseif outlist(ii) == 2 % BOTTOM SURFACE
            k = k+1;
            bc=[1.0 0.0 real(Ux(n)) imag(Ux(n)) 1.0 0.0 real(Uy(n)) imag(Uy(n)) 1.0 0.0 real(Uz(n)) imag(Uz(n)) 1.0 0.0 0.0];
            fprintf(fid,'%d %d  %.1f %e (%.7e,%.7e)  %.1f %e (%.7e,%.7e)  %.1f %e (%.7e,%.7e)  %.1f %e %e\n',k,n,bc);
        elseif outlist(ii) == 3 % TOP SURFACE
            k = k+1;
            bc=[1.0 0.0 real(Ux(n)) imag(Ux(n)) 1.0 0.0 real(Uy(n)) imag(Uy(n)) 1.0 0.0 real(Uz(n)) imag(Uz(n)) 1.0 0.0 0.0];
            fprintf(fid,'%d %d  %.1f %e (%.7e,%.7e)  %.1f %e (%.7e,%.7e)  %.1f %e (%.7e,%.7e)  %.1f %e %e\n',k,n,bc);
        end
    elseif type==5  % typ 2 pressure BCs on x=0 face
        if outlist(ii) == 1 % All SURFACES
            k = k+1;
            bc=[1.0 0.0 real(Ux(n)) imag(Ux(n)) 1.0 0.0 real(Uy(n)) imag(Uy(n)) 1.0 0.0 real(Uz(n)) imag(Uz(n)) 1.0 0.0 0.0];
            fprintf(fid,'%d %d  %.1f %e (%.7e,%.7e)  %.1f %e (%.7e,%.7e)  %.1f %e (%.7e,%.7e)  %.1f %e %e\n',k,n,bc);
        elseif outlist(ii) == 2 % BOTTOM SURFACE
            k = k+1;
            bc=[1.0 0.0 real(Ux(n)) imag(Ux(n)) 1.0 0.0 real(Uy(n)) imag(Uy(n)) 1.0 0.0 real(Uz(n)) imag(Uz(n)) 1.0 0.0 0.0];
            fprintf(fid,'%d %d  %.1f %e (%.7e,%.7e)  %.1f %e (%.7e,%.7e)  %.1f %e (%.7e,%.7e)  %.1f %e %e\n',k,n,bc);
        elseif outlist(ii) == 3 % TOP SURFACE
            k = k+1;
            bc=[1.0 0.0 real(Ux(n)) imag(Ux(n)) 1.0 0.0 real(Uy(n)) imag(Uy(n)) 1.0 0.0 real(Uz(n)) imag(Uz(n)) 1.0 0.0 0.0];
            fprintf(fid,'%d %d  %.1f %e (%.7e,%.7e)  %.1f %e (%.7e,%.7e)  %.1f %e (%.7e,%.7e)  %.1f %e %e\n',k,n,bc);
        elseif outlist(ii) == 4 % min x SURFACE
            k = k+1;
            bc=[1.0 0.0 real(Ux(n)) imag(Ux(n)) 1.0 0.0 real(Uy(n)) imag(Uy(n)) 1.0 0.0 real(Uz(n)) imag(Uz(n)) 2.0 0.0 0.0];
            fprintf(fid,'%d %d  %.1f %e (%.7e,%.7e)  %.1f %e (%.7e,%.7e)  %.1f %e (%.7e,%.7e)  %.1f %e %e\n',k,n,bc);

        end   
    else
        disp('What type of dataset is this?')
    end
    
    % finding type 1 and 2 bc nodes
    if outlist(ii) ~= 0 && bc(13) == 1
        kk=kk+1;        
        bc1(kk,:) = bnod(ii,2:4);
    elseif outlist(ii) ~= 0 && bc(13) == 2
        jj=jj+1;
        bc2(jj,:) = bnod(ii,2:4);
    end
end
fclose(fid);

% CALCULATE: Axis Limits
AXIS = ([0 max(nod(:,2)) 0 max(nod(:,3)) 0 max(nod(:,4))]);

hold on;
if type==1 && SIM==1
    figure
    plot3(bc1(:,1),bc1(:,2),bc1(:,3),'k.'); hold on;
    plot3(bc2(:,1),bc2(:,2),bc2(:,3),'r.')
    hold off;
    title('\bfBoundary Nodes for Simulated Phantom Data')
    xlabel('\bfX');
    ylabel('\bfY');
    zlabel('\bfZ');
    legend('Pressure BC Type 1','Pressure BC Type 2','Location',[0.15 0.025 0.1 0.1])    
    axis equal, axis(AXIS);
    grid;
elseif type==1 && SIM~=1
    figure
    plot3(bc1(:,1),bc1(:,2),bc1(:,3),'k.'); hold on;
    plot3(bc2(:,1),bc2(:,2),bc2(:,3),'r.')
    hold off;
    title('\bfBoundary Nodes for Real Phantom Data')
    xlabel('\bfX');
    ylabel('\bfY');
    zlabel('\bfZ');
    legend('Pressure BC Type 1','Pressure BC Type 2','Location',[0.15 0.025 0.1 0.1])    
    axis equal, axis(AXIS);
    grid;
elseif type==2 && SIM==1
    figure
    plot3(bc1(:,1),bc1(:,2),bc1(:,3),'k.'); hold on;
    plot3(bc2(:,1),bc2(:,2),bc2(:,3),'r.')
    hold off;
    title('\bfBoundary Nodes for Simulated Brain Data')
    xlabel('\bfX');
    ylabel('\bfY');
    zlabel('\bfZ');
    legend('Pressure BC Type 1','Pressure BC Type 2','Location',[0.15 0.025 0.1 0.1])    
    axis equal, axis(AXIS);
    grid;
elseif type==2 && SIM~=1
    plot3(bc1(:,1),bc1(:,2),bc1(:,3),'k.'); hold on;
    plot3(bc2(:,1),bc2(:,2),bc2(:,3),'r.')
    hold off;
    title('\bfBoundary Nodes for Real Brain Data')
    xlabel('\bfX');
    ylabel('\bfY');
    zlabel('\bfZ');
    legend('Pressure BC Type 1','Pressure BC Type 2','Location',[0.15 0.025 0.1 0.1])    
    axis equal, axis(AXIS);
    grid;
elseif type==3 && SIM~=1
    figure
    plot3(bc1(:,1),bc1(:,2),bc1(:,3),'k.'); hold on;
    plot3(bc2(:,1),bc2(:,2),bc2(:,3),'r.')
    hold off;
    title('\bfBoundary Nodes for Real IA Phantom Data')
    xlabel('\bfX');
    ylabel('\bfY');
    zlabel('\bfZ');
    legend('Pressure BC Type 1','Pressure BC Type 2','Location',[0.15 0.025 0.1 0.1])    
    axis equal, axis(AXIS);
    grid;
end

eval(['!mv ',name,'.inv.bcs MRPE/'])

%% Creating Pressure Boundary Condition FIle
pbcs=load(strcat(name,'.bnod'));

ptypein=input('Do you want Type 1 (1) or Type 2 (2) pressure boundary condition (Default=1)>>:');
if isempty(ptypein);
    ptypein=1;
end

pbcs(:,3)=ptypein; %Pressure Type
pbcs(:,4)=0; %Real Part of the pressure values are 0
pbcs(:,5)=0; %Imag Part of the pressure values are 0

fid=fopen([name,'.pbcs'],'wt');
fprintf(fid,'%7i %7i %3i %3i %3i \n',pbcs');
fclose(fid);

nod=load([name,'.nod']); 
nod(:,5)=1; 
fid=fopen([name,'.nod'],'wt');
fprintf(fid, '%7i %15.8e %15.8e %15.8e %15.8e \n',nod');
fclose(fid);     
end
