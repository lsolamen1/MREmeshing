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

function MRPE_BCs_FWD(name)

%Load in .bel, .nod, and .bnod  Files:
bel_name = ['~/meshes/',name,'.mesh/',name,'.bel'];
bel = load(bel_name);
nod_name = ['~/meshes/',name,'.mesh/',name,'.nod'];
nod = load(nod_name);
bnod_name = ['~/meshes/',name,'.mesh/',name,'.bnod'];
bnods = load(bnod_name);

% OTHER STUFF
% bel_name = [name,'.bel'];
% bel = load(bel_name);
% nod_name = [name,'.nod'];
% nod = load(nod_name);
% bnod_name = [name,'.bnod'];
% bnods = load(bnod_name);

% UC_name = [name,'.v3c'];
% UC = load(UC_name);
% Ux = UC(:,2)+1i.*UC(:,3);
% Uy = UC(:,4)+1i.*UC(:,5);
% Uz = UC(:,6)+1i.*UC(:,7);


nn = size(bnods,1);

% Distinguish Types of BC's:
outlist = ones(nn,1);

% SURFACES
ind = bnods(:,2);
bnod(:,1) = ind;
bnod(:,2:4) = nod(ind,2:4);

% all surfaces
s1 = bnod(:,2:4);

% bottom
ind = find(bnod(:,4)<=1.1*min(bnod(:,4)));
s2 = bnod(ind,2:4); outlist(ind) = 2;

% % top
% ind = find(bnod(:,4)>=(max(bnod(:,4))-0.0001));
% s3 = bnod(ind,2:4); outlist(ind) = 3;

% % sides
% ind = find(bnod(:,2)>=(min(bnod(:,2))+0.0001));
% s4 = bnod(ind,2:4); outlist(ind) = 3;

% simulation for column
% bottom
% ind = find(bnod(:,4)==min(bnod(:,4)));
% s2 = bnod(ind,2:4); outlist(ind) = 2;
% 
% % top
% ind = find(bnod(:,4)==max(bnod(:,4)));
% s3 = bnod(ind,2:4); outlist(ind) = 3;
% 
% % top edges
% ind = find((bnod(:,2)==min(bnod(:,2)) & bnod(:,4)==max(bnod(:,4))) | ...
%            (bnod(:,3)==min(bnod(:,3)) & bnod(:,4)==max(bnod(:,4))) | ...
%            (bnod(:,2)==max(bnod(:,2)) & bnod(:,4)==max(bnod(:,4))) | ...
%            (bnod(:,3)==max(bnod(:,3)) & bnod(:,4)==max(bnod(:,4))));
% s4 = bnod(ind,2:4); outlist(ind) = 4;    
% 
% % bottom edges
% ind = find((bnod(:,2)==min(bnod(:,2)) & bnod(:,4)==min(bnod(:,4))) | ...
%            (bnod(:,3)==min(bnod(:,3)) & bnod(:,4)==min(bnod(:,4))) | ...
%            (bnod(:,2)==max(bnod(:,2)) & bnod(:,4)==min(bnod(:,4))) | ...
%            (bnod(:,3)==max(bnod(:,3)) & bnod(:,4)==min(bnod(:,4))));
% s5 = bnod(ind,2:4); outlist(ind) = 5;
% 
% % side edges
% ind = find((bnod(:,2)==min(bnod(:,2)) & bnod(:,3)==min(bnod(:,3))) | ...
%            (bnod(:,2)==min(bnod(:,2)) & bnod(:,3)==max(bnod(:,3))) | ...
%            (bnod(:,2)==max(bnod(:,2)) & bnod(:,3)==min(bnod(:,3))) | ...
%            (bnod(:,2)==max(bnod(:,2)) & bnod(:,3)==max(bnod(:,3))));
% s6 = bnod(ind,2:4); outlist(ind) = 6;

%*************************************************************************
bc_name = [name,'.bcs'];
fid=fopen(bc_name,'wt');
%fprintf(fid,['1 ',num2str(nn),'\n']);
k = 0;
A = -100;
P = 0.0;
for ii = 1:nn
    n = bnod(ii,1);
    if outlist(ii) == 1 % ALL SURFACES
        k = k+1;
        bc=[2.0 0.0 0.0 2.0 0.0 0.0 2.0 0.0 0.0 1.0 0.0 0.0];
        fprintf(fid,'%d %d %d %.1f %.1f %.7e  %.1f %.1f %.7e  %.1f %.1f %.7e  %.1f %.1f %e\n',k,n,0,bc);

%         bc=[1.0 0.0 real(Ux(n)) imag(Ux(n)) 1.0 0.0 real(Uy(n)) imag(Uy(n)) 1.0 0.0 real(Uz(n)) imag(Uz(n)) 1.0 0.0 0.0];
%         bc=[2.0 0.0 real(0.0) imag(0.0) 2.0 0.0 real(0.0) imag(0.0) 2.0 0.0 real(0.0) imag(0.0) 1.0 0.0 0.0];        
%         fprintf(fid,'%d %d %d %.1f %.1f (%.7e,%.7e)  %.1f %.1f (%.7e,%.7e)  %.1f %.1f (%.7e,%.7e)  %.1f %.1f %.e\n',k,n,0,bc);          
    elseif outlist(ii) == 2 % BOTTOM SURFACE
        k = k+1;
        bc=[1.0 0.0 1e-6 1.0 0.0 0.0 1.0 0.0 0.0 2.0 0.0 0.0];
        fprintf(fid,'%d %d %d %.1f %.1f %.7e  %.1f %.1f %.7e  %.1f %.1f %.7e  %.1f %.1f %e\n',k,n,0,bc);
%         bc=[1.0 0.0 real(Ux(n)) imag(Ux(n)) 1.0 0.0 real(Uy(n)) imag(Uy(n)) 1.0 0.0 real(Uz(n)) imag(Uz(n)) 2.0 0.0 0.0];
%         fprintf(fid,'%d %d %d %.1f %.1f (%.7e,%.7e)  %.1f %.1f (%.7e,%.7e)  %.1f %.1f (%.7e,%.7e)  %.1f %.1f %.e\n',k,n,0,bc); 
%     elseif outlist(ii) == 2 % BOTTOM SURFACE
%         k = k+1;
%         bc=[1.0 0.0 0.0 1.0 0.0 0.0 1.0 0.0 1e-4 2.0 0.0 0.0];
%         fprintf(fid,'%d %d %d %.1f %.1f %.7e  %.1f %.1f %.7e  %.1f %.1f %.7e  %.1f %.1f %e\n',k,n,0,bc);        
%     elseif outlist(ii) == 3 % TOP SURFACE
%         k = k+1;
%         bc=[1.0 0.0 0.0 1.0 0.0 0.0 1.0 0.0 0.0 2.0 0.0 0.0];
%         fprintf(fid,'%d %d %d %.1f %.1f %.7e  %.1f %.1f %.7e  %.1f %.1f %.7e  %.1f %.1f %e\n',k,n,0,bc); 
%     elseif outlist(ii) == 4 % TOP EDGES
%         k = k+1;
%         bc=[1.0 0.0 0.0 1.0 0.0 0.0 2.0 0.0 A/2 1.0 0.0 P];
%         fprintf(fid,'%d %d %d %.1f %.1f %.7e  %.1f %.1f %.7e  %.1f %.1f %.7e  %.1f %.1f %e\n',k,n,0,bc);
%     elseif outlist(ii) == 5 % BOTTOM EDGES
%         k = k+1;
%         bc=[1.0 0.0 0.0 1.0 0.0 0.0 1.0 0.0 0.0 2.0 0.0 0.0];
%         fprintf(fid,'%d %d %d %.1f %.1f %.7e  %.1f %.1f %.7e  %.1f %.1f %.7e  %.1f %.1f %e\n',k,n,0,bc);
%     elseif outlist(ii) == 6 % SIDE EDGES
%         k = k+1;
%         bc=[1.0 0.0 0.0 1.0 0.0 0.0 2.0 0.0 0.0 2.0 0.0 0.0];
%         fprintf(fid,'%d %d %d %.1f %.1f %.7e  %.1f %.1f %.7e  %.1f %.1f %.7e  %.1f %.1f %e\n',k,n,0,bc);
    end
end
fclose(fid);

% %LOAD material file
% %mat_name = ['cyl.mat'];
% mat_name = [name,'.mat'];
% mat = load(mat_name,'-ASCII');
% 
% %LOAD Pressure Source Node File:
% %psn_name = ['cyl.psn'];
% psn_name = [name,'.psn'];
% psn = load(psn_name);
% 
% mtr = ones(nn,4);
% mtr(:,1) = 1:nn;
% mtr(:,4) = mat(1,5);
% for i = 1:length(psn)
%     mtr(psn(i,2),4) = mat(psn(i,3),5);
% end
% 
% %WRITE .mtr File:
% mtr_name = [name,'.inv.mtr'];
% fid=fopen(mtr_name,'wt');
% fprintf(fid,'%d %d %d %e\n',mtr');
% fclose(fid);

%CALCULATE: Axis Limits
AXIS = ([0 max(nod(:,2)) 0 max(nod(:,3)) 0 max(nod(:,4))]);

hold on;
plot3(s1(:,1),s1(:,2),s1(:,3),'b.')     % ALL
plot3(s2(:,1),s2(:,2),s2(:,3),'k.')     % BOTTOM
% plot3(s3(:,1),s3(:,2),s3(:,3),'r.')     % TOP
% plot3(s4(:,1),s4(:,2),s4(:,3),'y.')     % TOP EDGES
hold off;

title('\bfBoundary Nodes')
xlabel('\bfX');
ylabel('\bfY');
zlabel('\bfZ');
axis equal, axis(AXIS);
grid;
