function MRPE_BCs_sample_FWD(name,displacement)

%Load in .bel, .nod, and .bnod  Files:
bel_name = ['~/meshes/',name,'.mesh/',name,'.bel'];
bel = load(bel_name);
nod_name = ['~/meshes/',name,'.mesh/',name,'.nod'];
nod = load(nod_name);
bnod_name = ['~/meshes/',name,'.mesh/',name,'.bnod'];
bnods = load(bnod_name);

nn = size(bnods,1);
tn = size(nod,1);

% Distinguish Types of BC's:
outlist = ones(nn,1);

% SURFACES
ind = bnods(:,2);
bnod(:,1) = ind;
bnod(:,2:4) = nod(ind,2:4);

% all surfaces
s1 = bnod(:,2:4);

% bottom
ind = find(bnod(:,4)<=(1.01*min(bnod(:,4))));
s2 = bnod(ind,2:4); outlist(ind) = 2;

% top
ind = find(bnod(:,4)>=(max(bnod(:,4))-0.0001));
s3 = bnod(ind,2:4); outlist(ind) = 3;

% % corner 1
% ind = find(bnod(:,2)<=(min(bnod(:,2))+0.0001) & bnod(:,3)<=(min(bnod(:,3))+0.0001) & bnod(:,4)<=(min(bnod(:,4))+0.0001));
% s4 = bnod(ind,2:4); outlist(ind) = 4;
% 
% % corner 2
% ind = find(bnod(:,2)>=(max(bnod(:,2))-0.0001) & bnod(:,3)<=(min(bnod(:,3))+0.0001) & bnod(:,4)<=(min(bnod(:,4))+0.0001));
% s5 = bnod(ind,2:4); outlist(ind) = 5;
% 
% % corner 3
% ind = find(bnod(:,2)<=(min(bnod(:,2))+0.0001) & bnod(:,3)>=(max(bnod(:,3))-0.0001) & bnod(:,4)<=(min(bnod(:,4))+0.0001));
% s6 = bnod(ind,2:4); outlist(ind) = 6;

%*************************************************************************
bc_name = [name,'.bcs'];
fid=fopen(bc_name,'wt');
k = 0;
for ii = 1:nn
    n = bnod(ii,1);
    if outlist(ii) == 1 % ALL SURFACES
        k = k+1;
        bc=[2.0 0.0 0.0 2.0 0.0 0.0 2.0 0.0 0.0 1.0 0.0 0.0];
        fprintf(fid,'%d %d %d %.1f %.1f %.7e  %.1f %.1f %.7e  %.1f %.1f %.7e  %.1f %.1f %e\n',k,n,0,bc);    
    elseif outlist(ii) == 2 % BOTTOM SURFACE
        k = k+1;
        bc=[1.0 0.0 0.0 1.0 0.0 0.0 1.0 0.0 0.0 2.0 0.0 0.0];
        fprintf(fid,'%d %d %d %.1f %.1f %.7e  %.1f %.1f %.7e  %.1f %.1f %.7e  %.1f %.1f %e\n',k,n,0,bc);
    elseif outlist(ii) == 3 % TOP SURFACE
        k = k+1;
        bc=[1.0 0.0 0.0 1.0 0.0 0.0 1.0 0.0 displacement 2.0 0.0 0.0];
        fprintf(fid,'%d %d %d %.1f %.1f %.7e  %.1f %.1f %.7e  %.1f %.1f %.7e  %.1f %.1f %e\n',k,n,0,bc);
    elseif outlist(ii) == 4 % CORNER NODE 1
        k = k+1;
        bc=[1.0 0.0 0.0 1.0 0.0 0.0 1.0 0.0 0.0 2.0 0.0 0.0];
        fprintf(fid,'%d %d %d %.1f %.1f %.7e  %.1f %.1f %.7e  %.1f %.1f %.7e  %.1f %.1f %e\n',k,n,0,bc);
    elseif outlist(ii) == 5 % CORNER NODE 2
        k = k+1;
        bc=[2.0 0.0 0.0 1.0 0.0 0.0 1.0 0.0 0.0 2.0 0.0 0.0];
        fprintf(fid,'%d %d %d %.1f %.1f %.7e  %.1f %.1f %.7e  %.1f %.1f %.7e  %.1f %.1f %e\n',k,n,0,bc);
    elseif outlist(ii) == 6 % CORNER NODE 3
        k = k+1;
        bc=[1.0 0.0 0.0 2.0 0.0 0.0 1.0 0.0 0.0 2.0 0.0 0.0];
        fprintf(fid,'%d %d %d %.1f %.1f %.7e  %.1f %.1f %.7e  %.1f %.1f %.7e  %.1f %.1f %e\n',k,n,0,bc);       
    end
end
fclose(fid);

%CALCULATE: Axis Limits
AXIS = ([0 max(nod(:,2)) 0 max(nod(:,3)) 0 max(nod(:,4))]);

hold on;
plot3(s1(:,1),s1(:,2),s1(:,3),'b.')     % ALL
plot3(s2(:,1),s2(:,2),s2(:,3),'k.')     % BOTTOM
plot3(s3(:,1),s3(:,2),s3(:,3),'r.')     % TOP
% plot3(s4(:,1),s4(:,2),s4(:,3),'y.')     % CN1
% plot3(s5(:,1),s5(:,2),s5(:,3),'y.')     % CN2
% plot3(s6(:,1),s6(:,2),s6(:,3),'y.')     % CN3
hold off;

title('\bfBoundary Nodes')
xlabel('\bfX');
ylabel('\bfY');
zlabel('\bfZ');
axis equal, axis(AXIS);
grid;

% calculate node numbers for top face of sample
tpf_name = [name,'.tpf'];
fid=fopen(tpf_name,'wt');
k=0;

list=zeros(tn,1);
ind = find(nod(:,4)>=(max(nod(:,4))-0.0001));
list(ind) = 3;

for ii=1:tn
    if list(ii) == 3 % TOP SURFACE
        k = k+1;
        fprintf(fid,'%i %d\n',k,ii);
    end
end
fclose(fid);