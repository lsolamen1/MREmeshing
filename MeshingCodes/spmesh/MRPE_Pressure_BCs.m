%oooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooo
% Matt Mcgarry
% August 2013
%
% PURPOSE:
%   Define Pressure boundary conditions file for poroelastic reconstruction v11.
%   v11 does not take displacement BCs, it gets them from the measured displacement file
%   It only defines pressure BCs
%   Format:  BC# BCnode BCtype Re(BCval) Im(BCval)
%
% FILES:
%   Nodal Positions     (.nod)
%   Boundary Nodes      (.bnod)
%oooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooo

function MRPE_Pressure_BCs(nodf,btyp)
close all;

if (nargin==1)
    disp('Define BCs:')
    disp('Enter 1 for EA Phantom (type 2 on bottom face, type 1 everywhere else, [0 0 0 0 2 0 1])')
    disp('Enter 2 for IA Phantom (type 2 on bottom/top face, type 1 everywhere else, [0 0 0 0 2 2 1])')
    disp('Enter 3 for Brain (type 1 on top face/bottom face, type 2 everywhere else, [0 0 0 0 1 1 2]')
    disp('default: Type 1 everywhere [0 0 0 0 0 0 1]');    
    btyp=input('BC type, or values on faces [xmin xmax ymin ymax zmin zmax rest] >> ');
    if(isempty(btyp))
        btyp=[0 0 0 0 0 0 1];
    end
    btyp
    length(btyp)
    
    if(length(btyp)==1)
        if(btyp==1)
            btyp=[0 0 0 0 2 0 1];
        elseif(btyp==2)
            btyp=[0 0 0 0 2 2 1];
        elseif(btyp==3)
            btyp=[0 0 0 0 1 1 2];
        else
            error(['btyp = ' int2str(btyp) ' not recongnized'])  
        end
    elseif(length(btyp)==7)
        disp('Using supplied btyp array');
    else
        error('Unrecognized btyp array')        
    end
    
end

disp('Boundary conditions being applied to each face:')
disp(btyp)

junk = load(nodf);
xyz=junk(:,2:4);
junk=load([nodf(1:end-3) 'bnod']);
bnods=junk(:,2); 


nn = size(bnods,1);

xb=xyz(bnods,1);
yb=xyz(bnods,2);
zb=xyz(bnods,3);

rngx=range(xb);
rngy=range(yb);
rngz=range(zb);

tol=0.03; % tolerance for assigning nodes to faces/
face(1).nod=find(xb<(min(xb)+tol*rngx));
face(2).nod=find(xb>(max(xb)-tol*rngx));
face(3).nod=find(yb<(min(yb)+tol*rngy));
face(4).nod=find(yb>(max(yb)-tol*rngy));
face(5).nod=find(zb<(min(zb)+tol*rngz));
face(6).nod=find(zb>(max(zb)-tol*rngz));
face(7).nod=1:nn;


bcarray=zeros(nn,5);
bcarray(:,1)=1:nn;
bcarray(:,2)=bnods;

% Start with all nodes, the overwrite with xmin, xmax, ymin, ymax,zmin,zmax in that order.
for ii=[7 1 2 3 4 5 6];
    if(btyp(ii)~=0)
        for jj=1:length(face(ii).nod)
            bcarray(face(ii).nod(jj),3)=btyp(ii);
        end
    end
end

%*************************************************************************
% record BC type to make it clear which 
bc_name = [nodf(1:end-4),'.inv.pressure.bcs'];
bc_typ = [nodf(1:end-4),'.inv.bcfaces'];
fid=fopen(bc_typ,'wt');
  fprintf(fid,'BC type on faces [minx maxx miny maxy minz maxz rest] \n');
  fprintf(fid,'%i %i %i %i %i %i %i',btyp);
fclose(fid);
  
fid=fopen(bc_name,'wt');
  fprintf(fid,'%10i %10i %3i %.7e %.7e \n',bcarray')
fclose(fid);


% CALCULATE: Axis Limits
AXIS = ([min(xb)-tol*rngx max(xb)+tol*rngx  min(yb)-tol*rngy max(yb)+tol*rngy min(zb)-tol*rngz max(zb)+tol*rngz]);
figure
% Start with type 1 BCs
I1=find(bcarray(:,3)==1);
plot3(xb(I1),yb(I1),zb(I1),'b.');
hold on
I2=find(bcarray(:,3)==2);
plot3(xb(I2),yb(I2),zb(I2)','g.');
I0=find(bcarray(:,3)==0);

if(~isempty(I0)>0)
    warning('Some boundary nodes have not been assinged a BC type')
    plot3(xb(I0),yb(I0),zb(I0),'r.');
    legend('Pressure BC Type 1','Pressure BC Type 2','No pressure BC','Location',[0.15 0.025 0.1 0.1]) 
else
    legend('Pressure BC Type 1','Pressure BC Type 2','Location',[0.15 0.025 0.1 0.1]) 
end
xlabel('x');ylabel('y');zlabel('z');   
axis equal, axis(AXIS);


movefile(bc_name,'MRPE/')
movefile(bc_typ,'MRPE/')
