function softprior(fname,flag,nremf,nremb,spfile)

load MRE_3DMotionData.mat
%load Mask.mat
load HeaderData.mat

load(spfile)

[nx ny nsl ndir]=size(MagIm);
if(nargin<2)
    flag=input('Do you need to segment regions? (yes = 1, already done = 0, no regions needed = 2)>> ');
end

if(flag==1)
    %regs=mask; %This is the actual mask
    %SP_segment;
    load(spfile)
    seg=zeros(1,1);
    seg=regs;%This is the variable from the mask file that defines the various regions
    
    
    %load segmentation.mat
    regs=double(regs);
    %regs(seg==1)=2;
    %regs(seg==2)=3;
    regs(seg==0)=0;
    regs(seg==1)=1;
    regs(seg==2)=2; 
    
    mask=regs;

    !save regs.mat regs
    
    save regs.mat mask
    disp(sprintf('Soft Prior region added'));
elseif(flag==2)
    regs=mask;
    regs(:)=0;
    regs=double(regs);
    save regs.mat regs
    disp(sprintf('No SP regions defined'));
end

load regs.mat
flist=dir('filelist*');
if (size(flist,1)==1) 
    nod=load(['~/meshes/', fname, '.mesh/' fname '.nod']);
else
    nod=load([fname, '.nod']);
end

dx=DirIndex(4,1)*1e-3;
dy=DirIndex(4,2)*1e-3;

% dx=spacing(1)*1e-3;
% dy=spacing(2)*1e-3;

MINz=min(nod(:,4));
MAXz=max(nod(:,4));

% manually remove slices
v=nremf+1;
vv=nremb;
vs(1)=v; vs(2)=vv+1; vm=max(vs);
reg=regs(:,:,v:end-vv);
dz=(MAXz-MINz)./(nsl-vm);        %number of removed slices +1

if(flag==1 || flag==0)
    [xx,yy,zz]=meshgrid(dx:dx:nx*dx,dy:dy:ny*dy,MINz:dz:MAXz);
    regions=interp3(xx,yy,zz,reg,nod(:,2),nod(:,3),nod(:,4),'nearest',1);
    regions(regions==0)=1;
elseif(flag==2)
    regions=zeros(length(nod),1);
end

fid=fopen('regions.dat','w');
for ii=1:length(regions)
    fprintf(fid,'%i %i \n',ii,regions(ii));
end 
fclose(fid);

load regions.dat
x1=find(regions(:,2)==2);
x2=find(regions(:,2)==3);


end
