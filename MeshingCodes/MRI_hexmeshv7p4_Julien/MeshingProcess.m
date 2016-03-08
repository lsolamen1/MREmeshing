function [temps,node,deplacement]=MeshingProcess(interpo,deplacement,coord)
tmesh=tic;
%Start the meshing process:
dim=size(interpo.maskint);

% Assign node numbers to each interpolated voxel
nodnum=1:numel(deplacement.Urint(:,:,:,1));
nodnum=reshape(nodnum,size(deplacement.Urint(:,:,:,1)));

% Start in bottom corner, march elements along mesh, inserting them if they
% are inside the mask
intmp=zeros(prod(floor(dim/2)),27);
node.nel=0;
node.nodin=false(size(nodnum));
for ii=1:2:dim(1)-2
    for jj=1:2:dim(2)-2
        for kk=1:2:dim(3)-2
            if(sum(sum(sum(interpo.maskint(ii:ii+2,jj:jj+2,kk:kk+2))))==27) % Element is inside mask
                node.nel=node.nel+1;
                [intmp(node.nel,:)]=Hex27incidencelist(nodnum(ii:ii+2,jj:jj+2,kk:kk+2)); % Add element to mesh
                node.nodin(ii:ii+2,jj:jj+2,kk:kk+2)=true; % Tag these nodes as included
            end
        end
    end
end

intmp=intmp(1:node.nel,:);
temps.tin=toc(tmesh);
% Renumber Nodes (exclude unused nodes)
node.nod=zeros(prod(dim),3);
deplacement.uvwr=zeros(prod(dim),3);
deplacement.uvwi=zeros(prod(dim),3);
node.magimnod=zeros(prod(dim),1);
node.regnod=zeros(prod(dim),1);
node.idx=zeros(prod(dim),3);
old2newnod=zeros(size(nodnum));
%in=zeros(size(intmp));
node.nn=0;
node.nbnod=0;
node.bnod=zeros(prod(dim),1);

for ii=1:dim(1)
    for jj=1:dim(2)
        for kk=1:dim(3)
            if(node.nodin(ii,jj,kk))
                node.nn=node.nn+1;
                old2newnod(ii,jj,kk)=node.nn;
                %in(intmp==nodnum(ii,jj,kk))=nn;
                node.nod(node.nn,:)=[coord.xout(ii) coord.yout(jj) coord.zout(kk)];
                deplacement.uvwr(node.nn,:)=[deplacement.Urint(ii,jj,kk,1) deplacement.Urint(ii,jj,kk,2) deplacement.Urint(ii,jj,kk,3)];
                deplacement.uvwi(node.nn,:)=[deplacement.Uiint(ii,jj,kk,1) deplacement.Uiint(ii,jj,kk,2) deplacement.Uiint(ii,jj,kk,3)];
                node.idx(node.nn,:)=[ii jj kk];
                node.magimnod(node.nn)=interpo.MagIm_int(ii,jj,kk);
                node.regnod(node.nn)=interpo.regs_int(ii,jj,kk);
            end
        end
    end
end
node.in=old2newnod(intmp); % 330 times faster than in(intmp==nodnum(ii,jj,kk))=nn; inserted in nested loop above. Checked output files are the same with 'diff' command.
clear intmp

% Resize arrays to correct size.
node.nod=node.nod(1:node.nn,:);
deplacement.uvwr=deplacement.uvwr(1:node.nn,:);
deplacement.uvwi=deplacement.uvwi(1:node.nn,:);
node.magimnod=node.magimnod(1:node.nn);
node.regnod=node.regnod(1:node.nn);
node.idx=node.idx(1:node.nn,:);
node.bnod=node.bnod(1:node.nbnod);

disp('Building Element Connetivity')
% Build element connectivity number to determine boundary nodes
elcon=zeros(node.nn,1);
for ii=1:node.nel
    for jj=1:27
        elcon(node.in(ii,jj))=elcon(node.in(ii,jj))+1;
    end
end
disp('Finding Boundary nodes')
disp(['nn = ' int2str(node.nn)])
% Build Boundary node array based on nodal connectivity % THis step accounts for ~99% of the meshing time
bnum=1;


%Bnode finding - 2800x faster than old version below - bnodes come out
%in differet order but get all the same nodes.
bnodchecked=false([node.nn 1]);
for ii=1:node.nel
    for jj=1:27
        if(~bnodchecked(node.in(ii,jj)))
            if jj<=4 || (jj>=10 && jj<=13) % corner nodes =1,2,3,4,10,11,12,13
                if elcon(node.in(ii,jj))<8 %corner node is boundary node
                    node.bnod(bnum)=node.in(ii,jj);
                    bnum=bnum+1;
                end
            elseif (jj>=5 && jj<=8) || (jj>=14 && jj<=17) || (jj>=19 && jj<=22) % mid-edge nodes = 5,6,7,8,14,15,16,17,19,20,21,22
                if elcon(node.in(ii,jj))<4 %mid-edge node is boundary node
                    node.bnod(bnum)=node.in(ii,jj);
                    bnum=bnum+1;
                end
            elseif  jj==9 || jj==18 || (jj>=23 && jj<=26) %mid-fanodce nodes are 9,18,23,24,25,26
                if elcon(node.in(ii,jj))<2 %mid-face node is boundary node
                    node.bnod(bnum)=node.in(ii,jj);
                    bnum=bnum+1;
                end
            end
            bnodchecked(node.in(ii,jj))=true;
        end
    end
end

node.nbnod=bnum-1;
%tbnod=toc(tbni);
temps.tfemesh=toc(tmesh); % Time for full meshing process

end