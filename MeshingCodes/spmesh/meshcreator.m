function meshcreator(fname,inclusion)

% Location of mesh creator
meshex='/home/apattiso/bin/3d_box_mesh_gen.x';

if inclusion==2         % 3-inclusion
    X0=0; Xend=0.1; X_range=[X0 Xend];
    Y0=0; Yend=0.05; Y_range=[Y0 Yend];
    Z0=0; Zend=0.05; Z_range=[Z0 Zend];
    nodalspacing=[2.5/1000 2.5/1000 2.5/1000];
elseif inclusion==1     % 1-inclusion
    X0=0; Xend=0.04; X_range=[X0 Xend];
    Y0=0; Yend=0.04; Y_range=[Y0 Yend];
    Z0=0; Zend=0.04; Z_range=[Z0 Zend];    
    nodalspacing=[2/1000 2/1000 2/1000];
elseif inclusion==3
    X0=0; Xend=0.1; X_range=[X0 Xend];
    Y0=0; Yend=0.05; Y_range=[Y0 Yend];
    Z0=0; Zend=0.05; Z_range=[Z0 Zend];
    nodalspacing=[2.5/1000 2.5/1000 2.5/1000];  
else
    X0=0; Xend=296*0.002; X_range=[X0 Xend];
    Y0=0; Yend=8*0.002; Y_range=[Y0 Yend];
    Z0=0; Zend=8*0.002; Z_range=[Z0 Zend];    
    nodalspacing=[0.002 0.002 0.002];
end

Xnod = Xend/nodalspacing(1);
Ynod = Yend/nodalspacing(2);
Znod = Zend/nodalspacing(3);

%% Create a mesh creation input file in the folder
runfile=['./3Dmesh_input.inp'];
fid=fopen(runfile,'w');
fprintf(fid,[num2str(X_range) '\n']);
fprintf(fid,[num2str(Y_range) '\n']);
fprintf(fid,[num2str(Z_range) '\n']);
fprintf(fid,['./' fname '.nod\n']);
fprintf(fid,['./' fname '.elm\n']);
fprintf(fid,[num2str(Xnod) ',' num2str(Ynod) ',' num2str(Znod) '\n']);
fprintf(fid,['1' '\n']);

% get node and element files
eval(['!' meshex ' < ' runfile])

% get boundary element file
getBelFromTet2(fname);
BEL = load([fname '.bel']);
SZ = size(BEL);
BELF = zeros(SZ(1),6);
BELF(:,1:4) = BEL;
BELF(:,5) = 1;
save([fname '.bel'], 'BELF','-ascii')

%% 1-inclusion phantom
if inclusion == 1
    nod = load([fname '.nod']);
    Nsz = size(nod);
    radius = 0.01;      % meters
    cdist = 1;
    k=0;
    
    % find center-most node
    for ii=1:Nsz(1)
        dist = sqrt((nod(ii,2)-(Xend/2.)).^2 + (nod(ii,3)-(Yend/2.)).^2 + (nod(ii,4)-(Zend/2.)).^2);
        if dist <= cdist
            cdist = dist;
            cnod = ii;
        end
    end
    
    % find nodes in inclusion
    for ii=1:Nsz(1)
        dist = sqrt((nod(ii,2)-nod(cnod,2)).^2 + (nod(ii,3)-nod(cnod,3)).^2 + (nod(ii,4)-nod(cnod,4)).^2);
        if dist <= radius
            k = k+1;
            incl_nod(k,1) = ii;
            incl_nod(k,2) = 2;
        end
    end
    psnfile(fname,incl_nod)

%% 3-inclusion phantom
elseif inclusion == 2
    nod = load([fname '.nod']);
    nn = size(nod);
    radius = 0.05;

    % find inclusion-1 center node
    cdist = 1;
    k=0;
    for ii=1:nn(1)
        dist = sqrt((nod(ii,2)-(Xend/4.-radius/2)).^2 + (nod(ii,3)-(Yend/2.)).^2 + (nod(ii,4)-(Zend/2.)).^2);
        if dist <= cdist
            cdist = dist;
            cnod = ii;
        end
    end
    
    % find nodes in inclusion
    for ii=1:nn(1)
        dist = sqrt((nod(ii,2)-nod(cnod,2)).^2 + (nod(ii,3)-nod(cnod,3)).^2 + (nod(ii,4)-nod(cnod,4)).^2);
        if dist <= radius
            k = k+1;
            incl_nod(k,1) = ii;
            incl_nod(k,2) = 2;
        end
    end

    % find inclusion-2 center node
    cdist = 1;
    for ii=1:nn(1)
        dist = sqrt((nod(ii,2)-(Xend/2.)).^2 + (nod(ii,3)-(Yend/2.)).^2 + (nod(ii,4)-(Zend/2.)).^2);
        if dist <= cdist
            cdist = dist;
            cnod = ii;
        end
    end
    
    % find nodes in inclusion
    for ii=1:nn(1)
        dist = sqrt((nod(ii,2)-nod(cnod,2)).^2 + (nod(ii,3)-nod(cnod,3)).^2 + (nod(ii,4)-nod(cnod,4)).^2);
        if dist <= radius
            k = k+1;
            incl_nod(k,1) = ii;
            incl_nod(k,2) = 3;
        end
    end    
    
    % find inclusion-3 center node
    cdist = 1;
    for ii=1:nn(1)
        dist = sqrt((nod(ii,2)-(3*Xend/4.+radius/2)).^2 + (nod(ii,3)-(Yend/2.)).^2 + (nod(ii,4)-(Zend/2.)).^2);
        if dist <= cdist
            cdist = dist;
            cnod = ii;
        end
    end
    
    % find nodes in inclusion
    for ii=1:nn(1)
        dist = sqrt((nod(ii,2)-nod(cnod,2)).^2 + (nod(ii,3)-nod(cnod,3)).^2 + (nod(ii,4)-nod(cnod,4)).^2);
        if dist <= radius
            k = k+1;
            incl_nod(k,1) = ii;
            incl_nod(k,2) = 4;
        end
    end
    psnfile(fname,incl_nod)
    %% conical inclusions
elseif inclusion==3
    nod = load([fname '.nod']);
    nn = size(nod);
    radius = 0.01;
    lth=0.03;

    % find axes nodes
    cdist = 1;
    k=0;
    
    ind=find(nod(:,3)==Yend/2 & nod(:,4)==Zend/2 & nod(:,2)>=0.01 & nod(:,2)<=0.04);
    for ii=1:length(ind)
        radius = tan(0.463647609000806)*(0.04-nod(ind(ii),2));
        incl_nods = find((sqrt((nod(:,3)-nod(ind(ii),3)).^2+(nod(:,4)-nod(ind(ii),4)).^2) <= radius) & (nod(:,2)==nod(ind(ii),2)));
        for jj=1:length(incl_nods)
            k=k+1;
            incl_nod(k,1)=incl_nods(jj);
            incl_nod(k,2)=2;
        end
    end
        
    ind=find(nod(:,3)==Yend/2 & nod(:,4)==Zend/2 & nod(:,2)>=0.06 & nod(:,2)<=0.09);
    for ii=1:length(ind)
        radius = tan(0.463647609000806)*(nod(ind(ii),2)-0.06);
        incl_nods = find((sqrt((nod(:,3)-nod(ind(ii),3)).^2+(nod(:,4)-nod(ind(ii),4)).^2) <= radius) & (nod(:,2)==nod(ind(ii),2)));
        for jj=1:length(incl_nods)
            k=k+1;
            incl_nod(k,1)=incl_nods(jj);
            incl_nod(k,2)=3;
        end
    end      
    
%     for ii=1:nn(1)
%         if (nod(ii,3)==max(nod(:,3)/2) && nod(ii,4)==max(nod(:,4)/2))
%             if (nod(ii,2)>=0.01 && nod(ii,2)<=0.04) % found axis
%                 for jj=1:nn(1)
%                     dist = sqrt((nod(jj,3)-nod(incl_nod(ii,1),3)).^2 + (nod(jj,4)-nod(incl_nod(ii,1),4)).^2);
%                     if dist <=
%                     
%                     
%                     
%                 k=k+1;
%                 incl_nod(k,1)=ii;
%                 incl_nod(k,2)=2;
%             elseif (nod(ii,2)>=0.06 && nod(ii,2)<=0.09)
%                 k=k+1;
%                 incl_nod(k,1)=ii;
%                 incl_nod(k,2)=3;
%             end
%         end
%     end
        
%     for ii=1:nn(1)
%         if (incl_nod(ii,2)==2)
%             for jj=1:nn(1)
%                 dist = sqrt((nod(jj,3)-nod(incl_nod(ii,1),3)).^2 + (nod(jj,4)-nod(incl_nod(ii,1),4)).^2);
%                 if dist <= radius
%                     k = k+1;
%                     incl_nod(k,1) = ii;
%                     incl_nod(k,2) = 2;
%                 end
%             end
%         end
%     end
    psnfile(fname,incl_nod)
else
    disp('no inclusion');
    incl_nod=ones(100,2);
    psnfile(fname,incl_nod)
end
rm 3Dmesh_input.inp
end
    
function psnfile(fname,incl_nod)

psn_size = size(incl_nod);
psn = zeros(psn_size(1),3);
incl_nod = sortrows(incl_nod,1);

for ii = 1:psn_size(1)
    psn(ii,1) = ii;
    psn(ii,2:3) = incl_nod(ii,:);
end
save([fname '.psn'], 'psn','-ascii')

end

