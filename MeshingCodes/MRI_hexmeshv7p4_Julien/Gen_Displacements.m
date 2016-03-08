function [deplacement,interpo,coord,meshprop]...
    =Gen_Displacements(par,A,P,mask,MagIm,strat,vox,regs,dim,meshprop,noreg)
if(par)
    matlabpool('open',nlabs);
end
% Generate Real and Imag Displacements and apply mask
Ur=A.*cos(P);
Ui=A.*sin(P);
clear A P

% % Apply a masked selective median filter to clear up outliers
% thresh=0.2; % Thereshold to apply median filter
% disp(['Median filter with threshold ' num2str(thresh)])
% for ii=1:3
%     Ur(:,:,:,ii)=selectivemedianfilter(Ur(:,:,:,ii),mask,thresh);
%     Ui(:,:,:,ii)=selectivemedianfilter(Ui(:,:,:,ii),mask,thresh);
% end

Ur(repmat(mask,[1 1 1 3])==0)=nan;
Ui(repmat(mask,[1 1 1 3])==0)=nan;
coord.xin=0:vox(1):vox(1)*(dim(1)-1);
coord.yin=0:vox(2):vox(2)*(dim(2)-1);
coord.zin=0:vox(3):vox(3)*(dim(3)-1);

if(strat.mesh==1) % Each MR voxel is a FE node
    interpo.maskint=mask;
    deplacement.Urint=Ur;
    deplacement.Uiint=Ui;
    interpo.MagIm_int=MagIm;
    clear Ur Ui MagIm
    
    coord.xout=coord.xin;
    coord.yout=coord.yin;
    coord.zout=coord.zin;
    interpo.regs_int=regs;
    meshprop.meshres=vox;
elseif(strat.mesh==2)||(strat.mesh==3)  
    %  Process geometry and interpolate to get desired nodes per wavelength if meshstrat = 2.
    if(strat.mesh)==2
        meshprop.targetres(1:3) = meshprop.wl/meshprop.npw;
        % targetres already defined from inputs if meshstrat==3
    end
    disp(['Target Resolution: ' num2str(meshprop.targetres)])
    
    % Allow some deviation from the target resolution to fit an integer
    % number of Hex27 elements across the geometry. This will eliminate
    % wasting planes of data when there is an even number of slices which
    % occurs using a voxel based Hex27 meshing scheme.
    
    Imask=find(mask);
    [mi,mj,mk]=ind2sub(dim,Imask);
    ext(1)=(max(mi)-min(mi))*vox(1);
    ext(2)=(max(mj)-min(mj))*vox(2);
    ext(3)=(max(mk)-min(mk))*vox(3);
    
    % Number of targetres fitting across ext = ext(ii)/targetres. Want an
    % EVEN number of resolutions across the domain to fit 27 node elements
    % (Even number of resolutions = odd number of nodes)
    meshprop.meshres=zeros(1,3);
    for ii=1:3
        nres=floor(ext(ii)/meshprop.targetres(ii));
        if(mod(nres,2)==1)
            nres=nres+1;
        end
        if nres<2
            nres=2;
            disp('ERROR: Mesh Resolution Larger Than Data Set!!!');
        end
        meshprop.meshres(ii) = ext(ii)/nres;
    end
    disp(['Actual New Resolution(m): ' num2str(meshprop.meshres)])
    disp(['Data Resolution(m): ' num2str(vox)])
    
    % Warning if mesh is being interpolated to a lower resolution than the
    % data, this is not really a good idea.
    if(meshprop.meshres(1)>vox(1)||meshprop.meshres(2)>vox(2)||meshprop.meshres(3)>vox(3))
        disp([' >> Warning: Mesh is being interpolated to a lower resolution than the data'])
    end
    
    % Perform Interpolation
    
    xout_tmp= (min(mi)-1)*vox(1):meshprop.meshres(1):(max(mi)-1)*vox(1);
    yout_tmp= (min(mj)-1)*vox(2):meshprop.meshres(2):(max(mj)-1)*vox(2);
    zout_tmp= (min(mk)-1)*vox(3):meshprop.meshres(3):(max(mk)-1)*vox(3);
    
    % Include a buffer around the extents of the mask so that MagIm can be
    % interpolated, and the dataset extends past the boundaries of the
    % mask. e.g. for brain data, the buffer means the interpolated MagIm
    % shows the skull etc, not just the masked region of the brain.
    bufsiz=[4 4 0]; % Size of buffer in each direction (interpolated MR voxels) MUST BE EVEN
    if(sum(mod(bufsiz,2)>0))
        error('bufsiz must be even to fit hex27 elements efficiently inside mask extents!!')
    end
    
    coord.xout=zeros(1,2*bufsiz(1)+length(xout_tmp));
    coord.xout(bufsiz(1)+1:bufsiz(1)+length(xout_tmp))=xout_tmp;
    for ii=1:bufsiz(1)
        coord.xout(bufsiz(1)+1-ii)=coord.xout(bufsiz(1)+1)-meshprop.meshres(1)*ii;
        coord.xout(bufsiz(1)+length(xout_tmp)+ii)=coord.xout(bufsiz(1)+length(xout_tmp))+meshprop.meshres(1)*ii;
    end
    
    coord.yout=zeros(1,2*bufsiz(2)+length(yout_tmp));
    coord.yout(bufsiz(2)+1:bufsiz(2)+length(yout_tmp))=yout_tmp;
    for ii=1:bufsiz(2)
        coord.yout(bufsiz(2)+1-ii)=coord.yout(bufsiz(2)+1)-meshprop.meshres(2)*ii;
        coord.yout(bufsiz(2)+length(yout_tmp)+ii)=coord.yout(bufsiz(2)+length(yout_tmp))+meshprop.meshres(2)*ii;
    end
    
    coord.zout=zeros(1,2*bufsiz(3)+length(zout_tmp));
    coord.zout(bufsiz(3)+1:bufsiz(3)+length(zout_tmp))=zout_tmp;
    for ii=1:bufsiz(3)
        coord.zout(bufsiz(3)+1-ii)=coord.zout(bufsiz(3)+1)-meshprop.meshres(3)*ii;
        coord.zout(bufsiz(3)+length(zout_tmp)+ii)=coord.zout(bufsiz(3)+length(zout_tmp))+meshprop.meshres(3)*ii;
    end
    
    if(0==1) % Check locations of xin and xout
        figure
        vjunk=zeros(size(coord.xin));vjunk(min(mi):max(mi))=1;
        plot(coord.xin,vjunk,'r-',coord.xin,zeros(size(coord.xin)),'ro',coord.xout,zeros(size(coord.xout)),'b.')
        
        coord.xout(1)
        coord.xin(min(mi))
        
        coord.xout(end)
        coord.xin(max(mi))
        coord.xout;
    end
    
    % Interpolate Ur, Ui and MagIm with the function SplineInterp3D_withnans.
    deplacement.Urint=zeros([length(coord.xout) length(coord.yout) length(coord.zout) 3]);
    deplacement.Uiint=zeros([length(coord.xout) length(coord.yout) length(coord.zout) 3]);
    %MagIm_int=zeros([length(xout) length(yout) length(zout) 1]);
    if(par)
        n_int=8;
        needsint=zeros([size(MagIm) n_int]);
        needsint(:,:,:,1:3)=Ur;
        needsint(:,:,:,4:6)=Ui;
        needsint(:,:,:,7)=MagIm;
        needsint(:,:,:,8)=regs;
        
        intdata=zeros([length(coord.xout) length(coord.yout) length(coord.zout) n_int]);
        
        parfor ii=1:n_int
            intdata(:,:,:,ii)=SplineInterp3D_withnans(coord,needsint(:,:,:,ii),2,'no');
        end
        deplacement.Urint(:,:,:,:)=intdata(:,:,:,1:3);
        deplacement.Uiint(:,:,:,:)=intdata(:,:,:,4:6);
        interpo.MagIm_int=intdata(:,:,:,7);
        interpo.regs_int=intdata(:,:,:,8);
        clear intdata needsint
    else
        for ii=1:3
            deplacement.Urint(:,:,:,ii)=SplineInterp3D_withnans(coord,Ur(:,:,:,ii),2,'no');
        end
        interpo.MagIm_int=SplineInterp3D_withnans(coord,MagIm,2,'no');
        
        if(noreg)
            interpo.regs_int=zeros(size(interpo.MagIm_int));
        else
            interpo.regs_int=SplineInterp3D_withnans(coord,regs,2,'no');
        end
    end
    
    %     if(par)
    %         parfor ii=1:3
    %             Urint(:,:,:,ii)=SplineInterp3D_withnans(xin,yin,zin,Ur(:,:,:,ii),xout,yout,zout,2,'no');
    %         end
    %     else
    %         for ii=1:3
    %             Urint(:,:,:,ii)=SplineInterp3D_withnans(xin,yin,zin,Ur(:,:,:,ii),xout,yout,zout,2,'no');
    %         end
    %     end
    disp('Real Displacements Interpolated')
    % Interpolate Ui with the function SplineInterp3D_withnans.
    %    Uiint=zeros([length(xout) length(yout) length(zout) 3]);
    %     if(par)
    %         parfor ii=1:3
    %             Uiint(:,:,:,ii)=SplineInterp3D_withnans(xin,yin,zin,Ui(:,:,:,ii),xout,yout,zout,2,'no');
    %         end
    %     else
    %         for ii=1:3
    %             Uiint(:,:,:,ii)=SplineInterp3D_withnans(xin,yin,zin,Ui(:,:,:,ii),xout,yout,zout,2,'no');
    %         end
    %     end

    disp('Imag Displacements Interpolated')
    interpo.maskint=~isnan(deplacement.Urint(:,:,:,1));
    clear MagIm Ui Ur
    disp('MagIm Interpolated')
    if(0==1)
        for ii=1:3
            %montagestack(Ur(:,:,:,ii))
            montagestack(deplacement.Urint(:,:,:,ii))
        end
        montagestack(interpo.MagIm_int)
    end
    
end
if(par)
    matlabpool('close');
end
end