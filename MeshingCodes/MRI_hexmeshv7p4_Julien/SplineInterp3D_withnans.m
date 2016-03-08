function [ vout ] = SplineInterp3D_withnans(coord,vin,btyp,contspline)
%SplineInterp3D_withnans: Masked 3D cubic spline interpolation - does not uses NaN values for the interpolation. 
%   Interpolates a 3D stack of data using 1D cubic spline interpolation in
%   each of the 3 directions, one after the other. NaN values in the input
%   are not used in the interpolation, they are considered to be gaps
%   in the data. Interpolated values within these NaN sections are extrapolated 
%   from 'good' data and often highly inaccurate. 
%   2 options are possible when there is more than section of good data
%   across a line: 
%   1:  Continuous interpolation, where NaN's are considered bad
%       datapoints, but the end of one section and start of the next are
%       expected to be smoothly connected. Use contspline='continuous'
%   2:  Non-Continuous interpolation, where each section of data alond a line is treated
%       as independent, and indepenedent splines are fitted, the final 
%       solution is made up from a patchwork of these independent fits. This
%       is the default behavior.
%   The interpolated data is processed, setting any interpolated value which 
%   is closer to a NaN than a non-NaN datapoint to zero, because extrapolation 
%   using splines is extremly inaccurate. A mask of 'safe' values
%   can easily be created from the output using safemask=~isnan(vout).
%
%   Requires the function SplineInterp1D_withnans
%
%   INPUTS
%    xin: X locations of input data. X directionis assumed to be the 1st
%         matlab dimension, i.e. moving from vin(1,1,1) to vin(2,1,1) is a
%         movement in the x direction.
%    yin: Y locations of input data. Y directionis assumed to be the 2nd
%         matlab dimension, i.e. moving from vin(1,1,1) to vin(1,2,1) is a
%         movement in the y direction.
%    zin: Z locations of input data. Z directionis assumed to be the 3rd
%         matlab dimension, i.e. moving from vin(1,1,1) to vin(1,1,2) is a
%         movement in the z direction.
%    vin: Input data, containing NaN's where the measurements are not
%         relevant. Size=[length(xin) length(yin) length(zin)]
%    xout,yout,zout: X, Y and Z locations of interpolated data.
%    btyp(optional): Type of boundary condition applied for cubic spline fit. 1=
%          natural spline condition, v''(end)=0. 2= fixed gradient conditions,
%          v'(end) = gradient given by last 2 data points.
%          v'(start) = gradient given by first 2 data points.
%    contspline(optional): Set to 'continuous; to use a continuous spline
%          across a line of data with breaks in it. Otherwise, an indpenedent 
%          spline is built for eac continuous section of data.
%   OUTPUT
%    vout: Interpolated values of vin, at the points supplied by xout,yout,zout.
%          Values of vin which are NaN's are not used in the interpolation.
%          If [xout(i),yout(j),zout(k)] is closer to a NaN value in vin,
%          vout(i,j,k)=NaN. This is because extrapolation using splines is
%          very innacurate. Note that this distance measurement is not an
%          absolute 3D distance, because each dimension is processed
%          seperatly, it checks the z distance first, then the y distance,
%          then the x distance. 
%
%   Author: Matt Mcgarry  -matthew.d.mcgarry@darmtouth.edu
%           Thayer School of Engineering
%           Dartmouth College, Hanover NH.
%     Date: 22 Dec 2011.
%
%   Example:
%   load MRE_3DMotionData.mat;
%   load Mask.mat
%   Ur=A(:,:,:,1).*cos(P(:,:,:,1));
%   Ur(mask~=1)=nan; % <- Ur is just a masked 3D stack of data, values outside mask=NaN.
%   xin=1:size(Ur,1);yin=1:size(Ur,2);zin=1:size(Ur,3);
%   xout=1:0.7:size(Ur,1);yout=1:0.7:size(Ur,2);zout=1:0.7:size(Ur,3);
%   %EXAMPLE 1: natural spline BC's, process each continuous section along a line individually 
%   Ur_interp=SplineInterp3D_withnans(xin,yin,zin,Ur,xout,yout,zout,1,'no')
%   %EXAMPLE 2: Fixed gradient spline BC's, process each continuous section along a line individually. 
%   Ur_interp=SplineInterp3D_withnans(xin,yin,zin,Ur,xout,yout,zout,2,'no')
%   %EXAMPLE 3: natural spline BC's, process multiple continuous sections along a line as a continuous spline. 
%   Ur_interp=SplineInterp3D_withnans(xin,yin,zin,Ur,xout,yout,zout,1,'continuous')
%   safemask=~isnan(Ur_interp).

%% Process inputs
if(nargin<8) % neither btyp or contspline is supplied
    btyp=1; % 
    contspline='no';
elseif(nargin<9)% either btyp or contspline is supplied
    if(ischar(btyp)) % 4th argument is contspline
        contspline=btyp;
        btyp=1;
    else %4th argument is btyp;
        contspline='no';
    end
end

if(numel(coord.xin)>length(coord.xin))
    error('xin must be a 1D vector')
end
if(numel(coord.yin)>length(coord.yin))
    error('yin must be a 1D vector')
end
if(numel(coord.zin)>length(coord.zin))
    error('zin must be a 1D vector')
end
if(numel(coord.xout)>length(coord.xout))
    error('xout must be a 1D vector')
end
if(numel(coord.yout)>length(coord.yout))
    error('yout must be a 1D vector')
end
if(numel(coord.zout)>length(coord.zout))
    error('zout must be a 1D vector')
end
    

% ensure vin is correct size
if(size(vin,1)~=length(coord.xin))||(size(vin,2)~=length(coord.yin))||(size(vin,3)~=length(coord.zin))
    error('vin must be size [length(xin) length(yin) length(zin)]')
end

dim=size(vin);
dimout=[length(coord.xout) length(coord.yout) length(coord.zout)];

% Slice direction first
VoutZ=zeros([dim(1) dim(2) dimout(3)]);
for ii=1:dim(1)
    for jj=1:dim(2)
        VoutZ(ii,jj,:)= SplineInterp1D_withnans(coord.zin,vin(ii,jj,:),coord.zout,btyp,contspline);
        VoutZ(ii,jj,:)= Checkfornans(coord.zin,vin(ii,jj,:),coord.zout,VoutZ(ii,jj,:));
        
    end
end

% Now the second index
VoutYZ=zeros([dim(1) dimout(2) dimout(3)]);
for ii=1:dim(1)
    for jj=1:dimout(3)
        VoutYZ(ii,:,jj)=SplineInterp1D_withnans(coord.yin,VoutZ(ii,:,jj),coord.yout,btyp,contspline);
        VoutYZ(ii,:,jj)=Checkfornans(coord.yin,VoutZ(ii,:,jj),coord.yout,VoutYZ(ii,:,jj));
    end
end
clear VoutZ

% Now the final index
vout=zeros([dimout(1) dimout(2) dimout(3)]);
for ii=1:dimout(2)
    for jj=1:dimout(3)
        vout(:,ii,jj)=SplineInterp1D_withnans(coord.xin,VoutYZ(:,ii,jj),coord.xout,btyp,contspline);
        vout(:,ii,jj)=Checkfornans(coord.xin,VoutYZ(:,ii,jj),coord.xout,vout(:,ii,jj));
    end
end
clear VoutYZ

end

function vout2= Checkfornans(xin,vin,xout,vout)
%Checkfornans: Checks 1D interpolated values for proximity to NaN in the
%input. If the interpolated value is closer to a NaN than real data, sets
%it to NaN.

vout2=vout;
for ii=1:length(xout)
    [y,i]=min(abs(xout(ii)-xin));
    if(isnan(vin(i)))
        vout2(ii)=nan;
    end
end

end

