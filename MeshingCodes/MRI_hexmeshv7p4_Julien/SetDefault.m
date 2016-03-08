function [default]=SetDefault
default.mridef='MRE_3DMotionData.mat';
default.Hdrdef='HeaderData.mat';
default.mskdef='Mask.mat';
default.meshstratdef=1;
default.zonestratdef=2;
default.mudef=3300;             % Default shear modulus
default.npwdef=12;              % Nodes  per wavelength
default.rhoest=1000;
default.npzdef=3500;            % Node per zone
default.wlperzonedef=0.7;
default.resdef=[0.002 0.002 0.002];
default.DR=0.18; % Default damping ratio
%DR=0.10; % Default damping ratio
default.vincomp=0.499; % Default Incompressible Damping Ratio
default.vcomp=0.4; % Default compressible damping ratio
end