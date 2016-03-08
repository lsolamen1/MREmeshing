This folder contains codes to create tetrahedral mesheson a brain after a resection operation.

To use the code, you need to run 'ResectBrainUseMRStack.m' in Matlab.
Make sure this folder and its sub-directories are in Matlab's search path.

Within Matlab you can enter:

% [e p] = ResectBrainUseMRStack(stackfn, stackinfofn, ...
%                         stereovisionfn, resection_normal, ...
%                         outputfn, meshing_criteria)
% 
% This routine reads the segmented images of a brain (stored in stackfn and
% stackinfofn), creates a surface mesh and then removes part of the
% surface mesh using the 'stereovisionfn'. The resulting surface mesh is
% then used to create a tetrahedral mesh (saved in 'outputfn')
% 
% It will first compute the surface of the tetrahedral mesh and then
% perform a boolean operation between that and the stereovision surface
% using CGAL Nef polyhedrons. 

The first two arguments provide that path to .mat files containing the MR image stack.
Third argument is the file name containing the stereo-vision surface used for resection.

Each file is well documented so you can customize it to produce different mesh resolutions.

