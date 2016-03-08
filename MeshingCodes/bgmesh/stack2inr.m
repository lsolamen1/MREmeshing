function stackparam = stack2inr(stackfn,stackinfofn,savefn)
% Tries to read the segmented MR data from brain resection project and
% write them in .inr format
% (http://serdis.dis.ulpgc.es/~krissian/InrView1/IOformat.html)
% 
% This routine assumes that segmented images are saved in '.mat' file. Also
% the information about how MR was obtained (dicom headers) is given in
% 'stackinfofn' or 'stackInfo.mat'
% 
% It returns the following properties of the dataset in a structure
% variable:
%   nRow, nCol, nPln
%   xspacing, yspacing, zspacing
%   labels (different numbers used to segment regions of the images)
% 
% Written by:
%       Hamid Ghadyani Nov 2010

if ischar(stackfn)
    S = load(stackfn);
    n = fieldnames(S);
    mask = [];
    for i=1:length(n)
        mask = getfield(S,n{i});
        if ndims(mask) == 3
            break;
        end
    end
    if isempty(mask)
        error('Input Matlab file does not contain a field that is stack of 2D masks!');
    end
    [mydir stackfn] = fileparts(stackfn);
    if isempty(mydir)
        mydir = pwd;
    end
    if nargin == 2
        savefn = [mydir filesep stackfn '.inr'];
    end
else
    mask = stackfn;
    mydir = pwd;
    if nargin == 2
        savefn = [pwd filesep 'stack2inr.inr'];
    end
end

if nargin==1 % stack info file is assumed to be saved in stackInfo.mat
    stackinfofn = 'stackInfo.mat';
end
S = load(stackinfofn);
n = fieldnames(S);
if ~isfield(S,'stackInfo')
    error(' Provided stack info file does not have a ''stackInfo'' field in its structure');
end
for i=1:length(n)
    if strcmpi('stackInfo',n{i})
        S = S.stackInfo;
        break;
    end
end
rows = S.Rows;
cols = S.Columns;
if rows~=size(mask,1) || cols~=size(mask,2)
    error('The segmented image size does not correspond to stackInfo');
end
% Transpose the matrix since fwrite will output the 3D matrix in a
% column-based fashion
if islogical(mask)
    myclass = 'uint8';
else
    myclass = class(mask);
end

% mask = smooth3(mask, 'box',3);
mask = logical(mask);
foo = uint16(mask);
foo = permute(foo, [2,1,3]);
saveinr(foo,savefn,S);

stackparam.nRow = size(mask,1);
stackparam.nCol = size(mask,2);
stackparam.nPln = size(mask,3);
stackparam.xspacing = S.PixelSpacing(1);
stackparam.yspacing = S.PixelSpacing(2);
stackparam.zspacing = S.SliceThickness;
stackparam.labels = unique(mask(:));
    
    
    
    
    
    
    
    
