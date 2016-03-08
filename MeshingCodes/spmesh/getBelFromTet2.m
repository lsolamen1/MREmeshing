function BEL = getBelFromTet2 (seriesname)
%
% function BEL = getBelFromTet (seriesname)
%
% Pure matlab version that extracts the boundary elements from a given
% connectivity of tetrahedral elements TET (N-by-4).  Returns BEL (M-by-3),
% and IN(M-by-2). If IN(k,:) == (i,j), it indicates that the kth bel
% element comes from the ith tet element, and jth node arangement order (1:
% 123; 2: 234; 3:341; 4:412).
%
%
% See also getBelFromTri, getBelFromTri2, getBelFromTet
%
% Songbai Ji (09/01/2006).
% Songbai Ji (11/14/2006).
% Songbai Ji (8/7/2008). Bug with hist edge definition fixed.

name1 = load([seriesname,'.elm']);
name2 = load([seriesname,'.nod']);

TET_1 = name1(:,2:5);
NOD = name2(:,2:4);

faces = [ TET_1(:, [1 2 3]); TET_1(:,[2 3 4]); TET_1(:,[3 4 1]); TET_1(:,[4 1 2])];

%% create hash function
hash = sin(faces(:,1)) + sin(faces(:,2)) + sin(faces(:,3));

% the hash function must make sure shash is the same as uhash, because
% unique() is generally slow, we use sorted hash as the unique set.  Since
% hist() will sort the results, we need to pre-sort the edges/hash values
% so that we can get the correct indices for boundary edges (bel).

% uhash = unique(hash); 
[shash, I]=sort(hash);
faces = faces(I,:); % re-order faces

%%
N=hist( shash, shash+10*eps ); 
BEL = faces (N==1,:); %% those occur only once

A = linspace(1,length(BEL),length(BEL));
 
% BEL1 = zeros(length(BEL)+1,4);
% BEL1(1,:) = [1 3 2 length(BEL)]; 
% BEL1(2:length(BEL)+1,:) = [A' BEL];

BEL1 = zeros(length(BEL),4);
BEL1(1:length(BEL),:) = [A' BEL];

BNOD_index = unique(BEL1(2:end,2:4));

for i = 1:length(BNOD_index)
    BNOD_val(i,1:3) = name2(BNOD_index(i),1:3);
end


B = linspace(1,length(BNOD_index),length(BNOD_index));

% BNOD = zeros(length(BNOD_index)+1,4);
% BNOD(1,:) = [length(BNOD_index) 1 0 0];
% BNOD(2:length(BNOD_index)+1,:) = [B' BNOD_val];

BNOD = zeros(length(BNOD_index),4);
BNOD(1:length(BNOD_index),:) = [B' BNOD_val];


save([seriesname, '.bel'], 'BEL1', '-ascii');
save([seriesname, '.bnod'], 'BNOD', '-ascii');
