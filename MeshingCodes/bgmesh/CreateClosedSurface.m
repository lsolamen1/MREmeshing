function [ee pp] = CreateClosedSurface(es,ps,microscope_axis)
% Tries to create a closed surface based on the stereo vision camera
% surface.
% 'microscope_axis' is a 2x3 matrix containing the coordinates of the tip and hind of
% the microscope in mm. First row is assumed to be the 'tip'
% 
% Written by:
%  Hamid Ghadyani, Nov 2010

fprintf(' Creating a closed surface from stereovision mesh.\n');

[ee pp] = UseMicroscopeAxis(es,ps,microscope_axis);
% [ee pp] = UsePlaneFitting(es,ps);

ee = FixPatchOrientation(pp,ee,[],1);

input_args.verbose=0;
input_args.type=1;
[junk,junk,junk,myst] = CheckMesh3D(ee,pp,[],input_args);
if isfield(myst,'b') && myst.b~=0 && myst.b~=4
    error(' The stereovision surface is not closed, single material or manifold!\n');
end


function [ee pp] = UseMicroscopeAxis(es,ps,microscope_axis)
% microscope_axis contains information about the normal to brain surface
% before surgery and center of craniotomy

medialaxis_factor = 12;
% Compute the plane normal to 'microscope_axis'
n = microscope_axis(1,:);
n = n/norm(n);
A = microscope_axis(2,:);
A = A + n*medialaxis_factor;

% Get all the nodes on the boundary of stereovision surface and remove one
% ring of elements off of the edges of stereovision surface
edges = boundedges(ps,es);
bdynodes=unique(edges(:));
[tf idx]=ismember(es,bdynodes);
bf=sum(tf,2);
ee=es(bf==0,:);

es = ee;
edges = boundedges(ps,es);
bdynodes=unique(edges(:));
[tf idx]=ismember(es,bdynodes);
bf=sum(tf,2);
ee=es(bf==0,:);


% Renumber the nodes in original stereovision surface since we removed some
% of the bad elements.
nodes = unique(ee(:));
[tf ee] = ismember(ee,nodes);
pp = ps(nodes,:);
edges = boundedges(pp,ee);
bdynodes = unique(edges(:));

% Project the points on the normal plane
nlayers = 12;
te = ee; tp = pp;
tA = microscope_axis(2,:);
% figure; hold on

for i=1:nlayers
    
    nodes = unique(te(:));
    [tf te] = ismember(te,nodes);
    tp = tp(nodes,:);
    edges = boundedges(tp,te);
    bdynodes = unique(edges(:));
    
    bdyp = tp(bdynodes,:);
    PA = bdyp - repmat(tA,size(bdyp,1),1);
    proj = bdyp - repmat(dot(PA,repmat(n,size(bdyp,1),1),2),1,3).*repmat(n,size(bdyp,1),1);
    
    [tf idx] = ismember(edges,bdynodes);
    idx = idx + size(tp,1);
    newfacets = [edges idx(:,1)];
    newfacets = [newfacets; idx edges(:,2)];
    te = [te; newfacets];
    tp = [tp;proj];
    tA = tA + n*medialaxis_factor/nlayers;
%     trisurf(te,tp(:,1),tp(:,2),tp(:,3))
end

% bdyp = pp(bdynodes,:);
% PA = bdyp - repmat(A,size(bdyp,1),1);
% proj = bdyp - repmat(dot(PA,repmat(n,size(bdyp,1),1),2),1,3).*repmat(n,size(bdyp,1),1);
% 
% % Now add new facets
% [tf idx] = ismember(edges,bdynodes);
% idx = idx + size(pp,1);
% newfacets = [edges idx(:,1)];
% newfacets = [newfacets; idx edges(:,2)];
% ee = [ee;newfacets];
% pp = [pp;proj];

% Add a supernode to close the surface
supernode = tA+n*medialaxis_factor/nlayers;

pp = tp; ee = te;
edges = boundedges(pp,ee);
snid = size(pp,1)+1;
newfaces = [edges ones(size(edges,1),1)*snid];

pp = [pp; supernode];
ee = [ee; newfaces];
% trisurf(ee,pp(:,1),pp(:,2),pp(:,3),ones(size(pp,1),1),'FaceAlpha',0.8)
% writenodelm_surface_medit('closed.mesh',ee,pp);

function [ee pp] = UsePlaneFitting(es,ps)
n = 1:250;
m = (size(ps,1)-250):size(ps,1);
idx = [n' ; m'];
% idx = 1:size(ps,1);

% ps = ps - repmat(mean(ps),size(ps,1),1);
x = ps(idx,1);
y = ps(idx,2);
z = ps(idx,3);

% Fit a plane to these nodes so we can align the axis of 'resection' hole
% with 'z' axis;

ft = fittype( 'poly11' );
opts = fitoptions( ft );
opts.Weights = zeros(1,0);
[fitresult, gof] = fit( [x, y], z, ft, opts );

% Compute rotation angles
V = [fitresult.p10 fitresult.p01 -1.0]; D = fitresult.p00;
U=-V/norm(V);
a=U(1);b=U(2);c=U(3);
d=norm([0 b c]);
alpha=acos(c/d)*180/pi;
beta=asin(-a)*180/pi;

figure
h = trisurf(es,ps(:,1),ps(:,2),ps(:,3)); hold on
x1=40;y1=50;x2=60;y2=150;x3=40;y3=150;
z1=fitresult.p00 + fitresult.p10*x1 + fitresult.p01*y1;
z2=fitresult.p00 + fitresult.p10*x2 + fitresult.p01*y2;
z3=fitresult.p00 + fitresult.p10*x3 + fitresult.p01*y3;
patch([x1 x2 x3],[y1 y2 y3],[z1 z2 z3],0)

rotate(h,[1 0 0],alpha,[0 0 0]);
rotate(h,[0 1 0],beta,[0 0 0]);

% Calculate a node which is way above the 'flat' section of stereo surface
% move the poin
newp = get(h,'Vertices');
center=mean(newp);
z=abs(max(newp(:,3)))+200;
supernode = [center(1) center(2) z];


% Add new triangles to close the surface
edges = boundedges(ps,es);
snid = size(ps,1)+1;
newfaces = [edges ones(size(edges,1),1)*snid];

pp = [newp; supernode];
ee = [es; newfaces];
close(gcf);

h = trisurf(ee,pp(:,1),pp(:,2),pp(:,3));

% Rotate back to original location
rotate(h,[0 1 0],-beta,[0 0 0]);
rotate(h,[1 0 0],-alpha,[0 0 0]);

pp = get(h,'Vertices');
close(gcf);

% writeOFF('closed.off',ee,pp);

