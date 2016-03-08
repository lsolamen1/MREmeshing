function [e p] = bandwidth_reduction(oe,op)

fprintf('Optimizing Mesh...');

ni = zeros(9*size(oe,1),1);
nj = zeros(9*size(oe,1),1);
nk = 1;
for i = 1 : length(oe)
    for j = 1:3
        for k = 1:3
            ni(nk) = oe(i,k);
            nj(nk) = oe(i,j);
            nk = nk + 1;
        end
    end
end
gr = sparse(ni,nj,ones(size(ni,1),1),size(op,1),size(op,1));
% figure
% spy(gr)

NodeSort = symrcm(gr);
for i = 1 : length(op)
    invsort(NodeSort(1,i)) = i;
end

p = op(NodeSort,:);

elm_new = zeros(size(oe,1),size(oe,2));
for i = 1 : length(oe);
    elm_new(i,:) = invsort(oe(i,:));
end
e = elm_new;

ni = zeros(9*size(e,1),1);
nj = zeros(9*size(e,1),1);
nk = 1;
for i = 1 : length(e)
    for j = 1:3
        for k = 1:3
            ni(nk) = e(i,k);
            nj(nk) = e(i,j);
            nk = nk + 1;
        end
    end
end
gr = sparse(ni,nj,ones(size(ni,1),1),size(op,1),size(op,1));
% figure
% spy(gr)

fprintf(' Complete.\n');

