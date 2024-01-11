function [c1, n1] = get_centers_and_normals(source)

nodes = source.nodes;
elems = source.elems;

%% init
c1 = zeros([size(elems,1), 3]);
n1 = c1; 

for i = 1 : size(elems,1)
     nlist = elems(i,:);
     a = nodes(nlist(1),:);
     b = nodes(nlist(2),:);
     c = nodes(nlist(3),:);
     c1(i,:) = (a+b+c)./3;
     n1(i,:) = cross( b-a, c-a )/2;
    
end
