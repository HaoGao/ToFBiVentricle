function [mapEle, VerN, TMeshN] = nodeMapAdjust(Ver, TMesh)

%%%%adjus the node number in a natural order and same for mesh 

NoOfTotalEle = size(TMesh,1);
TMeshN = zeros(size(TMesh));

for i = 1 : NoOfTotalEle
    TMeshN(i,:)=[i TMesh(i,2:5)];
    mapEle(TMesh(i,1))=i;
end

NoOfTotalNode = size(Ver,1);
VerN = zeros(size(Ver));

for i = 1 : NoOfTotalNode
    VerN(i,:)=[i Ver(i,2:4)];
    mapNode(Ver(i,1))=i;
end

%%%adjust TMeshN's no ID
for i = 1 : NoOfTotalEle
    for j = 2 : size(TMeshN,2)
        newID = mapNode(TMeshN(i,j));
        if newID<1
            disp('mapping wrong: node ID < 1');
            pause;
        end 
        TMeshN(i,j)=newID;
    end
end
