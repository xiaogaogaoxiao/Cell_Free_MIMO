% 单用户情况下，找距离最近的K个RRH
function [A,D]=choose_nearest_rrh(RRH_matrix,USER_matrix,service_number)
rrh=size(RRH_matrix,1);
user=size(USER_matrix,1);
if user~=1
    warning('本函数只能用于单用户情况下')
end
x1=USER_matrix(1,1);
y1=USER_matrix(1,2);
d=zeros(1,rrh);
for i=1:size(RRH_matrix,1)
    d(1,i)=distance(x1,y1,RRH_matrix(i,1),RRH_matrix(i,2));
end
[sv,si]=sort(d(1,:),2,'ascend');
p=si(1,1:service_number);
D=sum(sv(1,1:service_number));
A=zeros(1,rrh);
for i=1:size(p,2)
    A(1,p(1,i))=1;
end
end

function d=distance(x1,y1,x2,y2)
d=((x1-x2)^2+(y1-y2)^2)^1/2;
% if d<1
%     d=1;
% end
% if d>50
%     d=50;
% end
end