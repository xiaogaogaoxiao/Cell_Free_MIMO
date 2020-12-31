function d_attenuation=fitness_based_distance_for_multiple_user(yita,x,service_number,RRH_matrix,USER_matrix)
if sum(ismember(sum(x,1),1))~=size(USER_matrix,1)*service_number
    d_attenuation=0;
    return
end
% 上一个if用于判断，每个天线仅可被选择一次，即将矩阵行与行相加，最后的向量仅有1和0元素
% 但是可能存在每个用户选的天线个数不是service_number的情况
if ~isequal(sum(x,2),service_number*ones(size(USER_matrix,1),1))
    d_attenuation=0;
    return
end
% 上一个if用于判断每个用户选择service_number个天线，集将矩阵列与列相加，最后表现为每个元素都为service_number的列矩阵
useful_distance=zeros(size(x,1),1);interference_distance=zeros(size(x,1),1);
for i=1:size(x,1)
    useful_distance(i,1)=yita*get_distance_attenuation(x,i,RRH_matrix,USER_matrix,1);
    interference_distance(i,1)=yita*get_distance_attenuation(x,i,RRH_matrix,USER_matrix,0);
end
d_attenuation=sum(useful_distance)/sum(interference_distance);
end

function d=get_distance_attenuation(x,row,RRH_matrix,USER_matrix,type)
d=0;
if type==1
    for j=1:size(x,2)
        if x(row,j)==1
            d=d+distance(USER_matrix(row,1),USER_matrix(row,2),RRH_matrix(j,1),RRH_matrix(j,2));
        end
    end
else
    for i=1:size(x,1)
        if i~=row
            for j=1:size(x,2)
                d=d+distance(USER_matrix(i,1),USER_matrix(i,2),RRH_matrix(j,1),RRH_matrix(j,2));
            end
        end
    end
end
end

function d=distance(x1,y1,x2,y2)
% 改造之后的距离：距离的-3.19次方
constant=-3.19;
d=((x1-x2)^2+(y1-y2)^2)^1/2;
d=d^constant;
end