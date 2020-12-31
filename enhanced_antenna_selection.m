% 本程序为天线选择算法的提升版本
% 在之前的天线选择算法中，用户的顺序固定，这样就会导致最后的用户总是需要作出妥协
% 天线选择算法的提升版本思路为，让每个用户的选择顺序不再固定，按照各个顺序都来一遍
% 并记录每个顺序所对应的总的距离，最后的输出结果为每个用户到其服务天线组的距离的总和为最小的那个结果
% 注意事项：基本算法还是antenna_selection中的算法，但是这次需要打乱用户的顺序
% 在高层函数中，USER和RRH之间的对应关系是按照矩阵的列定义的。在打乱用户顺序并生成更好的天线选择方案时，不能忽略索引的问题
function rrh_index=enhanced_antenna_selection(RRH_matrix,USER_matrix,service_number)
[distance,rrh_index]=unenhanced_antenna_selection(RRH_matrix,USER_matrix,service_number);
rand_array=zeros(5,size(USER_matrix,1));
for i=1:5
    rand_array(i,:)=randperm(size(USER_matrix,1));
end
temp=USER_matrix;
for i=1:size(rand_array,1)-1
    for j=1:size(rand_array,2)
       USER_matrix(rand_array(i,j),:)=temp(j,:); 
    end
    [temp_distance,temp_rrh_index]=unenhanced_antenna_selection(RRH_matrix,USER_matrix,service_number);
    if temp_distance<distance
        distance=temp_distance;
        rrh_index=temp_rrh_index;
    end
end
end

function [distance,rrh_index]=unenhanced_antenna_selection(RRH_matrix,USER_matrix,service_number)
% 这个函数跟antenna_selection一样，无非是多输出了一个distance_square标量，其意思是这种天线选择方案所有用户距离其服务天线组的总距离，其目的是评价的指标
rrh=size(RRH_matrix,1);
user=size(USER_matrix,1);
distance_square=zeros(rrh,user);
RRH_index=distance_square;
for i=1:user
    for j=1:rrh
        distance_square(j,i)=(USER_matrix(i,1)-RRH_matrix(j,1))^2+(USER_matrix(i,2)-RRH_matrix(j,2))^2;
        RRH_index(j,i)=j;
    end
end

for k=1:user
    for i=1:rrh-1
        for j=1:rrh-i
            if(distance_square(j,k)>distance_square(j+1,k))
                [distance_square(j,k),distance_square(j+1,k)]=swap(distance_square(j,k),distance_square(j+1,k));
                [RRH_index(j,k),RRH_index(j+1,k)]=swap(RRH_index(j,k),RRH_index(j+1,k));
            end
        end
    end
end
% 已经完成以用户为中心，计算RRH到用户之间的距离，目前也已经实现了RRH序号索引，确保对距离RRH排序过程中，RRH索引也能跟随变化，这样能进行之后的基于索引的RRH重复服务的检测。

% 设置服务于用户的RRH为可调参数，这样可以在今后进行扩展为多RRH服务于用户。
% 如果存在RRH重复使用，对distance_square矩阵和RRH_index进行变换。如果RRH_index在一行内出现重复，则对应位置距离远的选择下一行的index进行服务；如果index中被使用了，则选择下一行的index进行服务。
already_used=[];
new_RRH_index=RRH_index;
new_distance_square=distance_square;
for k=1:user*service_number
    for i=1:user
        if(ismember(new_RRH_index(k,i),already_used)==0)
            already_used(end+1)=new_RRH_index(k,i);
        else
            while(ismember(new_RRH_index(k,i),already_used)==1)
                temp=new_RRH_index(:,i);
                temp(k)=[];
                temp(end+1)=0;
                new_RRH_index(:,i)=temp;
                
                temp_d=new_distance_square(:,i);
                temp_d(k)=[];
                temp_d(end+1)=0;
                new_distance_square(:,i)=temp_d;
                if(ismember(0,already_used)==1)
                    break
                end
            end
            already_used(end+1)=new_RRH_index(k,i);
        end
    end
end
new_RRH_index=new_RRH_index(1:service_number,1:user);
% 上一行的语句可以看做是用户和RRH之间的对应关系。第k列的元素是第k个用户所选择的RRH编号。
distance_square=new_distance_square(1:service_number,1:user);
distance=sum(distance_square(:));
% 上一行的语句是用户到其选择的RRH之间的距离。在这个函数中，这个变量并未输出，但是后续会设计基于距离的优化算法，到那时需要将这个变量输出,用户优化
rrh_index=new_RRH_index;
end

function [a,b]=swap(x,y)
a=y;b=x;
end