% 本程序旨在分离downlink_model程序中的天线选择部分，使得之后的天线选择工作可以分隔完成
% 天线选择工作，输入变量应该是随机生成的AP位置矩阵和User位置矩阵，还有同时服务于某个用户的AP数目service_number
% 输出变量应该是用户与AP的对应选择矩阵new_RRH_index

function A=antenna_selection(RRH_matrix,USER_matrix,service_number)
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
new_distance_square=new_distance_square(1:service_number,1:user);
% 上一行的语句是用户到其选择的RRH之间的距离。在这个函数中，这个变量并未输出，但是后续会设计基于距离的优化算法，到那时需要将这个变量输出,用户优化
A=new_RRH_index;
B=sum(new_distance_square(:));
% new_RRH_index:无重复的服务RRH索引组成的矩阵，支持多RRH服务于一个用户。该矩阵每列即为服务于该用户的RRH编号。
% 存在的问题：单纯是使用用户的先后顺序处理RRH重复服务的问题，一开始想使用距离来处理冲突，未能实现，当很多RRH服务于用户的时候，性能可能会下降。少量RRH服务用户时暂不处理。
end

% 所有局部函数放在此行下面
function [a,b]=swap(x,y)
a=y;b=x;
end

% 2020.10.10 by LiJiaxiang