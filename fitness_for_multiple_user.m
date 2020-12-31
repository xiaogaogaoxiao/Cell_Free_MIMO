function R=fitness_for_multiple_user(yita,x,service_number,power_cell)
if max(max(x))>1                % x不再是0-1矩阵
    R=0;
    warning('迭代过程中x出现非0非1元素')
    return
end
% if max(sum(x,1))>Nrf || max(sum(x,2))>service_number || max(max(x))==0
if max(sum(x,1))>1 || max(sum(x,2))>service_number || max(max(x))==0
    R=0;
    return
end
N0=10^(-143/10)/1000;useful_power=zeros(size(x,1),1);interference_power=zeros(size(x,1),1);
for i=1:size(x,1)
    useful_power(i,1)=yita*get_power(x,i,power_cell,1);
    interference_power(i,1)=yita*get_power(x,i,power_cell,0);
end
SINR=sum(useful_power)/(N0+sum(interference_power));
R=log(1+SINR);
end

function power=get_power(x,row,power_cell,type)
% x:个体矩阵，某种天线选择方式
% row:行指示器，显示现在要算第几个用户的有用或干扰功率
% channel_cell:之前生成的user*rrh维度的收发之间的信道
% type:等于1时计算有用功率，等于0时计算干扰功率
power=0;
if type==1
    for j=1:size(x,2)
        if x(row,j)==1
            power=power+power_cell{row,j};
        end
    end
else
    for i=1:size(x,1)
        if i~=row
            for j=1:size(x,2)
                power=power+power_cell{i,j};
            end
        end
    end
end
end