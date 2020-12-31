% 本程序旨在绘制一个baseline情况下的系统下行总速率和yita之间关系的图像
% 考虑Nrf=1的情况，一个RRH至多服务于一个用户，一个用户可以选择service_number个RRH
% 当出现冲突时，RRH随机选择一个欲选择它的用户。没有被选择的用户继续进行下一次RRH选择
function R=baseline(yita,RRH_matrix,USER_matrix,service_number,power_cell)
rrh=size(RRH_matrix,1);
user=size(USER_matrix,1);
R=0;loop=20;    % loop用于打乱USER_matrix和RRH_matrix的行的顺序后取平均
for i=1:loop
    new_USER_matrix=USER_matrix(randperm(size(USER_matrix, 1)),:);
    A=antenna_selection(RRH_matrix,new_USER_matrix,service_number);
    R=R+norm(get_capacity(A,yita,rrh,user,power_cell));
end
R=R/loop;
end

function R=get_capacity(A,yita,rrh,user,power_cell)
x=zeros(user,rrh);  % 建立一个二进制的天线选择方案矩阵，好利用之前写的函数
for i=1:user
    for j=1:size(A,1)
        x(i,A(j,i))=1;
    end
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