RRH_matrix=rand(10,2)*100;
USER_matrix=rand(3,2)*40+ones(3,2)*30;  % 在（10,10）到（90，90）的范围内生成用户坐标
Nrf=1;service_number=2;Nt=4;Nr=1;
yita=0.1:0.1:50;
rrh=size(RRH_matrix,1);
user=size(USER_matrix,1);
LOS_matrix=zeros(user,rrh);     % 该矩阵用于记录第i个用户到第j个RRH之间有无直射路径
for i=1:user
    for j=1:rrh
        d=distance(USER_matrix(i,1),USER_matrix(i,2),RRH_matrix(j,1),RRH_matrix(j,2));
        LOS_matrix(i,j)=rand(1)<min(20/d,1)*(1-exp(-d/39))+exp(-d/39);
    end
end
channel_cell=cell(user,rrh);    % 该cell用于存储第i个用户到第j个RRH之间的平均信道，即平均掉小尺度的信道
precoder_cell=cell(user,rrh,2);
power_cell=cell(user,rrh);
loop=20;                        % 该参数用于更改平均计算所用的循环次数
for i=1:user
    for j=1:rrh
        for l=1:loop
            if size(channel_cell{i,j},1)==0||size(channel_cell{i,j},2)==0   % 新生成的cell行数和列数都为0，不进行重新赋值操作会无法与new_channel矩阵相加
                channel_cell{i,j}=zeros(Nt,Nr);
            end
            channel_cell{i,j}=channel_cell{i,j}+new_channel(USER_matrix(i,1),USER_matrix(i,2),RRH_matrix(j,1),RRH_matrix(j,2),Nt,Nr,LOS_matrix(i,j));
        end
        channel_cell{i,j}=channel_cell{i,j}/loop;   % 平均掉小尺度衰落，得到从第i个用户到第j个RRH之间的平均信道
        [Vrf,Vb]=precoder(channel_cell{i,j},Nrf,3);
        precoder_cell{i,j,1}=Vrf;
        precoder_cell{i,j,2}=Vb;
        power_cell{i,j}=channel_cell{i,j}'*Vrf*(Vb*Vb')*Vrf'*channel_cell{i,j};
    end
end


function d=distance(x1,y1,x2,y2)
d=((x1-x2)^2+(y1-y2)^2)^1/2;
end

function [Vrf,Vb]=precoder(H,Nrf,type)
if type==1 % 全随机波束成形
    Vrf=rand(size(H,1),Nrf);
    for i=1:size(Vrf,1)
        for j=1:size(Vrf,2)
            w=2*pi*rand;
            Vrf(i,j)=exp(w*1i);
        end
    end
    Vb=(H'*Vrf)'/norm(H'*Vrf,2);
elseif type==2 % 无约束的模拟预编码
    [~,~,V]=svd(H');
    transpose_V=V';
    Vrf=transpose_V(:,1:Nrf);
    Vrf=Vrf/norm(Vrf);
    Vb=(H'*Vrf)'/norm(H'*Vrf,2);
    Vb=Vb/norm(Vb);
elseif type==3 % 有约束的模拟预编码
    [~,~,V]=svd(H');
    transpose_V=V';
    Vrf=angle(transpose_V(:,1:Nrf));
    Vrf=Vrf/norm(Vrf);
    Vb=(H'*Vrf)'/norm(H'*Vrf,2);
    Vb=Vb/norm(Vb);
end
end