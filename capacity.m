% 本程序旨在生成service_number变量与第k个用户所获得的最大下行速率之间的关系
% service_number是被选中的用于服务与第k个用户的RRH数目
% 天线选择采取antenna_selection.m中的基于距离选择，且采用简单的队列法
% 信道模型采用channel.m文件中的s-v模型
% 其中，Vrf为随机矩阵且每个元素模长为1，也就是说由e的jw次方构成；Vb为随机矩阵
% 混合波束成形的两个矩阵待升级
% 基站端发射天线数目为4，用户接收端天线数目为1

function R_average=capacity(service_number,USER_matrix)
Ns=1;Nrf=2;Nt=4;Nr=1;N0=10^(-143/10)/1000;
%发送信号是1*1维向量，Ns表示发射功率，2个RF链路，4个发射天线，用于生成后面的随机发送信号；N0是噪声高斯分布方差
% 长宽为100的方形区域内有5个用户和200个AP
user=size(USER_matrix,1);
rrh=60;range=100;
R=zeros(1,user);
%service_number=1;
%service_number变量为一个用户需要几个RRH为其服务
loop=10;
R_total=zeros(loop,user);
% 矩阵中的每行即为RRH的坐标
% USER_matrix=rand(user,2)*range;
% 矩阵中的每行即为USER的坐标
for loop_time=1:loop
    Vb=randn(Nrf,Ns);
    Vrf=rand(Nt,Nrf);
    for i=1:size(Vrf,1)
        for j=1:size(Vrf,2)
            w=2*pi*rand;
            Vrf(i,j)=exp(w*1i);
        end
    end
    V=(Vrf*Vb)/norm(Vrf*Vb);
    RRH_matrix=rand(rrh,2)*range;
    for k=1:user
        % 研究第k个用户的下行容量
        % 先计算useful_power，这个量需要知道哪些天线被选中，然后生成被选中的天线到第k个用户之间的信道矩阵
        new_RRH_index=enhanced_antenna_selection(RRH_matrix,USER_matrix,service_number);
        useful_RRH=new_RRH_index(:,k); % 该矩阵所含元素为发射对第k个用户的有用信号的RRH编号
        interference_RRH=new_RRH_index;
        interference_RRH(:,k)=[];    % 该矩阵所含元素为发射对第k个用户的无用信号的RRH编号
        x_user=USER_matrix(k,1);
        y_user=USER_matrix(k,2);
        H_useful=cell(1,size(useful_RRH,1));
        useful_power=0;
        for i=1:size(useful_RRH,1)
            m=useful_RRH(i,1);
            x_rrh_useful=RRH_matrix(m,1);
            y_rrh_useful=RRH_matrix(m,2);
            H_useful{1,i}=channel(x_rrh_useful,y_rrh_useful,x_user,y_user,Nt,Nr);
            H=H_useful{1,i};
            useful_power=useful_power+Ns*H'*(V*V')*H;
        end
        H_interference=cell(1,size(interference_RRH,1)*size(interference_RRH,2));
        interference_power=0;
        for j=1:size(interference_RRH,2)
            for i=1:size(interference_RRH,1)
                m=interference_RRH(i,j);
                x_rrh_interference=RRH_matrix(m,1);
                y_rrh_interference=RRH_matrix(m,2);
                H_interference{1,i+size(interference_RRH,1)*(j-1)}=channel(x_rrh_interference,y_rrh_interference,x_user,y_user,Nt,Nr);
                H=H_interference{1,i+size(interference_RRH,1)*(j-1)};
                interference_power=interference_power+Ns*H'*(V*V')*H;
            end
        end
        SINR=useful_power/(N0+interference_power);
        R(1,k)=real(log2(1+SINR));
    end
    R_total(loop_time,:)=R;
end
R_average=mean(R_total,1);
end