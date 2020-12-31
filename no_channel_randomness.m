% 本程序旨在处理信道带来的随机性波动问题，单纯考虑不同预编码方案带来的问题
% 也就是说，将要生成一个RRH_matrix和USER_matrix并确保他们之间H的大小尺度衰落都相同，单纯研究不同预编码方案的下行速率
% 输入为test函数中规定的USER_matrix和service_number，这里再加两个，rrh和range，这样就不用来这里修改rrh和range了
% 输出为三种不同预编码方案经过多次循环取均值之后的下行速率。注意输出应该是3个
function [random_precoder,no_constraint,constraint]=no_channel_randomness(service_number,USER_matrix,rrh,range)
Ns=1;Nrf=2;Nt=4;Nr=1;N0=10^(-143/10)/1000;
user=size(USER_matrix,1);
R1=zeros(1,user);
R2=zeros(1,user);
R3=zeros(1,user);
loop=500;
R_total_1=zeros(loop,user);
R_total_2=zeros(loop,user);
R_total_3=zeros(loop,user);
for loop_time=1:loop
    RRH_matrix=rand(rrh,2)*range;
    for k=1:user
        % 研究第k个用户的下行容量,目前考虑单用户情况，这个k只会等于1，循环只进行一次
        % 先计算useful_power，这个量需要知道哪些天线被选中，然后生成被选中的天线到第k个用户之间的信道矩阵
        new_RRH_index=enhanced_antenna_selection(RRH_matrix,USER_matrix,service_number);
        useful_RRH=new_RRH_index(:,k); % 该矩阵所含元素为发射对第k个用户的有用信号的RRH编号
        interference_RRH=new_RRH_index;
        interference_RRH(:,k)=[];    % 该矩阵所含元素为发射对第k个用户的无用信号的RRH编号
        x_user=USER_matrix(k,1);
        y_user=USER_matrix(k,2);
        H_useful=cell(1,size(useful_RRH,1));
        useful_power_1=0;
        useful_power_2=0;
        useful_power_3=0;
        for i=1:size(useful_RRH,1)
            m=useful_RRH(i,1);
            x_rrh_useful=RRH_matrix(m,1);
            y_rrh_useful=RRH_matrix(m,2);
            H_useful{1,i}=channel(x_rrh_useful,y_rrh_useful,x_user,y_user,Nt,Nr);
            H=H_useful{1,i};
            [Vrf_1,Vb_1]=precoder(H,Nrf,1);
            [Vrf_2,Vb_2]=precoder(H,Nrf,2);
            [Vrf_3,Vb_3]=precoder(H,Nrf,3);
            % Nt*Nrf维矩阵，由信道矩阵H的SVD分解后的右酉矩阵的转置的前Nrf列构成
            % 注意，这里输入函数的H是Nt*Nr维，没有转置处理
            V_1=(Vrf_1*Vb_1)/norm(Vrf_1*Vb_1);
            useful_power_1=useful_power_1+Ns*H'*(V_1*V_1')*H;
            useful_power_2=useful_power_2+Ns*H'*Vrf_2*Vb_2*(Vb_2')*Vrf_2'*H;
            useful_power_3=useful_power_3+Ns*H'*Vrf_3*Vb_3*(Vb_3')*Vrf_3'*H;
        end
        H_interference=cell(1,size(interference_RRH,1)*size(interference_RRH,2));
        interference_power_1=0;
        interference_power_2=0;
        interference_power_3=0;
        for j=1:size(interference_RRH,2)
            for i=1:size(interference_RRH,1)
                m=interference_RRH(i,j);
                x_rrh_interference=RRH_matrix(m,1);
                y_rrh_interference=RRH_matrix(m,2);
                H_interference{1,i+size(interference_RRH,1)*(j-1)}=channel(x_rrh_interference,y_rrh_interference,x_user,y_user,Nt,Nr);
                H=H_interference{1,i+size(interference_RRH,1)*(j-1)};
                interference_power_1=interference_power_1+Ns*H'*Vrf_1*Vb_1*(Vb_1')*Vrf_1'*H;
                interference_power_2=interference_power_2+Ns*H'*Vrf_2*Vb_2*(Vb_2')*Vrf_2'*H;
                interference_power_3=interference_power_3+Ns*H'*Vrf_3*Vb_3*(Vb_3')*Vrf_3'*H;
            end
        end
        SINR_1=useful_power_1/(N0+interference_power_1);
        SINR_2=useful_power_2/(N0+interference_power_2);
        SINR_3=useful_power_3/(N0+interference_power_3);
        R1(1,k)=norm(log2(1+SINR_1));
        R2(1,k)=norm(log2(1+SINR_2));
        R3(1,k)=norm(log2(1+SINR_3));
    end
    R_total_1(loop_time,:)=R1;
    R_total_2(loop_time,:)=R2;
    R_total_3(loop_time,:)=R3;
end
random_precoder=mean(R_total_1,1);
no_constraint=mean(R_total_2,1);
constraint=mean(R_total_3,1);
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
    [U,S,V]=svd(H');
    transpose_V=V';
    Vrf=transpose_V(:,1:Nrf);
    Vrf=Vrf/norm(Vrf);
    Vb=(H'*Vrf)'/norm(H'*Vrf,2);
    Vb=Vb/norm(Vb);
elseif type==3 % 有约束的模拟预编码
    [U,S,V]=svd(H');
    transpose_V=V';
    Vrf=angle(transpose_V(:,1:Nrf));
    Vrf=Vrf/norm(Vrf);
    Vb=(H'*Vrf)'/norm(H'*Vrf,2);
    Vb=Vb/norm(Vb);
end
end