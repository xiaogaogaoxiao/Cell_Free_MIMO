% 根据文献“毫米波MIMO系统的混合预编码研究_张丹”,提取了无约束的模拟波束成形矩阵的相位信息作为有约束的模拟波束成形矩阵
% 单用户情况下，analog_precoder的性能应该介于无约束的模拟波束成形和全随机的预编码中间
% 经过测试也的确是达到了预想的情况
% 多用户情况下，analog_precoder和无约束的模拟预编码的情况不能得出明显规律

function R_average=analog_precoder(service_number,USER_matrix)
% service_number=2;
Ns=1;Nrf=2;Nt=4;Nr=1;N0=10^(-143/10)/1000;
user=size(USER_matrix,1);
rrh=60;range=100;R=zeros(1,user);
loop=10;
R_total=zeros(loop,user);
% USER_matrix=rand(user,2)*range;
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
        useful_power=0;
        for i=1:size(useful_RRH,1)
            m=useful_RRH(i,1);
            x_rrh_useful=RRH_matrix(m,1);
            y_rrh_useful=RRH_matrix(m,2);
            H_useful{1,i}=channel(x_rrh_useful,y_rrh_useful,x_user,y_user,Nt,Nr);
            H=H_useful{1,i};
            Vrf=get_analog_precoder(H,Nrf);
            Vb=get_digital_precoder(H,Vrf);
            % Nt*Nrf维矩阵，由信道矩阵H的SVD分解后的右酉矩阵的转置的前Nrf列构成
            % 注意，这里输入函数的H是Nt*Nr维，没有转置处理
            useful_power=useful_power+Ns*H'*Vrf*Vb*(Vb')*Vrf'*H;
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
                interference_power=interference_power+Ns*H'*Vrf*Vb*(Vb')*Vrf'*H;
            end
        end
        SINR=useful_power/(N0+interference_power);
        R(1,k)=norm(log2(1+SINR));
    end
    R_total(loop_time,:)=R;
end
R_average=mean(R_total,1);
end

function Vrf=get_analog_precoder(H,Nrf)
% 旨在对信道矩阵H进行SVD分解，取右酉矩阵转置的前Ns列组成不受恒模约束的模拟波束成形矩阵
% 信道矩阵需要传入，Nrf需要传入；处理后的模拟波束成形矩阵需要传出
[U,S,V]=svd(H');
transpose_V=V';
Vrf=angle(transpose_V(:,1:Nrf));
Vrf=Vrf/norm(Vrf);
end

function Vb=get_digital_precoder(H,Vrf)
% 本函数旨在生成数字波束成形矩阵，其等于等效矩阵的转置除以等效矩阵的2范数
% 等效矩阵由实际信道H'和模拟波束成形矩阵相乘得到，因此需要传入参数H和Vrf
% 传出参数为数字波束成形矩阵
Vb=(H'*Vrf)'/norm(H'*Vrf,2);
Vb=Vb/(norm(Vb));
end