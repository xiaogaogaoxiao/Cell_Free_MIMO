% 本程序旨在复现带有上行导频估计的下行信号传输，应该写成一个函数
% 在本程序中需要先建立一个符合要求的导频，供所有用户向所有AP发送
% 因此用户个数和AP个数是输入变量，信道矩阵H是输入变量
% 另外，AP准备向用户发送的信号，以及最后每用户接收到的信号都应该被输出
% 模型中所研究的发送信号、输出信号都是针对某一个AP或者用户的，这里需要以矩阵形式输出全部的信号，不应该以某一个来区分
% 通过对文献
% User-Centric 5G Cellular Networks: Resource Allocation and ComparisonWith the Cell-Free Massive MIMO Approach
% 的学习，决定采用PM信道估计法，并且按照文献给出的下行策略模拟，并生成速率表达
% 注意事项：目前还没有写混合波束成形的部分，因此在下行速率的表达式中，没有添加模拟波束成形矩阵和混合波束成形矩阵
% 当完成波束成形的代码后，需要回来补充

% 先考虑导频问题
% 导频需要具备的特点：导频phik是K*taop维矩阵，其中K是用户的个数，taop是小于信道相干时间样本长度taoc的量
% 如果要一次生成导频矩阵phi，这个矩阵应该是K行taop*K列，属于不同用户的导频横向排列
% phik与其转置的乘积为K*K维单位矩阵，phik为正交矩阵，phik与phij的转置的乘积为0
% 简便起见，设taop为1，这样导频矩阵就变成了K行K列，每列模长为1，每列与自己的转置相乘为单位阵，不同列相乘为0或者尽量接近于0

% 问题来了，这个导频矩阵不会生成...............
% 问了老师，导频设计可以省略掉，信道估计误差用服从某个分布的随机数目来确定

% 假设预编码矩阵Q已经生成完毕，接下来建立User-Centric下行传输模型
Ns=3;Nrf=2;Nt=4;Nr=1;
%发送信号Ns是3*1维向量，2个RF链路，4个发射天线，用于生成后面的随机发送信号
%注意：这里没有考虑波束成形矩阵的约束条件

% 长宽为100的方形区域内有5个用户和20个AP
user=5;rrh=20;range=100;service_number=2;
%该变量为一个用户需要几个RRH为其服务
RRH_matrix=rand(rrh,2)*range;
% 矩阵中的每行即为RRH的坐标
USER_matrix=rand(user,2)*range;
% 矩阵中的每行即为USER的坐标

% 以用户为中心，对每个用户找到距离其最近的几个RRH并生成一个服务矩阵
% 应该对用户距离每个RRH的距离生成一个排好序的序列以解决RRH重复服务的问题
new_RRH_index=antenna_selection(RRH_matrix,USER_matrix,service_number);
% 以上即为在100*100的方格中随机生成的RRH基于距离最近原则服务于随机生成的用户，并且不会存在一个RRH服务于多个用户的情况。

% 尝试生成一个信道细胞元，这个信道细胞元中的每个元素是某一个用户到某一个AP之间的信道，这个细胞元是一个user*rrh的元胞数组
H_cell=cell(user,rrh);
for i=1:user
    for j=1:rrh
        H_cell{i,j}=channel(USER_matrix(i,1),USER_matrix(i,2),RRH_matrix(j,1),RRH_matrix(j,2),Nt,Nr);
    end
end

% 尝试生成一个发送信号细胞元给每个AP分配上。如果这个AP被某User选中，那么在这个位置上就有发送信号，且发给选中它的用户k。如果没有被User选中，则置0.
% 细胞元的每一个元胞数组为Nt*1维矩阵，这个细胞元是1*rrh维
transmit_cell=cell(1,rrh);
for j=1:rrh
    if(ismember(j,new_RRH_index)==1)
        S=rand(Ns,1);
        Vb=randn(Nrf,Ns);
        Vrf=rand(Nt,Nrf);
        x=Vrf*Vb*S;
        % 以上的波束成形矩阵和预编码矩阵都是随便生成的，等完成后面的优化部分需要调用hybrid_beamforming之内的函数来生成这两个矩阵
        % 每个AP会被至多一个User选中，某个AP发送的下行信号为
        % 其中，S是发给第k个用户的向量。
        transmit_cell{1,j}=x;
    else
        transmit_cell{1,j}=zeros(Nt,1);
    end
end

% 下行链路模型中，第k个用户的接收信号应该是来自所有AP的发送信号之和，还需要在其中考虑天线的选择问题：哪些是有用信号，哪些是干扰，哪些是噪声。
% 尝试生成一个接收信号细胞元给每个User分配上
% 其有三部分构成：一是有用信号，即第k个用户所选择的AP给这个用户发送的信号；
% 二是干扰信号，即除了第k个用户所选择的AP之外的所有AP通过信道矩阵发送给第k个用户的信号；
% 三是噪声，服从标准正态分布。这三个量相加，构成接收信号细胞元(cell类型不能用加法，需要对对应元素分别相加，这种做法需要保证细胞元中每个元素可加)
received_useful_cell=cell(1,user);
% 因为下面用到了累加，需要设置初始值，每个细胞元元素均为0，维度是Nr*1维
for i=1:user
    received_useful_cell{1,i}=zeros(Nr,1);
end
% 初始化干扰细胞元和噪声细胞元
received_interference_cell=received_useful_cell;
received_noise_cell=received_useful_cell;
received_cell=received_useful_cell;
% 下一行循环为了生成有用信号。这需要先观察该user选择了哪些AP
for i=1:user
    for j=1:service_number
        received_useful_cell{1,i}=received_useful_cell{1,i}+H_cell{i,new_RRH_index(j,i)}.'*transmit_cell{1,new_RRH_index(j,i)};
    end
end
% 下一行循环为了生成干扰信号。对于每一个user，干扰信号等于所有AP发送的信号，经过该AP与这个user之间的信道被该user接收，减去该user的接收的有用接收信号
for i=1:user
    for j=1:rrh
        received_interference_cell{1,i}=received_interference_cell{1,i}+H_cell{i,j}.'*transmit_cell{1,j};
    end
    received_interference_cell{1,i}=received_interference_cell{1,i}-received_useful_cell{1,i};
end
% 下一行循坏为了生成噪声信号。噪声信号就是服从0,1分布的高斯白噪声
for i=1:user
   received_noise_cell{1,i}=randn(Nr,1); 
end
% 将上述三个变量相加得到最终的接收信号
for i=1:user
   received_cell{1,i}=received_useful_cell{1,i}+received_interference_cell{1,i}+received_noise_cell{1,i};
end



% 2020.10.7 by LiJiaxiang 