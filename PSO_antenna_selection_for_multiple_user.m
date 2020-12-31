% 根据信道容量，尝试写多用户的天线选择粒子群算法
% 跟基于距离的多用户天线选择粒子群算法不同的是，现在规定RRH可以服务rf个用户

function [A,y,power_cell]=PSO_antenna_selection_for_multiple_user(RRH_matrix,USER_matrix,service_number,Nrf)
rrh=size(RRH_matrix,1);
user=size(USER_matrix,1);
LOS_matrix=zeros(user,rrh);     % 该矩阵用于记录第i个用户到第j个RRH之间有无直射路径
Nt=4;Nr=1;                      % 用于生成channel_cell
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

NP = 100;                         % 种群个数
G = 1000;                          % 迭代次数
c1 = 2;                         % 学习因子
c2 = 2;
w_max = 1.5;                    % 惯性权重
w_min = 0.6;
v_max = 3;                      % 粒子的速度限制
v_min = -3;
x=zeros(user,rrh,NP);           % 初始化种群个体 每个个体是一个user行rrh列的矩阵，共有NP个，这里采用三维空间堆叠方案
% 需要对x每一层的个体进行处理，使个体满足行和列的限制条件
A=zeros(user,rrh,G);
y=zeros(1,G);
set=eye(user);
for i=1:NP
    a=1:rrh;
    for j=1:user
        b=randperm(size(a,2));
        b=b(1:service_number);
        for k=1:size(b,2)
            x(:,a(1,b(1,k)),i)=set(:,j);
        end
        a=delete_elements(a,b);
    end
end                             % 这个for循环，实际上只能生成每个rrh服务于1个用户的情况
% 按道理来说，生成每个rrh服务于Nrf个用户的情况应该有利于收敛
% 可是目前并没有想到这种个体应该具备什么一般特征
% 其实，可行集set应该是user*(nchoosek(user,Nrf)+user)维矩阵
% 前nchoosek(user,Nrf)是每列有Nrf个1，后user是每列有1个1
% 但是却不能再用service_number这个来筛选位置
% 因为前几列在相同的行位置会有1
% 如何解决？
v = v_min + rand(user,rrh,NP)*(v_max - v_min);
% 速度进行初始化
vs=v;
individual_best = x;            %  每个个体的历史最优
pbest = zeros(NP, 1);           %  个体最优位置对应的适应度值
for k=1:NP
    pbest(k, 1) = fitness_for_multiple_user(1,individual_best(:,:,k),service_number,power_cell);
end
% 初始化全局最优
global_best = zeros(user,rrh);
global_best_fit = 0;
for k=1:NP
    temp = fitness_for_multiple_user(1,individual_best(:,:,k),service_number,power_cell);
    if temp > global_best_fit
        global_best = individual_best(:,:,k);
        global_best_fit = temp;
    end
end
% 进行迭代
for gen = 1:G
    w = w_max - (w_max-w_min) * gen / G;
    for k=1:NP
        % 更新速度
        v(:,:,k) = w * v(:,:,k) + c1 * rand() * (individual_best(:,:,k) - x(:,:,k)) + c2 * rand() * (global_best - x(:,:,k));
        % 边界条件处理
        for i=1:user
            for j=1:rrh
                if v(i,j,k)>v_max
                    v(i,j,k)=v_max;
                end
                if v(i,j,k)<v_min
                    v(i,j,k)=v_min;
                end
            end
        end
        % 使用sigmoid函数对速度进行映射
        vs(:,:,k) = abs(1./(1+exp(-v(:,:,k))));
        % 更新粒子的位置
        for i=1:user
            for j=1:rrh
                if vs(i,j,k)>rand()
                    x(i,j,k) = 1;
                else
                    x(i,j,k) = 0;
                end
            end
        end
    end
    % 计算个体历史最优与全局最优
    % 个体历史最优
    for k=1:NP
        old_fitness = fitness_for_multiple_user(1,individual_best(:,:,k),service_number,power_cell);
        new_fitness = fitness_for_multiple_user(1,x(:,:,k),service_number,power_cell);
        if new_fitness > old_fitness
            individual_best(:,:,k) = x(:,:,k);
            pbest(k, 1) = new_fitness;
        end
    end
    % 全局最优
    for k=1:NP
        temp = fitness_for_multiple_user(1,individual_best(:,:,k),service_number,power_cell);
        if temp > global_best_fit
            global_best = individual_best(:,:,k);
            global_best_fit = temp;
        end
    end
    y(1,gen)=norm(global_best_fit);
    A(:,:,gen)=global_best;
end
end

function d=distance(x1,y1,x2,y2)
d=((x1-x2)^2+(y1-y2)^2)^1/2;
end

function a=delete_elements(original,pointer)
for i=1:size(pointer,2)
    original(pointer)=0;
end
original(original==0)=[];
a=original;
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