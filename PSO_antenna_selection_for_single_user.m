% 尝试用粒子群算法写一个出来，不管每个用户有没有最小容量限制
% 先尝试使用粒子群背包问题解决单用户的天线选择
% 用背包粒子群算法解单用户的天线选择，最后的结果应该跟基于最近的相同

function [A,y]=PSO_antenna_selection_for_single_user(RRH_matrix,USER_matrix,service_number)
tic
rrh=size(RRH_matrix,1);
user=size(USER_matrix,1);
LOS_matrix=zeros(user,rrh);     % 该矩阵用于记录第i个用户到第j个RRH之间有无直射路径
Nt=4;Nr=1;                      % 用于生成channel_cell
tic
for i=1:user
    for j=1:rrh
        d=distance(USER_matrix(i,1),USER_matrix(i,2),RRH_matrix(j,1),RRH_matrix(j,2));
        LOS_matrix(i,j)=rand(1)<min(20/d,1)*(1-exp(-d/39))+exp(-d/39);
    end
end
disp(['生成LOS_matrix所用时间：',num2str(toc)]);
% LOS_matrix
% 上面的语句可以选择性输出。因为直射路径是否存在可能会对rrh的选择有影响，导致实验组和对照组（基于距离选）不一致
channel_cell=cell(user,rrh);    % 该cell用于存储第i个用户到第j个RRH之间的平均信道，即平均掉小尺度的信道
loop=20;                        % 该参数用于更改平均计算所用的循环次数
tic
for i=1:user
    for j=1:rrh
        for l=1:loop
            if size(channel_cell{i,j},1)==0||size(channel_cell{i,j},2)==0   % 新生成的cell行数和列数都为0，不进行重新赋值操作会无法与new_channel矩阵相加
                channel_cell{i,j}=zeros(Nt,Nr);
            end
            channel_cell{i,j}=channel_cell{i,j}+new_channel(USER_matrix(i,1),USER_matrix(i,2),RRH_matrix(j,1),RRH_matrix(j,2),Nt,Nr,LOS_matrix(i,j));
        end
        channel_cell{i,j}=channel_cell{i,j}/loop;
    end
end
disp(['生成channel_cell所用时间：',num2str(toc)]);
if user~=1
    warning('本函数只能用于单用户情况下')
end
tic
NP = 100;      % 种群个数
D = rrh;        % 决策变量的维度
G = 1000;       % 迭代次数
c1 = 2;      % 学习因子
c2 = 2;
w_max = 1.5;   % 惯性权重
w_min = 0.6;
% t=c1+c2;
% w = 2/(2-t-(t^2-4*t)^0.5)^0.5;
v_max = 3;     % 粒子的速度限制
v_min = -3;
% 初始化种群个体
y=zeros(1,G);
A=zeros(G,D);
% x = (rand(NP,D)>0.5);                        % 产生均匀分布的二进制串  randn产生的是符合正态分布的随机数
% 需要考虑service_number个1元素
x=zeros(NP,D);                                 % 经过下面的for处理，可以选择每个粒子中的x都有service_number个位置置1
for i=1:NP
    a=randperm(D);
    a=a(1:service_number);
    for j=1:service_number
        x(i,a(j))=1;
    end
end
v = v_min + rand(NP,D)*(v_max - v_min);        % 速度进行初始化
vs=v;
% 初始化个体最优
individual_best = x;       %  每个个体的历史最优
pbest = zeros(NP, 1);      %  个体最优位置对应的适应度值
for k=1:NP
    pbest(k, 1) = func(individual_best(k, :),service_number,channel_cell);
end

% 初始化全局最优
global_best = zeros(1, D);
global_best_fit = eps;
for k=1:NP
    temp = func(individual_best(k, :),service_number,channel_cell);
    if temp > global_best_fit
        global_best = individual_best(k, :);
        global_best_fit = temp;
    end
end

% 进行迭代
for gen = 1:G
    w = w_max - (w_max-w_min) * gen / G;
    % w=0.9;
    for k=1:NP
        % 更新速度
        v(k, :) = w * v(k, :) + c1 * rand() * (individual_best(k, :) - x(k, :)) + c2 * rand() * (global_best - x(k, :));
        % 边界条件处理    % 边界吸收
        for t=1:D
            if v(k, t) > v_max
                v(k, t) = v_max;
            end
            if v(k, t) < v_min
                v(k, t) = v_min;
            end
        end
        % 使用sigmoid函数对速度进行映射
        vs(k, :) = abs(1./(1+exp(-v(k, :))));
        % 更新粒子的位置
        for t=1:D
            if vs(k, t)>rand()
                x(k, t) = 1;
            else
                x(k, t) = 0;
            end
        end
    end
    % 计算个体历史最优与全局最优
    % 个体历史最优
    for k=1:NP
        old_fitness = func(individual_best(k, :),service_number,channel_cell);
        new_fitness = func(x(k, :),service_number,channel_cell);
        if new_fitness > old_fitness
            individual_best(k, :) = x(k, :);
            pbest(k, 1) = new_fitness;
        end
    end
    % 全局最优
    for k=1:NP
        temp = func(individual_best(k, :),service_number,channel_cell);
        if temp > global_best_fit
            global_best = individual_best(k, :);
            global_best_fit = temp;
        end
    end
    y(1,gen)=global_best_fit;
    A(gen,:)=global_best;
end
disp(['循环迭代所用时间：',num2str(toc)]);
disp(['基于容量的函数的总时间：',num2str(toc)]);
end

function R=func(x,service_number,channel_cell)   % 适应度函数应该是选取固定个数的RRH使单用户下行信道容量最大
% 适应度函数的输入参数
% x: 可行解  二进制串
% service_number:
% RRH_index和USER_index用于算信道容量
if sum(x)~=service_number
    symbol=0;
else
    symbol=1;
end
N0=10^(-143/10)/1000;Nrf=2;useful_power=0;
for i=1:size(x,2)
    if x(1,i)==1
        H=channel_cell{1,i};
        [Vrf,Vb]=precoder(H,Nrf,3);
        useful_power=useful_power+H'*Vrf*(Vb*Vb')*Vrf'*H;
    end
end
R=abs(symbol*log(1+useful_power/N0));
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

function d=distance(x1,y1,x2,y2)
d=((x1-x2)^2+(y1-y2)^2)^1/2;
end