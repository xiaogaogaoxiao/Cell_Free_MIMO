% 尝试用粒子群算法写一个出来，不管每个用户有没有最小容量限制
% 先尝试使用粒子群背包问题解决单用户的天线选择
% 用背包粒子群算法解单用户的天线选择，最后的结果应该跟基于最近的相同

% 为什么要写基于距离的粒子群天线选择方案
% 之前写的基于信道容量的方案，无法克服信道随机性带来的问题，容易引起全局最优产生变化
% 在distance_channel文件中，生成了log(H'H)与距离的函数关系，呈明显的递减趋势
% 因此，考虑用距离来定义粒子群算法的适应度函数

function [A,y]=PSO_antenna_selection_based_distance_for_single_user(RRH_matrix,USER_matrix,service_number)
rrh=size(RRH_matrix,1);
user=size(USER_matrix,1);
if user~=1
    warning('本函数只能用于单用户情况下')
end

NP = 300;      % 种群个数
D = rrh;        % 决策变量的维度
G = 2000;       % 迭代次数
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
    pbest(k, 1) = func(individual_best(k, :),RRH_matrix,USER_matrix,service_number);
end

% 初始化全局最优
global_best = Inf(1, D);
global_best_fit = inf;
for k=1:NP
    temp = func(individual_best(k, :),RRH_matrix,USER_matrix,service_number);
    if temp < global_best_fit
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
    %         [sv,si]=sort(vs(k,:),2,'descend');
    %         for i=1:size(si,2)
    %             if i<=service_number
    %                 x(k,si(1,i))=1;
    %             else
    %                 x(k,si(1,i))=0;
    %             end
    %         end
    %
    %     end
    % 计算个体历史最优与全局最优
    % 个体历史最优
    for k=1:NP
        old_fitness = func(individual_best(k, :),RRH_matrix,USER_matrix,service_number);
        new_fitness = func(x(k, :),RRH_matrix,USER_matrix,service_number);
        if new_fitness < old_fitness                     % 这里因为距离越小越好，所以用了小于号
            individual_best(k, :) = x(k, :);
            pbest(k, 1) = new_fitness;
        end
    end
    % 全局最优
    for k=1:NP
        temp = func(individual_best(k, :),RRH_matrix,USER_matrix,service_number);
        if temp < global_best_fit                         % 这里因为距离越小越好，所以用了小于号
            global_best = individual_best(k, :);
            global_best_fit = temp;
        end
    end
    y(1,gen)=global_best_fit;
    A(gen,:)=global_best;
end
end

function R=func(x,RRH_matrix,USER_matrix,service_number)   % 适应度函数是选取固定个数的RRH使用户到RRH的总距离最小
% 适应度函数的输入参数
% x: 可行解  二进制串
% service_number:

if sum(x)~=service_number
    R=inf;
    return
end
if sum(x)<1
    R=inf;
    return
end
d=0;
for i=1:size(x,2)
    if x(1,i)==1
        d=d+distance(USER_matrix(1,1),USER_matrix(1,2),RRH_matrix(i,1),RRH_matrix(i,2));
    end
end
R=d;
end

function d=distance(x1,y1,x2,y2)
d=((x1-x2)^2+(y1-y2)^2)^1/2;
end