% 根据基于距离进行的单用户粒子群算法选择天线的函数，写一个基于距离进行的多用户粒子群天线选择算法
% 目标：每个用户根据距离远近，选择天线service_number根，每个天线只能服务于一个用户
% *被某个用户选中的天线不可再被其他用户选中*
% 选出的选择方案中，用户到天线的总距离应该是最小的

function [A,y]=PSO_antenna_selection_based_distance_for_multiple_user(RRH_matrix,USER_matrix,service_number)
user=size(USER_matrix,1);
rrh=size(RRH_matrix,1);
NP = 300;      % 种群个数
G = 2000;       % 迭代次数
c1 = 2;      % 学习因子
c2 = 2;
w_max = 1.5;   % 惯性权重
w_min = 0.6;
v_max = 3;     % 粒子的速度限制
v_min = -3;
x=zeros(user,rrh,NP); % 初始化种群个体 每个个体是一个user行rrh列的矩阵，共有NP个，这里采用三维空间堆叠方案
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
end
% 经过上面的循环，可以确保生成的种群中的每个个体都满足约束要求
v = v_min + rand(user,rrh,NP)*(v_max - v_min);        % 速度进行初始化
vs=v;
individual_best = x;       %  每个个体的历史最优
pbest = zeros(NP, 1);      %  个体最优位置对应的适应度值
for k=1:NP
    pbest(k, 1) = func(individual_best(:,:,k),RRH_matrix,USER_matrix,service_number);
end
% 初始化全局最优
global_best = Inf(user,rrh);
global_best_fit = inf;
for k=1:NP
    temp = func(individual_best(:,:,k),RRH_matrix,USER_matrix,service_number);
    if temp < global_best_fit
        global_best = individual_best(:,:,k);
        global_best_fit = temp;
    end
end

% 进行迭代
for gen = 1:G
    w = w_max - (w_max-w_min) * gen / G;
    % w=0.9;
    for k=1:NP
        % 更新速度
        v(:,:,k) = w * v(:,:,k) + c1 * rand() * (individual_best(:,:,k) - x(:,:,k)) + c2 * rand() * (global_best - x(:,:,k));
        % 边界条件处理    % 边界吸收
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
        % 多用户情况下，如果将选择矩阵中的每个元素进行粒子位置的更新，将会产生大量的不和规则的个体
        % 这种情况下，粒子将不会跟新的选择有什么关联
        % 导致pbest极大概率保持不变化
        % 因此，决定将vs矩阵更新为以单位阵中的列以及0列为元素的矩阵
        % 提取这个列矩阵中所有变化可能的列矩阵，给他们打乱索引顺序后，再放回列矩阵中
        % 这样能生成新的个体，且新个体满足模型所建立的规则。
        % 但是问题是，这样是不是就跟速度无关了？
        
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
        old_fitness = func(individual_best(:,:,k),RRH_matrix,USER_matrix,service_number);
        new_fitness = func(x(:,:,k),RRH_matrix,USER_matrix,service_number);
        if new_fitness < old_fitness                     % 这里因为距离越小越好，所以用了小于号
            individual_best(:,:,k) = x(:,:,k);
            pbest(k, 1) = new_fitness;
        end
    end
    % 全局最优
    for k=1:NP
        temp = func(individual_best(:,:,k),RRH_matrix,USER_matrix,service_number);
        if temp < global_best_fit                         % 这里因为距离越小越好，所以用了小于号
            global_best = individual_best(:,:,k);
            global_best_fit = temp;
        end
    end
    y(1,gen)=global_best_fit;
    A(:,:,gen)=global_best;
end
end

function R=func(x,RRH_matrix,USER_matrix,service_number)   
% 适应度函数的输入参数
% x: 可行解  是一个user行rrh列的二进制矩阵
% service_number:

% 这里的x将变为矩阵形式
% 第i行第j列为1，则说明第i个用户选择了第j个rrh
% 所以，适应度函数输入的可行解有明显的约束
% 首先，跟单用户情况一样，每行求和应该等于固定值service_number
% 其次，列和求和应该等于1
% 输出R应该所讨论的用户到其所选择的天线之间的距离
if sum(ismember(sum(x,1),1))~=size(USER_matrix,1)*service_number
    R=inf;
    return
end
% 上一个if用于判断，每个天线仅可被选择一次，即将矩阵行与行相加，最后的向量仅有1和0元素
% 但是可能存在每个用户选的天线个数不是service_number的情况

if ~isequal(sum(x,2),service_number*ones(size(USER_matrix,1),1))
    R=inf;
    return
end


% 上一个if用于判断每个用户选择service_number个天线，集将矩阵列与列相加，最后表现为每个元素都为service_number的列矩阵
d=0;
for j=1:size(x,1)
    for i=1:size(x,2)
        if x(j,i)==1
            d=d+distance(USER_matrix(j,1),USER_matrix(j,2),RRH_matrix(i,1),RRH_matrix(i,2));
        end
    end
end
R=d;
end

function a=delete_elements(original,pointer)
for i=1:size(pointer,2)
    original(pointer)=0;
end
original(original==0)=[];
a=original;
end

function d=distance(x1,y1,x2,y2)
d=((x1-x2)^2+(y1-y2)^2)^1/2;
end