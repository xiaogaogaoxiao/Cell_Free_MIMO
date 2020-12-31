% 本函数是PSO_antenna_celection_for_multiple_user的简化版，用于绘制yita_R函数用的，是test_yita_R_based_capacity的配套函数
function R=PSO_based_capacity(yita,RRH_matrix,USER_matrix,service_number,power_cell)
NP = 30;                         % 种群个数
G = 300;                          % 迭代次数
c1 = 2;                         % 学习因子
c2 = 2;
w_max = 1.5;                    % 惯性权重
w_min = 0.6;
v_max = 3;                      % 粒子的速度限制
v_min = -3;
rrh=size(RRH_matrix,1);
user=size(USER_matrix,1);
x=zeros(user,rrh,NP);           
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
v = v_min + rand(user,rrh,NP)*(v_max - v_min);
% 速度进行初始化
vs=v;
individual_best = x;            %  每个个体的历史最优
pbest = zeros(NP, 1);           %  个体最优位置对应的适应度值
for k=1:NP
    pbest(k, 1) = fitness_for_multiple_user(yita,individual_best(:,:,k),service_number,power_cell);
end
% 初始化全局最优
global_best = zeros(user,rrh);
global_best_fit = 0;
for k=1:NP
    temp = fitness_for_multiple_user(yita,individual_best(:,:,k),service_number,power_cell);
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
        old_fitness = fitness_for_multiple_user(yita,individual_best(:,:,k),service_number,power_cell);
        new_fitness = fitness_for_multiple_user(yita,x(:,:,k),service_number,power_cell);
        if new_fitness > old_fitness
            individual_best(:,:,k) = x(:,:,k);
            pbest(k, 1) = new_fitness;
        end
    end
    % 全局最优
    for k=1:NP
        temp = fitness_for_multiple_user(yita,individual_best(:,:,k),service_number,power_cell);
        if temp > global_best_fit
            global_best = individual_best(:,:,k);
            global_best_fit = temp;
        end
    end
end
R=norm(global_best_fit);
end

function a=delete_elements(original,pointer)
for i=1:size(pointer,2)
    original(pointer)=0;
end
original(original==0)=[];
a=original;
end