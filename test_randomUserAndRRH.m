% function average_capacity=randomUserAndRRH(limit)
% for limit=0:0.1:0.8
% 100*100的方格，RRH和USER随机生成位置
user=20;
rrh=200;
range=100;
service_number=1;
loop=10;
power_matrix=zeros(user,loop);
capacity_matrix=zeros(user,loop);
sum_capacity_matrix=zeros(1,loop);
for time=1:loop
    %该变量为一个用户需要几个RRH为其服务
    RRH_matrix=rand(rrh,2)*range;
    % 矩阵中的每行即为RRH的坐标
    USER_matrix=rand(user,2)*range;
    % 矩阵中的每行即为USER的坐标
    
    % 以用户为中心，对每个用户找到距离其最近的几个RRH并生成一个服务矩阵
    % 应该对用户距离每个RRH的距离生成一个排好序的序列以解决RRH重复服务的问题
    distance_square=zeros(rrh,user);
    RRH_index=distance_square;
    for i=1:user
        for j=1:rrh
            distance_square(j,i)=(USER_matrix(i,1)-RRH_matrix(j,1))^2+(USER_matrix(i,2)-RRH_matrix(j,2))^2;
            RRH_index(j,i)=j;
        end
    end
    
    for k=1:user
        for i=1:rrh-1
            for j=1:rrh-i
                if(distance_square(j,k)>distance_square(j+1,k))
                    [distance_square(j,k),distance_square(j+1,k)]=swap(distance_square(j,k),distance_square(j+1,k));
                    [RRH_index(j,k),RRH_index(j+1,k)]=swap(RRH_index(j,k),RRH_index(j+1,k));
                end
            end
        end
    end
    % 已经完成以用户为中心，计算RRH到用户之间的距离，目前也已经实现了RRH序号索引，确保对距离RRH排序过程中，RRH索引也能跟随变化，这样能进行之后的基于索引的RRH重复服务的检测。
    
    % 设置服务于用户的RRH为可调参数，这样可以在今后进行扩展为多RRH服务于用户。
    % 如果存在RRH重复使用，对distance_square矩阵和RRH_index进行变换。如果RRH_index在一行内出现重复，则对应位置距离远的选择下一行的index进行服务；如果index中被使用了，则选择下一行的index进行服务。
    already_used=[];
    new_RRH_index=RRH_index;
    new_distance_square=distance_square;
    for k=1:user*service_number
        for i=1:user
            if(ismember(new_RRH_index(k,i),already_used)==0)
                already_used(end+1)=new_RRH_index(k,i);
            else
                while(ismember(new_RRH_index(k,i),already_used)==1)
                    temp=new_RRH_index(:,i);
                    temp(k)=[];
                    temp(end+1)=0;
                    new_RRH_index(:,i)=temp;
                    
                    temp_d=new_distance_square(:,i);
                    temp_d(k)=[];
                    temp_d(end+1)=0;
                    new_distance_square(:,i)=temp_d;
                    if(ismember(0,already_used)==1)
                        break
                    end
                end
                already_used(end+1)=new_RRH_index(k,i);
            end
        end
    end
    new_RRH_index=new_RRH_index(1:service_number,1:user);
    new_distance_square=new_distance_square(1:service_number,1:user);
    % new_RRH_index:无重复的服务RRH索引组成的矩阵，支持多RRH服务于一个用户。该矩阵每列即为服务于该用户的RRH编号。
    % 存在的问题：单纯是使用用户的先后顺序处理RRH重复服务的问题，一开始想使用距离来处理冲突，未能实现，当很多RRH服务于用户的时候，性能可能会下降。少量RRH服务用户时暂不处理。
    
    % 以上即为在100*100的方格中随机生成的RRH基于距离最近原则服务于随机生成的用户，并且不会存在一个RRH服务于多个用户的情况。
    
    % 接下来需要寻找合适的大小尺度衰落来建立容量计算公式。大尺度衰落选择反比于1/d^2，小尺度选择服从瑞利衰落raylrnd(0.5),二者相加得到a。
    N0=1;
    nt=2;nr=1; % 两个发射天线，一个接受天线
    a=1./new_distance_square;
    n=user;
    ones_vector=ones(service_number,n);
    sum_p=20;  % 分配给系统的总的功率
    k=0;
    while(1==1)
        k=k+1;
        a=a+raylrnd(0.5);   % a是一个行向量
        cvx_begin quiet
        variable p(n,service_number)
        maximize ones_vector*(((a.^2).').*p)    % 存在的问题：当服务天线大于1时，目标函数不是标量。
        subject to
        ones_vector*p<=sum_p
        min(log(1+N0^(-1)*nt*nr*((a.^2).').*p))>=0.2
        cvx_end
        if(isnan(p(1,1))~=1)
            break
        end
        if(k==10)
            disp('The capacity limiter is too high.');
            break
        end
    end  % 其实存在随机的小尺度衰落使得凸优化问题无解的情况，这里的代码意思是，当出现无解，重新生成小尺度衰落再解一遍凸优化问题
    C_matrix=log2(1+N0^(-1)*nt*nr*((a.^2).').*p);
    % 以上考虑了大小尺度衰落，辅以凸优化进行功率分配。容量计算公式为C=log（1+P*a^2*nt*nr/N0）
    % 总功率定位20w，每用户进行通信的最低容量限度为0.8bit/s/Hz
    % 现在输出为每个用户分配的功率，每用户通信容量与系统总容量
    
    if(isnan(p(1))~=1)
        capacity_sum=0;
        for i=1:user
            %disp(['The power and capacity of user ' num2str(i) ' are ' num2str(p(i)) 'W and ' num2str(C_matrix(i)) 'bit/s/Hz.']);
            capacity_sum=capacity_sum+log2(1+N0^(-1)*nt*nr*(a(i)^2)*p(i));
        end
        %disp(['The total capacity of the system is ' num2str(capacity_sum) 'bit/s/Hz.']);
    end
    power_matrix(:,time)=p;
    capacity_matrix(:,time)=C_matrix;
    sum_capacity_matrix(time)=capacity_sum;
end

% 现在需要进行循环，取得性能指标的统计量。
% 从随机生成用户和基站位置开始循环，但需要记录每一次的系统输出。
% 每次输出为功率分配阵p，容量阵C_matrix。二者均为列向量.
% 循环已经添加完成

% 接下来需要算系统容量平均值来衡量指标。
average_capacity=mean(sum_capacity_matrix)
% end

% 告一段落。本程序的凸优化基本上没什么用，但是编写了一个不科学的非最优的天线选择算法。
% 每个用户挑选最近的一个RRH为其服务，且一个RRH最多同时服务于一个用户。
% 2020年10月3日还是，准备做新的模型。