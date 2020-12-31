% function average_capacity=randomUserAndRRH(limit)
% for limit=0:0.1:0.8
% 100*100�ķ���RRH��USER�������λ��
user=20;
rrh=200;
range=100;
service_number=1;
loop=10;
power_matrix=zeros(user,loop);
capacity_matrix=zeros(user,loop);
sum_capacity_matrix=zeros(1,loop);
for time=1:loop
    %�ñ���Ϊһ���û���Ҫ����RRHΪ�����
    RRH_matrix=rand(rrh,2)*range;
    % �����е�ÿ�м�ΪRRH������
    USER_matrix=rand(user,2)*range;
    % �����е�ÿ�м�ΪUSER������
    
    % ���û�Ϊ���ģ���ÿ���û��ҵ�����������ļ���RRH������һ���������
    % Ӧ�ö��û�����ÿ��RRH�ľ�������һ���ź���������Խ��RRH�ظ����������
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
    % �Ѿ�������û�Ϊ���ģ�����RRH���û�֮��ľ��룬ĿǰҲ�Ѿ�ʵ����RRH���������ȷ���Ծ���RRH��������У�RRH����Ҳ�ܸ���仯�������ܽ���֮��Ļ���������RRH�ظ�����ļ�⡣
    
    % ���÷������û���RRHΪ�ɵ����������������ڽ�������չΪ��RRH�������û���
    % �������RRH�ظ�ʹ�ã���distance_square�����RRH_index���б任�����RRH_index��һ���ڳ����ظ������Ӧλ�þ���Զ��ѡ����һ�е�index���з������index�б�ʹ���ˣ���ѡ����һ�е�index���з���
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
    % new_RRH_index:���ظ��ķ���RRH������ɵľ���֧�ֶ�RRH������һ���û����þ���ÿ�м�Ϊ�����ڸ��û���RRH��š�
    % ���ڵ����⣺������ʹ���û����Ⱥ�˳����RRH�ظ���������⣬һ��ʼ��ʹ�þ����������ͻ��δ��ʵ�֣����ܶ�RRH�������û���ʱ�����ܿ��ܻ��½�������RRH�����û�ʱ�ݲ�����
    
    % ���ϼ�Ϊ��100*100�ķ�����������ɵ�RRH���ھ������ԭ�������������ɵ��û������Ҳ������һ��RRH�����ڶ���û��������
    
    % ��������ҪѰ�Һ��ʵĴ�С�߶�˥���������������㹫ʽ����߶�˥��ѡ�񷴱���1/d^2��С�߶�ѡ���������˥��raylrnd(0.5),������ӵõ�a��
    N0=1;
    nt=2;nr=1; % �����������ߣ�һ����������
    a=1./new_distance_square;
    n=user;
    ones_vector=ones(service_number,n);
    sum_p=20;  % �����ϵͳ���ܵĹ���
    k=0;
    while(1==1)
        k=k+1;
        a=a+raylrnd(0.5);   % a��һ��������
        cvx_begin quiet
        variable p(n,service_number)
        maximize ones_vector*(((a.^2).').*p)    % ���ڵ����⣺���������ߴ���1ʱ��Ŀ�꺯�����Ǳ�����
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
    end  % ��ʵ���������С�߶�˥��ʹ��͹�Ż������޽�����������Ĵ�����˼�ǣ��������޽⣬��������С�߶�˥���ٽ�һ��͹�Ż�����
    C_matrix=log2(1+N0^(-1)*nt*nr*((a.^2).').*p);
    % ���Ͽ����˴�С�߶�˥�䣬����͹�Ż����й��ʷ��䡣�������㹫ʽΪC=log��1+P*a^2*nt*nr/N0��
    % �ܹ��ʶ�λ20w��ÿ�û�����ͨ�ŵ���������޶�Ϊ0.8bit/s/Hz
    % �������Ϊÿ���û�����Ĺ��ʣ�ÿ�û�ͨ��������ϵͳ������
    
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

% ������Ҫ����ѭ����ȡ������ָ���ͳ������
% ����������û��ͻ�վλ�ÿ�ʼѭ��������Ҫ��¼ÿһ�ε�ϵͳ�����
% ÿ�����Ϊ���ʷ�����p��������C_matrix�����߾�Ϊ������.
% ѭ���Ѿ�������

% ��������Ҫ��ϵͳ����ƽ��ֵ������ָ�ꡣ
average_capacity=mean(sum_capacity_matrix)
% end

% ��һ���䡣�������͹�Ż�������ûʲô�ã����Ǳ�д��һ������ѧ�ķ����ŵ�����ѡ���㷨��
% ÿ���û���ѡ�����һ��RRHΪ�������һ��RRH���ͬʱ������һ���û���
% 2020��10��3�ջ��ǣ�׼�����µ�ģ�͡�