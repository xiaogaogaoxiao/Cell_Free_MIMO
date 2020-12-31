x=1:1:20; % 横坐标为对某一个用户同时服务的RRH数目
y=zeros(1,20);
user=1;rrh=20;range=100;
USER_matrix=[range/2,range/2];
% for i=1:size(x,2)
%     a=no_beamforming(x(i),USER_matrix);
%     y(i)=a(1);% 这里选择对第一个用户进行分析
% end
% plot(x,y,'y-o')
% hold on
for i=1:size(x,2)
    a=capacity(x(i),USER_matrix);
    y(i)=a(1);% 这里选择对第一个用户进行分析
end
plot(x,y,'r-o')
hold on
for i=1:size(x,2)
    b=no_constraint_analog_precoder(x(i),USER_matrix);
    y(i)=b(1);
end
plot(x,y,'g-o');
hold on
for i=1:size(x,2)
    c=analog_precoder(x(i),USER_matrix);
    y(i)=c(1);
end
plot(x,y,'b-o');
title('单个用户在不同波束成形情况下的速率');
xlabel('服务于该用户的天线数目');
ylabel('下行可达速率');
legend('全随机波束成形','无约束的最佳模拟预编码','提取相位信息的恒模约束模拟预编码','location','northwest');
