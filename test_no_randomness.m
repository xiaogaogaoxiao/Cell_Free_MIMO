% 本程序是用来测试无信道随机的三种预编码方案的测试程序
x=1:1:20; % 横坐标为对某一个用户同时服务的RRH数目
y=zeros(1,20);
w=zeros(1,20);
v=zeros(1,20);
user=3;rrh=60;range=100;
USER_matrix=rand(user,2)*range;
for i=1:size(x,2)
    [a,b,c]=no_channel_randomness(x(i),USER_matrix,rrh,range);
    y(i)=mean(a);
    w(i)=mean(b);
    v(i)=mean(c);
end
plot(x,y,'r-o')
hold on
plot(x,w,'g-o')
hold on
plot(x,v,'b-o')
title('多个用户在不同波束成形情况下的速率');
xlabel('服务于该用户的天线数目');
ylabel('下行可达速率');
legend('全随机波束成形','考虑无约束的最佳模拟预编码','提取相位信息的恒模约束模拟预编码','location','northeast');
