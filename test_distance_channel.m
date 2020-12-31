% 本程序旨在测试信道矩阵H与距离的关系
% 信道矩阵H为Nt*Nr维，H'*H为Nr*Nr维，距离为1维
% 跟模型设计一样，Nt=4,Nr=1
 Nt=4;Nr=1;
 x=0.5:0.5:200;
 p=0.5:0.5:200;
 base=p;
 N0=10^(-143/10)/1000;
 for i=1:size(x,2)
     temp=0;
     for j=1:200
         H=channel(0,0,x(i),0,Nt,Nr);
         temp=temp+H'*H;
     end
     p(i)=temp/200;
     p(i)=log(p(i));
     base(i)=log(i^(-4));
 end
 plot(x,p,'g')
 hold on
 plot(x,base,'m')
 hold on
 title("距离和信道矩阵之间的关系");
 xlabel("收发天线之间的距离(m)")
 ylabel("信道矩阵转置和其本身的乘积取对数(log(H'*H))");
 
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