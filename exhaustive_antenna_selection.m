% 分布式天线选择的最优穷举搜索法
% 仿照antenna_selection来写

% 原来的antenna_selection的输入为RRH_matrix,USER_matrix,service_number
% 输出为服务RRH矩阵编号
% 若输出矩阵的第一列为1,3,5,7,9，则说明服务于第一个用户的RRH编号为第1,3,5,7,9号
% 今天写的穷举算法也尽量保持这个输入和输出

% 共有R个rrh，K个用户，每个用户选择N个天线
% 复杂度为C(R,NK)*[C(NK,N)*C(N(k-1),N)*C(N(k-2),N)*...*C(2N,N)*C(N,N)]/A(k,k)

% 穷举法可以换个思路，生成一个R行K列的矩阵，每行至多有一个1，其余全是0，不考虑用户最小下行速率约束
% 用户选择天线的数目，即每列1的数目。
% 复杂度为C(R,N)C(R-N,N)...C(R-(K-1)N,N)，共有K个C相乘

function A=exhaustive_antenna_selection(RRH_matrix.USER_matrix.service_number)
N0=10^(-143/10)/1000;
Ns=1;Nt=4;Nr=1;Nrf=2;
rrh=size(RRH_matrix,1);
user=size(USER_matrix,1);
A=zeros(service_number,user);
s1=nchoosek(1:rrh,service_number);
for i=1:size(s1,1)
    % 需要解决k次（k是不确定的）循环嵌套的问题。
    % 怎么在matlab里写？
end
end