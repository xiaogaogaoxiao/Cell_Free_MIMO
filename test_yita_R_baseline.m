% script_generate_matrix;
R_baseline=zeros(1,size(yita,2));
loop=100;
for i=1:size(yita,2)
    for j=1:loop
        R_baseline(1,i)=R_baseline(1,i)+baseline(yita(1,i),RRH_matrix,USER_matrix,service_number,power_cell);
    end
end
R_baseline=R_baseline/loop;
plot(R_baseline,'b')
hold on