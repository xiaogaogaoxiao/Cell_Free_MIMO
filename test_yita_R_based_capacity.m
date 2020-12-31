% script_generate_matrix;
R_capacity=zeros(1,size(yita,2));
loop=3;
for i=1:size(yita,2)
    for j=1:loop
        R_capacity(1,i)=R_capacity(1,i)+PSO_based_capacity(yita(1,i),RRH_matrix,USER_matrix,service_number,power_cell);
    end
end
R_capacity=R_capacity/loop;
plot(R_capacity,'r')
hold on
