% script_generate_matrix;
R_distance=zeros(1,size(yita,2));
loop=10;
for i=1:size(yita,2)
    for j=1:loop
        R_distance(1,i)=R_distance(1,i)+PSO_based_distance(yita(1,i),RRH_matrix,USER_matrix,service_number,power_cell);
    end
end
R_distance=R_distance/loop;
plot(R_distance,'g')
hold on