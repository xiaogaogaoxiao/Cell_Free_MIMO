% script_generate_matrix;
R_baseline=zeros(1,size(yita,2));
loop=100;
for i=1:size(yita,2)
    for j=1:loop
        R_baseline(1,i)=R_baseline(1,i)+baseline(yita(1,i),RRH_matrix,USER_matrix,service_number,power_cell);
    end
end
R_baseline=R_baseline/loop;
plot(yita,R_baseline,'b','LineWidth',1)
hold on
R_distance=zeros(1,size(yita,2));
loop=30;
for i=1:size(yita,2)
    for j=1:loop
        R_distance(1,i)=R_distance(1,i)+PSO_based_distance(yita(1,i),RRH_matrix,USER_matrix,service_number,power_cell);
    end
end
R_distance=R_distance/loop;
plot(yita,R_distance,'g','LineWidth',2)
hold on
R_capacity=zeros(1,size(yita,2));
loop=10;
for i=1:size(yita,2)
    for j=1:loop
        R_capacity(1,i)=R_capacity(1,i)+PSO_based_capacity(yita(1,i),RRH_matrix,USER_matrix,service_number,power_cell);
    end
end
R_capacity=R_capacity/loop;
plot(yita,R_capacity,'r','LineWidth',1)
hold on