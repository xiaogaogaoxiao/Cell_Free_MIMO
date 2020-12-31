function position_ploter(matrix,RRH_matrix,USER_matrix)
rrh=size(matrix,2);user=size(matrix,1);
rrh_x=zeros(1,rrh);
rrh_y=rrh_x;
user_x=zeros(1,user);
user_y=user_x;
for i=1:rrh
    rrh_x(1,i)=RRH_matrix(i,1);
    rrh_y(1,i)=RRH_matrix(i,2);
end
for i=1:user
    user_x(1,i)=USER_matrix(i,1);
    user_y(1,i)=USER_matrix(i,2);
end
plot(rrh_x,rrh_y,'b*')
hold on
plot(user_x,user_y,'ro')
hold on
for i=1:user
    for j=1:rrh
        if matrix(i,j)==1
            plot([user_x(1,i),rrh_x(1,j)],[user_y(1,i),rrh_y(1,j)],'g-')
        end
    end
end