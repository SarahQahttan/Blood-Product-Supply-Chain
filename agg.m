function [AG]=agg(C)
[n m]=size(C);
for i=1:n
    s1=1;
    s2=1;
    s3=1;
    for j=1:m
        s1=s1*(1-C{i,j}(1))^(1/n);
        s2=s2*(C{i,j}(2))^(1/n);
        s3=s3*(C{i,j}(3))^(1/n);
    end
    AG{i}=[1-s1 s2 s3];
end