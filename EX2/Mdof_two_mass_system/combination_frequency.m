% clear;clc;close all;
w1=0.2319;w2=1.416;
N_dof=1;
index_global=[1 0;0 1];
for N=3:2:9
    temp=[N 0];
    for n=N-1:-1:1
        m=N-n;
        temp=[temp;n m;n -m];
    end
    temp=[temp;0 N];
    index_global=[index_global;temp];
end
fundamental_w=[w1,w2];%fundamental_w=[fundamental_w,fundamental_w,fundamental_w];
for i=1:N_dof
    vector_w(:,i)=index_global(:,2*i-1:2*i)*fundamental_w(1,2*i-1:2*i)';
end
