%此程序加入了能量约束，位移和速度初值约束
clear;clc;close all;
global Tdata parameter_a N_dof N_harm M K a b N_w0 k1 k2 k12 k22 m1 m2 %index_global
% global N_min alpha temp_zeros
%% different initial values
w1=0.1587*2*pi;w2=0.8133*2*pi;
% w1=0.2319*2*pi;w2=1.416*2*pi;
N_w0=2;
% m1=1;m2=1.1;k1=1.5;k2=3;k12=1.3;k22=2;a=5;b=10;
m1=1;m2=1;k1=1;k2=1;k12=1;k22=1;a=5;b=1;
M=[m1,0;0,m2];K=[k1+k12,-k12;-k12,k2+k12];
N_dof=2;N_harm=13;N_dof1=7;

% 其他参数
Tdata=0:0.05:200;
%% 第一行存储频率，后面存储谐波系数,每个自由度两列
C0_1=[3,2]';S0_1=[0,0]';
C0_2=[3,-2]';S0_2=[0,0.0]';
% C0_1=[1,4]';S0_1=[2,2]';
% C0_2=[2,8]';S0_2=[-5,0.5]';
parameter_a=zeros(N_harm+1,2*N_dof);
parameter_a(1,:)=[w1,w2,0,0];
parameter_a(2:3,:)=[C0_1,S0_1,C0_2,S0_2];   % 6行2列
ini_parameter_a=parameter_a;%b
iteration=length(Tdata);
% to indicate where a is in the parametric space
parameter_a_judge=@(parameter_a)(abs(parameter_a)<10);
% load simple_fre_data.mat; % load observed data --- Tdata and Xdata
gammaT=1.414;rhob=0.5; % parameter for trust-region algorithm
parameter_a_record=parameter_a; % To record the values of parameters during iteration, and for each iteration, a is recorded in a single row of a_record
TR_record=[];  % recording the parameters during trust region
%% 本算例采用速度和加速度响应识别，且三个自由度噪声等级不同，线加速度1%，角加速度2%和5%
%% Response sensitivity iteration
Nmax=1000;   % maximum number for response sensitivity iteration
Ntr=20;      % maximum number for trust region iteration
%% response sensitivity Solution by ode45
% NT=length(residual(:,1));
for iii=1:Nmax
    % compute response and response sensitivity for each incremental
    Etol=1e-10;  % Relative error tolerance for convergence of the RS algorithm
    %%
    residual_iden=cal_residual(parameter_a);
    %% SSS为位移响应灵敏度矩阵，第一列和第二列为残差响应，频率灵敏度从第三列开始
    % 计算的灵敏度矩阵中包含S_11，但是这里要去掉
    SSS=reshape(residual_iden(:,N_dof1+1:2*N_dof1),N_dof1*length(Tdata),1);%w0
    for i=1:2*N_harm*N_dof-1+N_w0
        SSS=[SSS,reshape(residual_iden(:,N_dof1*(i+1)+1:N_dof1*(i+2)),N_dof1*length(Tdata),1)];
    end
    
    dR=-[residual_iden(:,1);residual_iden(:,2);residual_iden(:,3);residual_iden(:,4);residual_iden(:,5);residual_iden(:,6);residual_iden(:,7)];
    
    [U,s,V]=csvd(SSS);
    lambda_inverse=l_curve(U,s,dR);
    atemp=parameter_a(1:N_harm+1,1:2*N_dof);
    % trust-region algorithm
    for trust=1:Ntr
        %% 计算的da排序，为w0,C_11,S_11(缺,不参与迭代),C_12,S_12,按照参数矩阵的形式排列da
        real_da=tikhonov(U,s,V,dR,lambda_inverse);
%         temp_da=zeros(2*N_dof*N_harm,1);
        temp_real_w0=zeros(1,2*N_dof);% 参数矩阵的第一行
        temp_real_w0(1,1)=real_da(1,1);temp_real_w0(1,2)=real_da(2,1);
        %         temp_da(1,1)=real_da(3,1);temp_da(2,1)=0;
        %         temp_da(3:end,1)=real_da(4:end,1);
        da=reshape(real_da(3:end,1),2,N_dof*N_harm);da=da';
        sensitivity_parameter_da=da(1:N_harm,1:2);
        for num_dof=1:N_dof-1
            sensitivity_parameter_da=[sensitivity_parameter_da,da(num_dof*N_harm+1:(num_dof+1)*N_harm,1:2)];%da(N_harm+1:2*N_harm,1:2);
        end
        sensitivity_parameter_da=[temp_real_w0;sensitivity_parameter_da];
        %%
        if ~parameter_a_judge(atemp+sensitivity_parameter_da)      % if updated a is not out of the parametric space, then, lambda should be increased until updated a is in ...
            lambda_inverse=lambda_inverse*gammaT; %  update of lambda
            continue;
        end
        %% 计算的da排序，为w0,C_11,S_11(缺,不参与迭代),C_12,S_12,按照参数矩阵的形式排列da
        real_da=tikhonov(U,s,V,dR,lambda_inverse);
        ini_da=real_da;
%         temp_da=zeros(2*N_dof*N_harm,1);
        temp_real_w0=zeros(1,2*N_dof);% 参数矩阵的第一行
        temp_real_w0(1,1)=real_da(1,1);temp_real_w0(1,2)=real_da(2,1);
        %         temp_da(1,1)=real_da(3,1);temp_da(2,1)=0;
        %         temp_da(3:end,1)=real_da(4:end,1);
        da=reshape(real_da(3:end,1),2,N_dof*N_harm);da=da';
        sensitivity_parameter_da=da(1:N_harm,1:2);
        for num_dof=1:N_dof-1
            sensitivity_parameter_da=[sensitivity_parameter_da,da(num_dof*N_harm+1:(num_dof+1)*N_harm,1:2)];%da(N_harm+1:2*N_harm,1:2);
        end
        sensitivity_parameter_da=[temp_real_w0;sensitivity_parameter_da];
        
        %%
        parameter_a=atemp+sensitivity_parameter_da;
        %% 用新的parameter_a=atemp+da计算响应
        residual_da=cal_residual(parameter_a);

%         dRtemp=-reshape(residual_da(:,1:N_dof),N_dof*length(Tdata),1);
        dRtemp=-[residual_da(:,1);residual_da(:,2);residual_da(:,3);residual_da(:,4);residual_da(:,5);residual_da(:,6);residual_da(:,7)];
        LdR=SSS*ini_da-dR;
        rhos=(dR'*dR-dRtemp'*dRtemp)/(dR'*dR-LdR'*LdR);  % agreement indicator
        if rhos>=rhob
            break;
        end
        lambda_inverse=lambda_inverse*gammaT;
    end
    tolt=norm(da)/norm(parameter_a);
    parameter_a_record=[parameter_a_record,parameter_a];
    TR_record=[TR_record;lambda_inverse];
    parameter_a
    if tolt<=Etol
        break;
    end
    every_a(iii).parameter_a=parameter_a;
    iii
end

residua=cal_residual(parameter_a);

% figure;
% plot(Tdata,residua(:,1),'r-','LineWidth',1);
% hold on;
% plot(Tdata,residua(:,2),'k-','LineWidth',1);
% hold on;
% plot(Tdata,residua(:,3),'b-','LineWidth',1);
% % h1=legend('$$h$$');
% % set(h1,'Interpreter','latex','FontSize',15);
% set(gca,'FontName','Times New Roman','FontSize',15,'LineWidth',1.5);

Tdata=0:0.01:2000;
Harm_parameter_a=parameter_a(2:end,:);
%% 计算两个基频率的组合，注意，频率组合要去掉负数频率
fundamental_w=parameter_a(1,1:N_w0);fundamental_w=[fundamental_w,fundamental_w];
index_global=[1,0;0,1;1,-2;1,2;0,3;4,-1;1,-4;1,4;0,5;6,-1;1,-6;1,6;0,7];
% index_global=[1,0;0,1;0,3;0,5;0,7];
index_global=[index_global,index_global];
for i=1:N_dof
    vector_w(:,i)=index_global(:,2*i-1:2*i)*fundamental_w(1,2*i-1:2*i)';
end
%% 计算方程残差
x=zeros(N_dof,length(Tdata));dx=zeros(N_dof,length(Tdata));ddx=zeros(N_dof,length(Tdata));
% 有三个自由度
for j=1:N_dof
    for i=1:N_harm   % i=1,3,5
        x(j,:)=x(j,:)+Harm_parameter_a(i,2*j-1)*cos(vector_w(i,j)*Tdata)+Harm_parameter_a(i,2*j)*sin(vector_w(i,j)*Tdata);
        dx(j,:)=dx(j,:)-vector_w(i,j)*Harm_parameter_a(i,2*j-1)*sin(vector_w(i,j)*Tdata)+vector_w(i,j)*Harm_parameter_a(i,2*j)*cos(vector_w(i,j)*Tdata);
        ddx(j,:)=ddx(j,:)-(vector_w(i,j))^2*Harm_parameter_a(i,2*j-1)*cos(vector_w(i,j)*Tdata)-(vector_w(i,j))^2*Harm_parameter_a(i,2*j)*sin(vector_w(i,j)*Tdata);
    end
end
figure;
plot(Tdata,x(1,:),'k--','LineWidth',1);
hold on;
plot(Tdata,x(2,:),'r--','LineWidth',1);
% h1=legend('$$h$$');
% set(h1,'Interpreter','latex','FontSize',15);
set(gca,'FontName','Times New Roman','FontSize',15,'LineWidth',1.5);

% figure;
% plot(Tdata,dx(1,:),'k-','LineWidth',1);
% hold on;
% plot(Tdata,dx(2,:),'r-','LineWidth',1);
% % h1=legend('$$h$$');
% % set(h1,'Interpreter','latex','FontSize',15);
% set(gca,'FontName','Times New Roman','FontSize',15,'LineWidth',1.5);



ini_parameter_a
