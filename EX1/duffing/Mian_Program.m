clear;clc;close all;
global tf h Tdata parameter_a N_dof N_harm epsilon x_0 dx_0
%% different initial values
% case 1
N_dof=1;N_harm=20;
N_dof1=2;
% w0=1.8;epsilon=8;
% x_0=1; dx_0=0;
% w0=18;epsilon=2;
% x_0=15; dx_0=0;

w0=7;epsilon=80;
x_0=1; dx_0=0;



% ��������
tf=2*pi/w0;h=2*pi/(w0*1000);Tdata=(0:h:tf);
%% ��һ�д洢Ƶ�ʣ�����洢г��ϵ��,ÿ�����ɶ�����
% parameter_a=[w0,  0,    0,   0;
%             C_11,S_11,C_21,S_21;
%             C_12,S_12,C_22,S_22;
%             C_13,S_13,C_23,S_23;...];
parameter_a=zeros(N_harm+1,2*N_dof);
parameter_a(1,1)=w0;
parameter_a(2,:)=[-0.9,0.5];% һ��г����ֵ
ini_parameter_a=parameter_a;%b
iteration=length(Tdata);
% to indicate where a is in the parametric space
parameter_a_judge=@(parameter_a)(abs(parameter_a)<10);
% load simple_fre_data.mat; % load observed data --- Tdata and Xdata
gammaT=1.414;rhob=0.5; % parameter for trust-region algorithm
parameter_a_record=parameter_a; % To record the values of parameters during iteration, and for each iteration, a is recorded in a single row of a_record
TR_record=[];  % recording the parameters during trust region
%% �����������ٶȺͼ��ٶ���Ӧʶ�����������ɶ������ȼ���ͬ���߼��ٶ�1%���Ǽ��ٶ�2%��5%
%% Response sensitivity iteration
Nmax=1000;   % maximum number for response sensitivity iteration
Ntr=20;      % maximum number for trust region iteration
%% response sensitivity Solution by ode45
% NT=length(residual(:,1));
for iii=1:Nmax
    %% ��������w0,�˴�Ϊһ�������ھ���ѡȡ1K����
    w0=parameter_a(1,1);
    tf=2*pi/w0;h=2*pi/(w0*1000);Tdata=(0:h:tf);
    % compute response and response sensitivity for each incremental
    Etol=1e-10;  % Relative error tolerance for convergence of the RS algorithm
    %%
    residual_iden=cal_residual(parameter_a);
    %% SSSΪλ����Ӧ�����Ⱦ��󣬵�һ�к͵ڶ���Ϊ�в���Ӧ��Ƶ�������ȴӵ����п�ʼ
    % ����������Ⱦ����а���S_11����������Ҫȥ��
    temp_SSS=reshape(residual_iden(:,N_dof1+1:2*N_dof1),N_dof1*length(Tdata),1);%w0
    for i=1:2*N_harm*N_dof
        temp_SSS=[temp_SSS,reshape(residual_iden(:,N_dof1*(i+1)+1:N_dof1*(i+2)),N_dof1*length(Tdata),1)];
    end
    SSS(:,1:2)=temp_SSS(:,1:2);SSS(:,3:2*N_harm*N_dof)=temp_SSS(:,4:end);
%      SSS=temp_SSS;
%     dR=-reshape(residual_iden(:,1:N_dof),N_dof*length(Tdata),1);
    dR=-[residual_iden(:,1);residual_iden(:,2)];
    [U,s,V]=csvd(SSS);
    lambda_inverse=l_curve(U,s,dR);
    atemp=parameter_a;
    % trust-region algorithm
    for trust=1:Ntr
        %% �����da����Ϊw0,C_11,S_11(ȱ,���������),C_12,S_12,���ղ����������ʽ����da
        real_da=tikhonov(U,s,V,dR,lambda_inverse);
        temp_real_w0=zeros(1,2*N_dof);% ��������ĵ�һ��
        temp_real_w0(1,1)=real_da(1,1);
        real_da(1,1)=real_da(2,1);real_da(2,1)=0;
        da=reshape(real_da,2,N_dof*N_harm);da=da';
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
        %% �����da����Ϊw0,C_11,S_11(ȱ,���������),C_12,S_12,���ղ����������ʽ����da
        real_da=tikhonov(U,s,V,dR,lambda_inverse);
        ini_da=real_da;
        temp_real_w0=zeros(1,2*N_dof);% ��������ĵ�һ��
        temp_real_w0(1,1)=real_da(1,1);real_da(1,1)=real_da(2,1);real_da(2,1)=0;
        da=reshape(real_da,2,N_dof*N_harm);da=da';
        sensitivity_parameter_da=da(1:N_harm,1:2);
        for num_dof=1:N_dof-1
            sensitivity_parameter_da=[sensitivity_parameter_da,da(num_dof*N_harm+1:(num_dof+1)*N_harm,1:2)];%da(N_harm+1:2*N_harm,1:2);
        end
        sensitivity_parameter_da=[temp_real_w0;sensitivity_parameter_da];
        %%
        parameter_a=atemp+sensitivity_parameter_da;
        %% ���¼���Tdata
        w0=parameter_a(1,1);
        tf=2*pi/w0;h=2*pi/(w0*1000);Tdata=(0:h:tf);
        %% ���µ�parameter_a=atemp+da������Ӧ
        residual_da=cal_residual(parameter_a);

%         dRtemp=-reshape(residual_da(:,1:N_dof),N_dof*length(Tdata),1);
        dRtemp=-[residual_da(:,1);residual_da(:,2)];
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

figure;
plot(Tdata,residua(:,2),'r-','LineWidth',1);
h1=legend('$$h$$');
set(h1,'Interpreter','latex','FontSize',15);
set(gca,'FontName','Times New Roman','FontSize',15,'LineWidth',1.5);

Tdata=0:0.001:200;
w=parameter_a(1,1);
Harm_parameter_a=parameter_a(2:end,:);
%% ���㷽�̲в�
x=zeros(N_dof,length(Tdata));dx=zeros(N_dof,length(Tdata));ddx=zeros(N_dof,length(Tdata));
% ��X,Y�������ɶ�
% for i=1:N_harm   % i=1,3,5
for j=1:N_dof
    for i=1:N_harm   % i=1,3,5
        x(j,:)=x(j,:)+Harm_parameter_a(i,2*j-1)*cos((2*i-1)*w0*Tdata)+Harm_parameter_a(i,2*j)*sin((2*i-1)*w0*Tdata);
        dx(j,:)=dx(j,:)-w0*(2*i-1)*Harm_parameter_a(i,2*j-1)*sin((2*i-1)*w0*Tdata)+w0*(2*i-1)*Harm_parameter_a(i,2*j)*cos((2*i-1)*w0*Tdata);
        ddx(j,:)=ddx(j,:)-(w0*(2*i-1))^2*Harm_parameter_a(i,2*j-1)*cos((2*i-1)*w0*Tdata)-(w0*(2*i-1))^2*Harm_parameter_a(i,2*j)*sin((2*i-1)*w0*Tdata);
    end
end
figure;
plot(Tdata,x(1,:),'k-','LineWidth',1);
h1=legend('$$h$$');
set(h1,'Interpreter','latex','FontSize',15);
set(gca,'FontName','Times New Roman','FontSize',15,'LineWidth',1.5);

% figure;
% plot(x(1,:),dx(1,:),'k-','LineWidth',1);
% h1=legend('$$h$$');
% set(h1,'Interpreter','latex','FontSize',15);
% set(gca,'FontName','Times New Roman','FontSize',15,'LineWidth',1.5);



ini_parameter_a
