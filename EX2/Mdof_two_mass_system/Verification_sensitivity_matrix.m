clc;clear;close all;
global tf h N_dof N_harm M K Tdata a b N_w0 k2 k k1 m1 m2
%% parameters
w1=5;w2=2;
N_w0=2;
m1=1.1;m2=1.2;k=1.3;k1=1.4;k2=1.5;a=5;b=1;
M=[m1,0;0,m2];K=[k+k1,-k1;-k1,k+k1];
N_dof=2;N_harm=5;N_dof1=7;
tf=2*pi/w1;h=2*pi/(w1*1000);
Tdata=(0:0.01:10);
% Tdata=0:h:tf;
C0_1=[0.342172594555505,-0.000133005831228,0.000115963368886]';S0_1=[0.1,-0.000712982230221,-0.000108784126763]';
C0_2=[0.122467788388526,0.001563371451967,-0.000646804589735]';S0_2=[0.072813577124383,0.011674778196277,0.000571448020811]';
parameter_a=zeros(N_harm+1,2*N_dof);
parameter_a(1,:)=[w1,w2,0,0];
parameter_a(2:4,:)=[C0_1,S0_1,C0_2,S0_2];   % 6��2��
%% ����в�
residual=cal_residual(parameter_a);
%% ���Ʋв�����
% figure;
% plot(Tdata,residual(:,1),'r-','LineWidth',1);
% h1=legend('$$h$$');
% set(h1,'Interpreter','latex','FontSize',15);
% set(gca,'FontName','Times New Roman','FontSize',15,'LineWidth',1.5);

%% �˴�Ϊ��������֤���򣬱���Ϊ�ò�ִ��������ȣ�������֤�������Ƿ����
% ����������֮�󣬽�ĳ��������ȥһ��С�������²�������һ��
% �������ý�������ٳ���С����Ϊ�����ȣ���֤��ֵ������Ⱥ�ֱ�Ӽ���������������Ƿ��غ�
% ����Ҫע�⣬��ֽ���϶��ǶԵ�(�����Ĳ������)�����ܳ������ԭ��������ݣ�x_cal������
% ������ڶ�Ӧ����ֽ��(parameter_a1)�е�λ��(x1(1,:))���ٶȣ����ٶȶ�Ӧ��ԭ����(parameter_a)�е�����������x_cal(2,:)
% 
% for i=1:2*N_harm*N_dof
% % for i=1:3
%     ddt=0.000001;
%     %��ʼ������
%     parameter_a=zeros(N_harm+1,2*N_dof);
%     parameter_a(1,:)=[w1,w2,0,0];
%     parameter_a(2:4,:)=[C0_1,S0_1,C0_2,S0_2];   % 6��2��
%     
%     % ����ϵ������1������˴������Ĳ��
%     sensitivity_parameter_a1=zeros(2*N_harm*N_dof,1);
%     sensitivity_parameter_a1(i,1)=1;
%     sensitivity_parameter_a1=reshape(sensitivity_parameter_a1,2,N_harm*N_dof);sensitivity_parameter_a1=sensitivity_parameter_a1';
%     sensitivity_parameter_a=sensitivity_parameter_a1(1:N_harm,1:2);
%     for num_dof=1:N_dof-1
%         sensitivity_parameter_a=[sensitivity_parameter_a,sensitivity_parameter_a1(num_dof*N_harm+1:(num_dof+1)*N_harm,1:2)];%da(N_harm+1:2*N_harm,1:2);
%     end
%     temp_real_w0=zeros(1,2*N_dof);
%     sensitivity_parameter_a=[temp_real_w0;sensitivity_parameter_a];
%     parameter_a=parameter_a+ddt*sensitivity_parameter_a;% ��������ȥ��
%     
%     residual1=cal_residual(parameter_a);
%     aaaa=(residual1-residual)/ddt;
%     figure;
%     plot(Tdata,residual(:,N_dof1*(i+N_w0)+1),'r-')
%     hold on
%     plot(Tdata,residual(:,N_dof1*(i+N_w0)+2),'r-')
%     hold on
%     plot(Tdata,residual(:,N_dof1*(i+N_w0)+3),'r-')
%     hold on
%     plot(Tdata,residual(:,N_dof1*(i+N_w0)+4),'r-')
%     hold on
%     plot(Tdata,residual(:,N_dof1*(i+N_w0)+5),'r-')
%     hold on
%     plot(Tdata,residual(:,N_dof1*(i+N_w0)+6),'r-')
%     hold on
%     plot(Tdata,residual(:,N_dof1*(i+N_w0)+7),'r-')
%     hold on
%     
%     plot(Tdata,aaaa(:,1),'k-')
%     hold on
%     plot(Tdata,aaaa(:,2),'k-')
%     hold on
%     plot(Tdata,aaaa(:,3),'k-')
%     hold on
%     plot(Tdata,aaaa(:,4),'k-')
%     hold on
%     plot(Tdata,aaaa(:,5),'k-')
%     hold on
%     plot(Tdata,aaaa(:,6),'k-')
%     hold on
%     plot(Tdata,aaaa(:,7),'k-')
% end

% 
% 
for i=1:2
    % i=8;
    ddt=0.0001;
    w0=zeros(1,2*N_dof);
    w0(1,1:N_w0)=[w1,w2];
    w0(1,i)=w0(1,i)+ddt;
    parameter_a(1,:)=w0;
    parameter_a(2:4,:)=[C0_1,S0_1,C0_2,S0_2];   % 6��2��
    
    residual1=cal_residual(parameter_a);
    aaaa=(residual1-residual)/ddt;
    figure; 
    plot(Tdata,residual(:,i*N_dof1+1),'r-')
    hold on
    plot(Tdata,residual(:,i*N_dof1+2),'r-')
    hold on
    plot(Tdata,residual(:,i*N_dof1+3),'r-')
    hold on
    plot(Tdata,residual(:,i*N_dof1+4),'r-')
    hold on
    plot(Tdata,residual(:,i*N_dof1+5),'r-')
    hold on
    plot(Tdata,residual(:,i*N_dof1+6),'r-')
    hold on
    plot(Tdata,residual(:,i*N_dof1+7),'r-')
    hold on
    
    plot(Tdata,aaaa(:,1),'k-')
    hold on
    plot(Tdata,aaaa(:,2),'k-')
    hold on
    plot(Tdata,aaaa(:,3),'k-')
    hold on
    plot(Tdata,aaaa(:,4),'k-')
    hold on
    plot(Tdata,aaaa(:,5),'k-')
    hold on
    plot(Tdata,aaaa(:,6),'k-')
    hold on
    plot(Tdata,aaaa(:,7),'k-')
    
end
















