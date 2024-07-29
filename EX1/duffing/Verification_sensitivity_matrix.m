clc;clear;close all;
global tf h N_dof N_harm epsilon Tdata x_0 dx_0
%% parameters
w0=5;epsilon=0.8;
x_0=1; dx_0=0;
N_dof=1;N_harm=10;
tf=2*pi/w0;h=2*pi/(w0*1000);Tdata=(0:0.01:10);
% Tdata=0:h:tf;
C0_1=[0.342172594555505,-0.000133005831228,0.000115963368886]';S0_1=[0.1,-0.000712982230221,-0.000108784126763]';
C0_2=[0.122467788388526,0.001563371451967,-0.000646804589735]';S0_2=[0.072813577124383,0.011674778196277,0.000571448020811]';
parameter_a=zeros(N_harm+1,2*N_dof);
parameter_a(1,1)=w0;
parameter_a(2:4,:)=[C0_1,S0_1];   % 6��2��
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
for i=1:2*N_harm*N_dof
    ddt=0.000001;
    %��ʼ������
    parameter_a=zeros(N_harm+1,2*N_dof);
    parameter_a(1,1)=w0;
    parameter_a(2:4,:)=[C0_1,S0_1];   % 6��2��
    
    %% ����ϵ������1������˴������Ĳ��
    sensitivity_parameter_a1=zeros(2*N_harm*N_dof,1);
    sensitivity_parameter_a1(i,1)=1;
    sensitivity_parameter_a1=reshape(sensitivity_parameter_a1,2,N_harm*N_dof);sensitivity_parameter_a1=sensitivity_parameter_a1';
    sensitivity_parameter_a=sensitivity_parameter_a1(1:N_harm,1:2);
    for num_dof=1:N_dof-1
        sensitivity_parameter_a=[sensitivity_parameter_a,sensitivity_parameter_a1(num_dof*N_harm+1:(num_dof+1)*N_harm,1:2)];%da(N_harm+1:2*N_harm,1:2);
    end
    temp_real_w0=zeros(1,2*N_dof);
    sensitivity_parameter_a=[temp_real_w0;sensitivity_parameter_a];
    parameter_a=parameter_a+ddt*sensitivity_parameter_a;% ��������ȥ��
    
    residual1=cal_residual(parameter_a);
    aaaa=(residual1-residual)/ddt;
    figure; 
    plot(Tdata,residual(:,N_dof*(i+1)*2+1),'r-')
    hold on
    plot(Tdata,residual(:,N_dof*(i+1)*2+2),'r-')
    hold on
    plot(Tdata,aaaa(:,1),'k-')
    hold on
    plot(Tdata,aaaa(:,2),'k-')
end
% 
% 
% for i=1
%     % i=8;
%     ddt=0.000001;
%     parameter_a=zeros(N_harm+1,2*N_dof);
%     parameter_a(1,1)=w0+ddt;
%     parameter_a(2:4,:)=[C0_1,S0_1];   % 6��2��
%     residual1=cal_residual(parameter_a);
%     aaaa=(residual1-residual)/ddt;
%     figure; 
%     plot(Tdata,residual(:,3),'r-')
%     hold on  
%     plot(Tdata,residual(:,4),'r-')
%     hold on
%     plot(Tdata,aaaa(:,1),'k-')
%     hold on
%     plot(Tdata,aaaa(:,2),'k-')
% end
















