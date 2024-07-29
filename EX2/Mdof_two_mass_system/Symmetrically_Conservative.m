% clear;clc;close all;
global A k22
m1=1;m2=1;k1=1;k2=1;k12=1;k22=1;a=5;b=1;
% m1=1;m2=1.1;k1=1.5;k2=3;k12=1.3;k22=2;a=5;b=10;
M=[m1,0;0,m2];K=[k1+k12,-k12;-k12,k2+k12];
A=[zeros(2),eye(2);-M\K,zeros(2)];

t=0:0.01:2000;
% ini_x=[x(1,1);x(2,1);0;0];
ini_x=[a;b;0;0];
options=odeset('RelTol',1e-10,'AbsTol',1e-10);
[tt,num]=ode45('sub_Symmetrically_Conservative',t,ini_x,options);
figure;
plot(t,num(:,1),'r-','LineWidth',1);
hold on;
plot(t,num(:,2),'k-','LineWidth',1);
% h1=legend('$$h$$');
% set(h1,'Interpreter','latex','FontSize',15);
set(gca,'FontName','Times New Roman','FontSize',15,'LineWidth',1.5);

figure;
plot(t,[x(1,:)]'-num(:,1),'r-','LineWidth',1);
hold on;
plot(t,[x(2,:)]'-num(:,2),'k-','LineWidth',1);
% h1=legend('$$h$$');
% set(h1,'Interpreter','latex','FontSize',15);
set(gca,'FontName','Times New Roman','FontSize',15,'LineWidth',1.5);




% fs=100;%����Ƶ��
% % ����Ƶ����ʱ����֮��Ĺ�ϵ�� fs=1/dt
% % ��������������ǣ�����Ƶ��Ҫ�����ź�Ƶ�ʵ������� 
% N=2^16;  %��������2^17
% % N�������㣬����FFT֮�󣬾Ϳ��Եõ�N�����FFT�����Ϊ�˷������FFT���㣬ͨ��Nȡ2�������η���
% % Ҫ��ȷ��xHz������Ҫ��������Ϊ1/x����źţ�����FFT��
% % Ҫ���Ƶ�ʷֱ��ʣ�����Ҫ���Ӳ�������
% n=0:N-1;
% t=n/fs;  % dt=1/fs ��ʾʱ����   fs=1/dt
% y=fft(num(1:end,1),N);  % ����fft�任
% % �������Ƶ��ΪFs���ź�Ƶ��F����������ΪN����ôFFT֮��������һ��ΪN��ĸ�����
% % ÿһ����Ͷ�Ӧ��һ��Ƶ�ʵ㡣������ģֵ�����Ǹ�Ƶ��ֵ�µķ������ԡ�
% % y % ���y����fft֮��Ľ����
% m=abs(y(1:N/2))*2/N; % ���źŵ���ʵ��ֵ
% f=2*pi*n*fs/N;  % ע�� ��2*pi �Ͳ���2*pi������%m=log10(m);
% figure;
% plot(f(1:N/2)./(2*pi),m(1:N/2),'k-','LineWidth',1);
% set(gca,'FontName','Times New Roman','FontSize',15,'LineWidth',1.5);
