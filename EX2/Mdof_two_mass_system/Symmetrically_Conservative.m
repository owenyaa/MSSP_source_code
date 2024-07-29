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




% fs=100;%采样频率
% % 采样频率与时间间隔之间的关系： fs=1/dt
% % 采样定理告诉我们，采样频率要大于信号频率的两倍。 
% N=2^16;  %采样点数2^17
% % N个采样点，经过FFT之后，就可以得到N个点的FFT结果。为了方便进行FFT运算，通常N取2的整数次方。
% % 要精确到xHz，则需要采样长度为1/x秒的信号，并做FFT。
% % 要提高频率分辨率，就需要增加采样点数
% n=0:N-1;
% t=n/fs;  % dt=1/fs 表示时间间隔   fs=1/dt
% y=fft(num(1:end,1),N);  % 进行fft变换
% % 假设采样频率为Fs，信号频率F，采样点数为N。那么FFT之后结果就是一个为N点的复数。
% % 每一个点就对应着一个频率点。这个点的模值，就是该频率值下的幅度特性。
% % y % 输出y看看fft之后的结果。
% m=abs(y(1:N/2))*2/N; % 求信号的真实幅值
% f=2*pi*n*fs/N;  % 注意 乘2*pi 和不乘2*pi的区别%m=log10(m);
% figure;
% plot(f(1:N/2)./(2*pi),m(1:N/2),'k-','LineWidth',1);
% set(gca,'FontName','Times New Roman','FontSize',15,'LineWidth',1.5);
