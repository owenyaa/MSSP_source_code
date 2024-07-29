% clear;clc;close all;
global A epsilon
m=1;k=1;epsilon=80;
A=[0,1;-m\k,0];
t=0:0.01:400;
ini_x=[x(1,1);dx(1,1)];
% ini_x=[1;0];
options=odeset('RelTol',1e-10,'AbsTol',1e-10);
[tt,num]=ode45('sub_duffing',t,ini_x,options);
figure;
plot(t,num(:,1),'r-')

% ¸Ä½øµÄFFt
    
% t=t(900000:end);
% 
% xx=num(900000:end,1);
% 
% dt=t(2)-t(1);
% data=xx;
% N=length(data);
% %tt=0:dt:(N-1)*dt;tt=tt';
% N_fft=2^16;
% Y=fft(data,N_fft);
% Pyy=Y.*conj(Y)/N_fft;
% f=1/dt*(0:N_fft/2)/N_fft;
% figure
% plot(f,(Pyy(1:(N_fft/2+1))))
% title('Frequency content of y')
% xlabel('frequency (Hz)')
