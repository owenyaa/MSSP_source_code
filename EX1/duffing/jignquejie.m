% clear;
% m=1;k=1;epsilon=0.8;
m=1;k=0.5;epsilon=80;
syms fi
% x_0=1;dx_0=0;
x_0=1;dx_0=0;
mu=0.5*epsilon*x_0^2/(1+epsilon*x_0^2);
% x=0:0.01:10;
% k1=1.5;k2=5;
% y=k1*x+k2.*x.^3;
% plot(x,y)
T=4/sqrt(1+epsilon*x_0^2)*(double(int(1/sqrt(1-mu*(sin(fi))^2),0,2\pi)))

(T-2*pi/(parameter_a(1,1)))/T*100



