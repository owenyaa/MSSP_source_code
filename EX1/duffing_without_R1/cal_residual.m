function residual=cal_residual(parameter_a)
global epsilon N_dof N_harm Tdata x_0 dx_0
w0=parameter_a(1,1);
Harm_parameter_a=parameter_a(2:end,:);
%% 自由度数目，谐波数目，参数个数2*N_harm*N_dof
% N_dof=2;N_harm=3;
M=1;
C=0;
K=1;%epsilon=0.8;
%% 计算方程残差
x=zeros(N_dof,length(Tdata));dx=zeros(N_dof,length(Tdata));ddx=zeros(N_dof,length(Tdata));
% 有X,Y两个自由度
% for i=1:N_harm   % i=1,3,5
for j=1:N_dof
    for i=1:N_harm   % i=1,3,5
        x(j,:)=x(j,:)+Harm_parameter_a(i,2*j-1)*cos((2*i-1)*w0*Tdata)+Harm_parameter_a(i,2*j)*sin((2*i-1)*w0*Tdata);
        dx(j,:)=dx(j,:)-w0*(2*i-1)*Harm_parameter_a(i,2*j-1)*sin((2*i-1)*w0*Tdata)+w0*(2*i-1)*Harm_parameter_a(i,2*j)*cos((2*i-1)*w0*Tdata);
        ddx(j,:)=ddx(j,:)-(w0*(2*i-1))^2*Harm_parameter_a(i,2*j-1)*cos((2*i-1)*w0*Tdata)-(w0*(2*i-1))^2*Harm_parameter_a(i,2*j)*sin((2*i-1)*w0*Tdata);
    end
end
% residual(1:N_dof,:)=M*ddx+C*dx+K*x+epsilon*x.^3;
% residual(N_dof+1:2*N_dof,:)=1/2*M*dx.^2+1/2*K*x.^2+1/4*epsilon*x.^4-(1/2*M*dx_0^2+1/2*K*x_0^2+1/4*epsilon*x_0^4);

% residual(1:N_dof,:)=M*ddx+C*dx+K*x+epsilon*x.^3;
residual(1:N_dof,:)=1/2*M*dx.^2+1/2*K*x.^2+1/4*epsilon*x.^4-(1/2*M*dx_0^2+1/2*K*x_0^2+1/4*epsilon*x_0^4);
%% 计算频率的灵敏度
x_w0=zeros(N_dof,length(Tdata));dx_w0=zeros(N_dof,length(Tdata));ddx_w0=zeros(N_dof,length(Tdata));
% y_w0=0;dy_w0=0;ddy_w0=0;
for k=1:N_dof
    for i=1:N_harm   % i=1,3,5
        x_w0(k,:)=x_w0(k,:)-Harm_parameter_a(i,2*k-1)*(2*i-1)*Tdata.*sin((2*i-1)*w0*Tdata)+Harm_parameter_a(i,2*k)*(2*i-1)*Tdata.*cos((2*i-1)*w0*Tdata);
        dx_w0(k,:)=dx_w0(k,:)-((2*i-1)*Harm_parameter_a(i,2*k-1)*sin((2*i-1)*w0*Tdata)+w0*Tdata*(2*i-1)^2*Harm_parameter_a(i,2*k-1).*cos((2*i-1)*w0*Tdata))+...
            ((2*i-1)*Harm_parameter_a(i,2*k)*cos((2*i-1)*w0*Tdata)-w0*Tdata*(2*i-1)^2*Harm_parameter_a(i,2*k).*sin((2*i-1)*w0*Tdata));
        ddx_w0(k,:)=ddx_w0(k,:)-(2*w0*(2*i-1)^2*Harm_parameter_a(i,2*k-1)*cos((2*i-1)*w0*Tdata)-w0^2*Tdata*(2*i-1)^3*Harm_parameter_a(i,2*k-1).*sin((2*i-1)*w0*Tdata))-...
            (2*w0*(2*i-1)^2*Harm_parameter_a(i,2*k)*sin((2*i-1)*w0*Tdata)+w0^2*Tdata*(2*i-1)^3*Harm_parameter_a(i,2*k).*cos((2*i-1)*w0*Tdata));
    end
end
% residual(2*N_dof+1:3*N_dof,:)=M*ddx_w0+C*dx_w0+K*x_w0+3*epsilon*(x.^2.*x_w0);
% residual(3*N_dof+1:4*N_dof,:)=M*dx.*dx_w0+K*x.*x_w0+epsilon*x.^3.*x_w0;

% residual(2*N_dof+1:3*N_dof,:)=M*ddx_w0+C*dx_w0+K*x_w0+3*epsilon*(x.^2.*x_w0);
residual(N_dof+1:2*N_dof,:)=M*dx.*dx_w0+K*x.*x_w0+epsilon*x.^3.*x_w0;
%% 计算谐波系数的灵敏度
for i=1:2*N_harm*N_dof
    sensitivity_parameter_a1=zeros(2*N_harm*N_dof,1);
    x_a=zeros(N_dof,length(Tdata));dx_a=zeros(N_dof,length(Tdata));ddx_a=zeros(N_dof,length(Tdata));
    %     y_a=0;dy_a=0;ddy_a=0;
    sensitivity_parameter_a1(i,1)=1;
    sensitivity_parameter_a1=reshape(sensitivity_parameter_a1,2,N_harm*N_dof);sensitivity_parameter_a1=sensitivity_parameter_a1';
    sensitivity_parameter_a=sensitivity_parameter_a1(1:N_harm,1:2);
    for num_dof=1:N_dof-1
        sensitivity_parameter_a=[sensitivity_parameter_a,sensitivity_parameter_a1(num_dof*N_harm+1:(num_dof+1)*N_harm,1:2)];%da(N_harm+1:2*N_harm,1:2);
    end
    for k=1:N_dof
        for j=1:N_harm
            x_a(k,:)=x_a(k,:)+sensitivity_parameter_a(j,2*k-1)*cos((2*j-1)*w0*Tdata)+sensitivity_parameter_a(j,2*k)*sin((2*j-1)*w0*Tdata);
            dx_a(k,:)=dx_a(k,:)-w0*(2*j-1)*sensitivity_parameter_a(j,2*k-1)*sin((2*j-1)*w0*Tdata)+w0*(2*j-1)*sensitivity_parameter_a(j,2*k)*cos((2*j-1)*w0*Tdata);
            ddx_a(k,:)=ddx_a(k,:)-(w0*(2*j-1))^2*sensitivity_parameter_a(j,2*k-1)*cos((2*j-1)*w0*Tdata)-(w0*(2*j-1))^2*sensitivity_parameter_a(j,2*k)*sin((2*j-1)*w0*Tdata);
        end
    end
    %     residual(N_dof*(i+1)*2+1:N_dof*(i+2)*2-1,:)=M*ddx_a+C*dx_a+K*x_a+3*epsilon*(x.^2.*x_a);
    %     residual(N_dof*(i+1)*2+2:N_dof*(i+2)*2,:)=M*dx.*dx_a+K*x.*x_a+epsilon*x.^3.*x_a;
    
    %     residual(N_dof*(i+1)*2+1:N_dof*(i+2)*2-1,:)=M*ddx_a+C*dx_a+K*x_a+3*epsilon*(x.^2.*x_a);
    residual(N_dof*(i+1)*1+1,:)=M*dx.*dx_a+K*x.*x_a+epsilon*x.^3.*x_a;
end

residual=residual';




