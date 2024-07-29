function residual=cal_residualnew(parameter_a)
global r1 r4 N_dof N_harm Tdata 
global M C K f
w0=parameter_a(1,1);
Harm_parameter_a=parameter_a(2:end,:);
%% 计算方程残差
x=zeros(N_dof,length(Tdata));dx=zeros(N_dof,length(Tdata));ddx=zeros(N_dof,length(Tdata));
for j=1:N_dof
    for i=1:N_harm   % i=1,3,5
        x(j,:)=x(j,:)+Harm_parameter_a(i,2*j-1)*cos(i*w0*Tdata)+Harm_parameter_a(i,2*j)*sin(i*w0*Tdata);
        dx(j,:)=dx(j,:)-w0*i*Harm_parameter_a(i,2*j-1)*sin(i*w0*Tdata)+w0*i*Harm_parameter_a(i,2*j)*cos(i*w0*Tdata);
        ddx(j,:)=ddx(j,:)-(w0*i)^2*Harm_parameter_a(i,2*j-1)*cos(i*w0*Tdata)-(w0*i)^2*Harm_parameter_a(i,2*j)*sin(i*w0*Tdata);
    end
end
non_matrix=zeros(N_dof,N_dof);
non_matrix(1,1)=r1;non_matrix(4,4)=r4;
f_matrix=zeros(N_dof,length(Tdata));f_matrix(3,:)=f*sin(w0*Tdata);
residual(1:N_dof,:)=M*ddx+C*dx+K*x+non_matrix*x.^3-f_matrix;
% residual(1:N_dof,:)=M*ddx+C*dx+K*x+[r,0;0,beta]*x.^3;
%% 计算频率的灵敏度
x_w0=zeros(N_dof,length(Tdata));dx_w0=zeros(N_dof,length(Tdata));ddx_w0=zeros(N_dof,length(Tdata));
% y_w0=0;dy_w0=0;ddy_w0=0;
for k=1:N_dof
    for i=1:N_harm   % i=1,3,5
        x_w0(k,:)=x_w0(k,:)-Harm_parameter_a(i,2*k-1)*i*Tdata.*sin(i*w0*Tdata)+Harm_parameter_a(i,2*k)*i*Tdata.*cos(i*w0*Tdata);
        dx_w0(k,:)=dx_w0(k,:)-(i*Harm_parameter_a(i,2*k-1)*sin(i*w0*Tdata)+w0*Tdata*i^2*Harm_parameter_a(i,2*k-1).*cos(i*w0*Tdata))+...
            (i*Harm_parameter_a(i,2*k)*cos(i*w0*Tdata)-w0*Tdata*i^2*Harm_parameter_a(i,2*k).*sin(i*w0*Tdata));
        ddx_w0(k,:)=ddx_w0(k,:)-(2*w0*i^2*Harm_parameter_a(i,2*k-1)*cos(i*w0*Tdata)-w0^2*Tdata*i^3*Harm_parameter_a(i,2*k-1).*sin(i*w0*Tdata))-...
            (2*w0*i^2*Harm_parameter_a(i,2*k)*sin(i*w0*Tdata)+w0^2*Tdata*i^3*Harm_parameter_a(i,2*k).*cos(i*w0*Tdata));
    end
end
residual(N_dof+1:2*N_dof,:)=M*ddx_w0+C*dx_w0+K*x_w0+3*non_matrix*(x.^2.*x_w0);
%% 计算谐波系数的灵敏度
for i=1:2*N_harm*N_dof
    sensitivity_parameter_a1=zeros(2*N_harm*N_dof,1);
    x_a=zeros(N_dof,length(Tdata));dx_a=zeros(N_dof,length(Tdata));ddx_a=zeros(N_dof,length(Tdata));
    sensitivity_parameter_a1(i,1)=1;
    sensitivity_parameter_a1=reshape(sensitivity_parameter_a1,2,N_harm*N_dof);sensitivity_parameter_a1=sensitivity_parameter_a1';
    sensitivity_parameter_a=sensitivity_parameter_a1(1:N_harm,1:2);
    for num_dof=1:N_dof-1
        sensitivity_parameter_a=[sensitivity_parameter_a,sensitivity_parameter_a1(num_dof*N_harm+1:(num_dof+1)*N_harm,1:2)];%da(N_harm+1:2*N_harm,1:2);
    end
    for k=1:N_dof
        for j=1:N_harm
            x_a(k,:)=x_a(k,:)+sensitivity_parameter_a(j,2*k-1)*cos(j*w0*Tdata)+sensitivity_parameter_a(j,2*k)*sin(j*w0*Tdata);
            dx_a(k,:)=dx_a(k,:)-w0*j*sensitivity_parameter_a(j,2*k-1)*sin(j*w0*Tdata)+w0*j*sensitivity_parameter_a(j,2*k)*cos(j*w0*Tdata);
            ddx_a(k,:)=ddx_a(k,:)-(w0*j)^2*sensitivity_parameter_a(j,2*k-1)*cos(j*w0*Tdata)-(w0*j)^2*sensitivity_parameter_a(j,2*k)*sin(j*w0*Tdata);
        end
    end
    residual(N_dof*(i+1)+1:N_dof*(i+2),:)=M*ddx_a+C*dx_a+K*x_a+3*non_matrix*(x.^2.*x_a);
end

residual=residual';




