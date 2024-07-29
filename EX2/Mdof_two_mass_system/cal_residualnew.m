function residual=cal_residualnew(parameter_a)
global M K N_dof N_harm Tdata a b N_w0 k1 k2 k12 k22 m1 m2 
% global index_global
N_dof1=5;
Harm_parameter_a=parameter_a(2:end,:);
%% ����������Ƶ�ʵ���ϣ�ע�⣬Ƶ�����Ҫȥ������Ƶ��
fundamental_w=parameter_a(1,1:N_w0);fundamental_w=[fundamental_w,fundamental_w];
index_global=[1,0;0,1;1,-2;1,2;0,3;4,-1;1,-4;1,4;0,5;6,-1;1,-6;1,6;0,7];
index_global=[index_global,index_global];
for i=1:N_dof
    vector_w(:,i)=index_global(:,2*i-1:2*i)*fundamental_w(1,2*i-1:2*i)';
end
%% ���㷽�̲в�
x=zeros(N_dof,length(Tdata));dx=zeros(N_dof,length(Tdata));ddx=zeros(N_dof,length(Tdata));
% ���������ɶ�
for j=1:N_dof
    for i=1:N_harm   % i=1,3,5
        x(j,:)=x(j,:)+Harm_parameter_a(i,2*j-1)*cos(vector_w(i,j)*Tdata)+Harm_parameter_a(i,2*j)*sin(vector_w(i,j)*Tdata);
        dx(j,:)=dx(j,:)-vector_w(i,j)*Harm_parameter_a(i,2*j-1)*sin(vector_w(i,j)*Tdata)+vector_w(i,j)*Harm_parameter_a(i,2*j)*cos(vector_w(i,j)*Tdata);
        ddx(j,:)=ddx(j,:)-(vector_w(i,j))^2*Harm_parameter_a(i,2*j-1)*cos(vector_w(i,j)*Tdata)-(vector_w(i,j))^2*Harm_parameter_a(i,2*j)*sin(vector_w(i,j)*Tdata);
    end
end
residual(1:N_dof,:)=M*ddx+K*x+k22*[(x(1,:)-x(2,:)).^3;(x(2,:)-x(1,:)).^3];
residual(N_dof+1,:)=1/2*m1*dx(1,:).^2+1/2*m2*dx(2,:).^2+1/2*k1*x(1,:).^2+1/2*k2*x(2,:).^2+1/2*k12*(x(1,:)-x(2,:)).^2+1/4*k22*(x(1,:)-x(2,:)).^4-...
(1/2*k1*a^2+1/2*k2*b^2+1/2*k12*(a-b)^2+1/4*k22*(a-b)^4);
residual(N_dof+2,:)=sum(Harm_parameter_a(:,1))-a;
residual(N_dof+3,:)=sum(Harm_parameter_a(:,3))-b;
%% ����Ƶ�ʵ�������  ֻ��Ҫ�����Ƶ�ľͿ���
for j=1:N_w0
    x_w0=zeros(N_dof,length(Tdata));dx_w0=zeros(N_dof,length(Tdata));ddx_w0=zeros(N_dof,length(Tdata));
    for ij=1:N_dof
        for i=1:N_harm   % i=1,3,5
            x_w0(ij,:)=x_w0(ij,:)-index_global(i,j)*Harm_parameter_a(i,2*ij-1)*Tdata.*sin(vector_w(i,ij)*Tdata)+index_global(i,j)*Harm_parameter_a(i,2*ij)*Tdata.*cos(vector_w(i,ij)*Tdata);
            dx_w0(ij,:)=dx_w0(ij,:)-(index_global(i,j)*Harm_parameter_a(i,2*ij-1)*sin(vector_w(i,ij)*Tdata)+vector_w(i,ij)*index_global(i,j)*Tdata*Harm_parameter_a(i,2*ij-1).*cos(vector_w(i,ij)*Tdata))+...
                (index_global(i,j)*Harm_parameter_a(i,2*ij)*cos(vector_w(i,ij)*Tdata)-vector_w(i,ij)*index_global(i,j)*Tdata*Harm_parameter_a(i,2*ij).*sin(vector_w(i,ij)*Tdata));
            ddx_w0(ij,:)=ddx_w0(ij,:)-(2*vector_w(i,ij)*index_global(i,j)*Harm_parameter_a(i,2*ij-1)*cos(vector_w(i,ij)*Tdata)-vector_w(i,ij)^2*index_global(i,j)*Tdata*Harm_parameter_a(i,2*ij-1).*sin(vector_w(i,ij)*Tdata))-...
                (2*vector_w(i,ij)*index_global(i,j)*Harm_parameter_a(i,2*ij)*sin(vector_w(i,ij)*Tdata)+vector_w(i,ij)^2*index_global(i,j)*Tdata*Harm_parameter_a(i,2*ij).*cos(vector_w(i,ij)*Tdata));
        end
    end 
    residual(N_dof1*j+1:(1+j)*N_dof1-3,:)=M*ddx_w0+K*x_w0+3*k22*[(x(1,:)-x(2,:)).^2.*(x_w0(1,:)-x_w0(2,:));(x(2,:)-x(1,:)).^2.*(x_w0(2,:)-x_w0(1,:))];
    residual((1+j)*N_dof1-2,:)=m1*dx(1,:).*dx_w0(1,:)+m2*dx(2,:).*dx_w0(2,:)+k1*x(1,:).*x_w0(1,:)+k2*x(2,:).*x_w0(2,:)+k12*(x(1,:)-x(2,:)).*(x_w0(1,:)-x_w0(2,:))+k22*(x(1,:)-x(2,:)).^3.*(x_w0(1,:)-x_w0(2,:));
    residual((1+j)*N_dof1-1,:)=0;
    residual((1+j)*N_dof1,:)=0;
end
% ����г��ϵ����������
for i=1:2*N_harm*N_dof
    sensitivity_parameter_a1=zeros(2*N_harm*N_dof,1);
    x_a=zeros(N_dof,length(Tdata));dx_a=zeros(N_dof,length(Tdata));ddx_a=zeros(N_dof,length(Tdata));
    sensitivity_parameter_a1(i,1)=1;
    sensitivity_parameter_a1=reshape(sensitivity_parameter_a1,2,N_harm*N_dof);sensitivity_parameter_a1=sensitivity_parameter_a1';
    sensitivity_parameter_a=sensitivity_parameter_a1(1:N_harm,1:2);
    for num_dof=1:N_dof-1
        sensitivity_parameter_a=[sensitivity_parameter_a,sensitivity_parameter_a1(num_dof*N_harm+1:(num_dof+1)*N_harm,1:2)];%da(N_harm+1:2*N_harm,1:2);
    end
    for k0=1:N_dof
        for j=1:N_harm
            x_a(k0,:)=x_a(k0,:)+sensitivity_parameter_a(j,2*k0-1)*cos(vector_w(j,k0)*Tdata)+sensitivity_parameter_a(j,2*k0)*sin(vector_w(j,k0)*Tdata);
            dx_a(k0,:)=dx_a(k0,:)-vector_w(j,k0)*sensitivity_parameter_a(j,2*k0-1)*sin(vector_w(j,k0)*Tdata)+vector_w(j,k0)*sensitivity_parameter_a(j,2*k0)*cos(vector_w(j,k0)*Tdata);
            ddx_a(k0,:)=ddx_a(k0,:)-(vector_w(j,k0))^2*sensitivity_parameter_a(j,2*k0-1)*cos(vector_w(j,k0)*Tdata)-(vector_w(j,k0))^2*sensitivity_parameter_a(j,2*k0)*sin(vector_w(j,k0)*Tdata);
        end
    end
    residual(N_dof1*(i+N_w0)+1:N_dof1*(i+N_w0+1)-3,:)=M*ddx_a+K*x_a+3*k22*[(x(1,:)-x(2,:)).^2.*(x_a(1,:)-x_a(2,:));(x(2,:)-x(1,:)).^2.*(x_a(2,:)-x_a(1,:))];
    residual(N_dof1*(i+N_w0+1)-2,:)=m1*dx(1,:).*dx_a(1,:)+m2*dx(2,:).*dx_a(2,:)+k1*x(1,:).*x_a(1,:)+k2*x(2,:).*x_a(2,:)+k12*(x(1,:)-x(2,:)).*(x_a(1,:)-x_a(2,:))+k22*(x(1,:)-x(2,:)).^3.*(x_a(1,:)-x_a(2,:));
    if i<N_harm*N_dof
        residual(N_dof1*(i+N_w0+1)-1,:)=mod(i,2);
        residual(N_dof1*(i+N_w0+1),:)=0;
    else
        residual(N_dof1*(i+N_w0+1)-1,:)=0;
        residual(N_dof1*(i+N_w0+1),:)=mod(i,2);
    end
end

residual=residual';




