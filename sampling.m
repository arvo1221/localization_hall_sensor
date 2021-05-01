function errs = sampling()
%% global parameter
global mu_0 mu_r n;

%% properties of a magnet inside the capsule
Br=1;
mu_r = 1;
mu_0 = 4*pi*10^-7; % magnetic permeability of air
V=(10e-3)^3;
m_abs=Br/mu_0*V;

%% magnetic sensor
h=5e-2;
n=8;
[X,Y,Z]=meshgrid(linspace(-h*n/2,h*n/2,n),linspace(-h*n/2,h*n/2,n),1e-2);

%% Actual position (r_c) and orientation (m_c) of the capsule
% magnetization direction
s = rng;
m_th=rand(1)*pi;
m_fi=rand(1)*2*pi;
m_c=m_abs*[sin(m_th)*cos(m_fi) sin(m_th)*sin(m_fi) cos(m_th)];

% position
r_c=[randi(10)*100e-3 randi(10)*100e-3 randi(10)*100e-3];

%% measued B feild at sensors (B_s)
k=0;
for i=1:n
    for j=1:n
        k=k+1;
        r_s(k,:)=[X(i,j) Y(i,j) Z(i,j)];
        B_s(k,:)=dipole_field(m_c,r_c-r_s(k,:),3);
        B_s2(k,:)=dipole_field(m_c,r_c-r_s(k,:),2);
        B_s1(k,:)=dipole_field(m_c,r_c-r_s(k,:),1);
    end
end

%% find sol
fun = @(x) Berr(m_abs,x(1),x(2),x(3),x(4),x(5),r_s,B_s,3);
lb=[0 0 0 0 0];%lower bound
ub=[+100e-2 +100e-2 +100e-2 pi 2*pi];%upper bound
A=[];b=[];
Aeq=[];beq=[];
nonlcon=[];

x0=[rand(1)*100e-3 rand(1)*100e-3 rand(1)*100e-3 rand(1)*pi rand(1)*2*pi];
%x0=[0 0 0 0 0];
options = optimoptions(@fmincon,'Algorithm','interior-point','MaxIterations',1e+06,'OptimalityTolerance',1e-08,'StepTolerance',1e-08);
[x,fval] = fmincon(fun,x0,A,b,Aeq,beq,lb,ub,nonlcon);

%% result check
x_estimated=[x(1)*1000 x(2)*1000 x(3)*1000 x(4) x(5)];
x_exact=[r_c*1000 m_th m_fi];
err(1,:) = x_exact - x_estimated;
err_norm(1) = norm(err(1,:));

%% Levenberg-Marquardt Method
fun = @(x) Berr(m_abs,x(1),x(2),x(3),x(4),x(5),r_s,B_s,3);
x0=[rand(1)*100e-3 rand(1)*100e-3 rand(1)*100e-3 rand(1)*pi rand(1)*2*pi];
%x0=[0 0 0 0 0];
options = optimoptions(@lsqnonlin,'Algorithm','Levenberg-Marquardt','MaxIterations',1e+06,'OptimalityTolerance',1e-08,'StepTolerance',1e-08);
x = lsqnonlin(fun,x0,[],[],options);

%% result check
x_estimated_LM=[x(1)*1000 x(2)*1000 x(3)*1000 x(4) x(5)];
err(2,:) = x_exact - x_estimated_LM;
err_norm(2) = norm(err(2,:));
err_norm(1) - err_norm(2);
%dipole_field(m_c,r_c-r_s(k,:),2);

%% find sol with 2axis
fun = @(x) Berr(m_abs,x(1),x(2),x(3),x(4),x(5),r_s,B_s2,2);

x0=[rand(1)*100e-3 rand(1)*100e-3 rand(1)*100e-3 rand(1)*pi rand(1)*2*pi];
options = optimoptions(@fmincon,'Algorithm','interior-point','MaxIterations',1e+06,'OptimalityTolerance',1e-08,'StepTolerance',1e-08);
[x,fval] = fmincon(fun,x0,A,b,Aeq,beq,lb,ub,nonlcon);

%% result check
x_estimated_2=[x(1)*1000 x(2)*1000 x(3)*1000 x(4) x(5)];
err(3,:) = x_exact - x_estimated_2;
err_norm(3) = norm(err(3,:));

%% Levenberg-Marquardt Method with 2axis
fun = @(x) Berr(m_abs,x(1),x(2),x(3),x(4),x(5),r_s,B_s2,2);
x0=[rand(1)*100e-3 rand(1)*100e-3 rand(1)*100e-3 rand(1)*pi rand(1)*2*pi];
options = optimoptions(@lsqnonlin,'Algorithm','Levenberg-Marquardt','MaxIterations',1e+06,'OptimalityTolerance',1e-08,'StepTolerance',1e-08);
x = lsqnonlin(fun,x0,[],[],options);

%% result check
x_estimated_LM2=[x(1)*1000 x(2)*1000 x(3)*1000 x(4) x(5)];
err(4,:) = x_exact - x_estimated_LM2;
err_norm(4) = norm(err(4,:));

%% find sol with 1axis
fun = @(x) Berr(m_abs,x(1),x(2),x(3),x(4),x(5),r_s,B_s1,1);

x0=[rand(1)*100e-3 rand(1)*100e-3 rand(1)*100e-3 rand(1)*pi rand(1)*2*pi];
options = optimoptions(@fmincon,'Algorithm','interior-point','MaxIterations',1e+06,'OptimalityTolerance',1e-08,'StepTolerance',1e-08);
[x,fval] = fmincon(fun,x0,A,b,Aeq,beq,lb,ub,nonlcon);

%% result check
x_estimated_1=[x(1)*1000 x(2)*1000 x(3)*1000 x(4) x(5)];
err(5,:) = x_exact - x_estimated_1;
err_norm(5) = norm(err(5,:));

%% Levenberg-Marquardt Method with 1axis
fun = @(x) Berr(m_abs,x(1),x(2),x(3),x(4),x(5),r_s,B_s1,1);
x0=[rand(1)*100e-3 rand(1)*100e-3 rand(1)*100e-3 rand(1)*pi rand(1)*2*pi];
options = optimoptions(@lsqnonlin,'Algorithm','Levenberg-Marquardt','MaxIterations',1e+06,'OptimalityTolerance',1e-08,'StepTolerance',1e-08);
x = lsqnonlin(fun,x0,[],[],options);

%% result check
x_estimated_LM1=[x(1)*1000 x(2)*1000 x(3)*1000 x(4) x(5)];
err(6,:) = x_exact - x_estimated_LM1;
err_norm(6) = norm(err(6,:));

errs = abs(err);

end

