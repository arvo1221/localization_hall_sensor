function state = sampling(x,y,z,theta,phi,num,method)
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
n=num;
[X,Y,Z]=meshgrid(linspace(-h*n/2,h*n/2,n),linspace(-h*n/2,h*n/2,n),1e-2);
%real = [x y z theta phi num method];
%% Actual position (r_c) and orientation (m_c) of the capsule
% magnetization direction
s = rng;
m_th=theta;
m_fi=phi;
m_c=m_abs*[sin(m_th)*cos(m_fi) sin(m_th)*sin(m_fi) cos(m_th)];

r_c=[x*1e-3 y*1e-3 z*1e-3];

%% measued B feild at sensors (B_s)
%Method 1 : interior-point + 3axis
%Method 2 : Levenberg-Marquardt + 3axis
%Method 3 : interior-point + 2axis
%Method 4 : Levenberg-Marquardt + 2axis
%Method 5 : interior-point + 1axis
%Method 6 : Levenberg-Marquardt + 1axis
mt = mod(7-method,3);
if mt == 0
    mt = 3;
end
k=0;
for i=1:n
    for j=1:n
        k=k+1;
        r_s(k,:)=[X(i,j) Y(i,j) Z(i,j)];
        B_s(k,:)=dipole_field(m_c,r_c-r_s(k,:),mt);
    end
end

%% find sol
fun = @(x) Berr(m_abs,x(1),x(2),x(3),x(4),x(5),r_s,B_s,mt);
lb=[-h*n/2 -h*n/2 1e-2 0 0];%lower bound
ub=[h*n/2 h*n/2 1e-1 pi 2*pi];%upper bound
A=[];b=[];
Aeq=[];beq=[];
nonlcon=[];
x0=[rand(1)*100e-3 rand(1)*100e-3 rand(1)*100e-3 rand(1)*pi rand(1)*2*pi];
if mod(method,2) == 1
    options = optimoptions(@fmincon,'Algorithm','interior-point','MaxIterations',1e+06,'OptimalityTolerance',1e-08,'StepTolerance',1e-08);
    [x,fval] = fmincon(fun,x0,A,b,Aeq,beq,lb,ub,nonlcon);
else
    options = optimoptions(@lsqnonlin,'Algorithm','Levenberg-Marquardt','MaxIterations',1e+06,'OptimalityTolerance',1e-08,'StepTolerance',1e-08);
    x = lsqnonlin(fun,x0,[],[],options);
end

%% result check
estimated=[x(1)*1000 x(2)*1000 x(3)*1000 x(4) x(5)];
state = estimated;

end
