function B = dipole_field(m,r,num)
%% Global Parameter
global mu_r mu_0;
%% Main
B = mu_r*mu_0./(4*pi*norm(r)^3).*(3*dot(m,r)*r/norm(r)^2-m);
if num==2
    B(:,3) = [];
end
if num==1
   B(:,3) = [];
   B(:,2) = [];
end
