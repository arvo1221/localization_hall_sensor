function error=Berr(m_abs,x1,x2,x3,x4,x5,r_s,B_s,num)
global n
k=0;
error_local=zeros(n,n);
for i=1:n
    for j=1:n
        k=k+1;
        error_local(i,j)=norm(dipole_field(m_abs*[sin(x4)*cos(x5) sin(x4)*sin(x5) cos(x4)],[x1-r_s(k,1) x2-r_s(k,2) x3-r_s(k,3)],num)-B_s(k,:));
    end
end
error=sum(sum(error_local));
