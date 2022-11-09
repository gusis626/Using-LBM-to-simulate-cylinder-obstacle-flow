function [re]= feq(k,rho,u,v)
    global w;
    global c;
    cu= u*c(k,1)+v*c(k,2);
    uv= u.*u+v.*v;
    re= w(k)*rho.*(1+3*cu+4.5*cu.*cu-1.5*uv);
end