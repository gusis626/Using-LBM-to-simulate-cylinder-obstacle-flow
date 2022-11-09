clear;
clc;
%setting global variables
global w;
w= [4/9,1/9,1/9,1/9,1/9,1/36,1/36,1/36,1/36];
global c;
c= [0,0;1,0;0,1;-1,0;0,-1;1,1;-1,1;-1,-1;1,-1];

%setting constant variables
NX= 500;NY= 100;
Lx= 500;Ly= 100;
obstX= Lx/5;obstY= Ly/2;obstR= 100/10;
[x,y]= meshgrid(1:NX+1,1:NY+1);
obst= (x-obstX).^2+(y-obstY).^2 <= obstR.^2;
obst([1,NY+1],:)= 1;
[row,col]= find(obst);
rho= zeros(NY+1,NX+1);
u= zeros(NY+1,NX+1);
v= zeros(NY+1,NX+1);
u0= zeros(NY+1,NX+1);
v0= zeros(NY+1,NX+1);
f= zeros(NY+1,NX+1,9);
F= zeros(NY+1,NX+1,9);
rho0= 1;
Re= 1500;
U= 0.02;
nu= U*2*obstR/Re;
tau= 3*nu+0.5;

%initialize
rho(:,:)= rho0;
u(:,1)= U;
for k=1:9
    f(:,:,k)= feq(k,rho,u,v);
end
tic;
n= 0;
objV= VideoWriter('example5');
objV.FrameRate= round(180);
open(objV)
H=600*3;
W=1920*3;

%main 
for a= 1:40000
    n= n+1;
    %collision plus streaming
    F(2:NY,2:NX,1) = f(2:NY,2:NX,1)+ (feq(1,rho(2:NY,2:NX),u(2:NY,2:NX),v(2:NY,2:NX))-f(2:NY,2:NX,1))/tau;
    F(2:NY,2:NX,3) = f(1:NY-1,2:NX,3)+ (feq(3,rho(1:NY-1,2:NX),u(1:NY-1,2:NX),v(1:NY-1,2:NX))-f(1:NY-1,2:NX,3))/tau;
    F(2:NY,2:NX,2) = f(2:NY,1:NX-1,2)+ (feq(2,rho(2:NY,1:NX-1),u(2:NY,1:NX-1),v(2:NY,1:NX-1))-f(2:NY,1:NX-1,2))/tau;
    F(2:NY,2:NX,5) = f(3:NY+1,2:NX,5)+ (feq(5,rho(3:NY+1,2:NX),u(3:NY+1,2:NX),v(3:NY+1,2:NX))-f(3:NY+1,2:NX,5))/tau;
    F(2:NY,2:NX,4) = f(2:NY,3:NX+1,4)+ (feq(4,rho(2:NY,3:NX+1),u(2:NY,3:NX+1),v(2:NY,3:NX+1))-f(2:NY,3:NX+1,4))/tau;
    F(2:NY,2:NX,6) = f(1:NY-1,1:NX-1,6)+ (feq(6,rho(1:NY-1,1:NX-1),u(1:NY-1,1:NX-1),v(1:NY-1,1:NX-1))-f(1:NY-1,1:NX-1,6))/tau;
    F(2:NY,2:NX,9) = f(3:NY+1,1:NX-1,9)+ (feq(9,rho(3:NY+1,1:NX-1),u(3:NY+1,1:NX-1),v(3:NY+1,1:NX-1))-f(3:NY+1,1:NX-1,9))/tau;
    F(2:NY,2:NX,8) = f(3:NY+1,3:NX+1,8)+ (feq(8,rho(3:NY+1,3:NX+1),u(3:NY+1,3:NX+1),v(3:NY+1,3:NX+1))-f(3:NY+1,3:NX+1,8))/tau;
    F(2:NY,2:NX,7) = f(1:NY-1,3:NX+1,7)+ (feq(7,rho(1:NY-1,3:NX+1),u(1:NY-1,3:NX+1),v(1:NY-1,3:NX+1))-f(1:NY-1,3:NX+1,7))/tau;
   for i=1:size(row)
        F(row(i),col(i),2)= f(row(i),col(i),4);
        F(row(i),col(i),3)= f(row(i),col(i),5);
        F(row(i),col(i),4)= f(row(i),col(i),2);
        F(row(i),col(i),5)= f(row(i),col(i),3);
        F(row(i),col(i),6)= f(row(i),col(i),8);
        F(row(i),col(i),7)= f(row(i),col(i),9);
        F(row(i),col(i),8)= f(row(i),col(i),6);
        F(row(i),col(i),9)= f(row(i),col(i),7);
    end

    %macro value
    u0= u;
    v0= v;
    u= zeros(NY+1,NX+1);
    v= zeros(NY+1,NX+1);
    rho(2:NY,2:NX)= 0;
    for k=1:9
         f(2:NY,2:NX,k) = F(2:NY,2:NX,k);
         u(2:NY,2:NX)=u(2:NY,2:NX)+f(2:NY,2:NX,k)*c(k,1);
         v(2:NY,2:NX)=v(2:NY,2:NX)+f(2:NY,2:NX,k)*c(k,2);
    end
    rho1= sum(f,3);
    rho(2:NY,2:NX)= rho1(2:NY,2:NX);
     u(2:NY,2:NX) = u(2:NY,2:NX)./rho(2:NY,2:NX);
     v(2:NY,2:NX) = v(2:NY,2:NX)./rho(2:NY,2:NX);

    %boundary
    %left/right bounce back
    rho(:,1)= rho(:,2);
    u(:,1)= U;
    rho(:,NX+1)= rho(:,NX);
    u(:,NX+1)= u(:,NX);
    v(:,NX+1)= 0;
    for k=1:9
        f(:,1,k) = feq(k,rho(:,1),u(:,1),v(:,1))+ f(:,2,k)- feq(k,rho(:,2),u(:,2),v(:,2));
        f(:,NX+1,k) = feq(k,rho(:,NX+1),u(:,NX+1),v(:,NX+1))+ f(:,NX,k)- feq(k,rho(:,NX),u(:,NX),v(:,NX));

    %convergence
  if mod(n,10)==0
     %drawing
    u_norm=sqrt(u(:,:).^2+v(:,:).^2); % 速度标量
startx=0:10:100;
starty=0:10:500;
[startx1,starty1]=meshgrid(startx,starty);
x1=1:101;
y1=1:501;
    frame = getframe(gcf);
    frame.cdata= imresize(frame.cdata,[H W]);
    subplot(2,1,1)
    pcolor(x,y,u_norm);  shading interp; 
    title(n);   
    axis equal off
        subplot(2,1,2)
    yy=plot(0,0)
delete(yy)
  yy= streamline(y1,x1,u(:,:),v(:,:),starty1,startx1)
axis equal off
  title('RE=1500')
axis equal off
    writeVideo(objV,frame);
    end
    end
end
toc;
close(objV)