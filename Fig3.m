clear all
close all
clc



nu    = 0.03;
T_end = 400;


figure()
hold on



y_vec = 0:0.5:10;
z_vec = y_vec;
[Y0,Z0] = meshgrid(y_vec,z_vec);
[row,column] = size(Y0);
Y0 = Y0 + rand(row,column)/1e2;
Z0 = Z0 + rand(row,column)/1e2;
ti = 0:T_end/1e3:T_end;





alpha  = 200;
beta   = 9.8;
yi = [];
zi = [];
for r = 1:row
    for c = 1:column
        x0 = [Y0(r,c) Z0(r,c)];
        [t,x] = ode45(@(t,x) ODE_leaky(t,x,alpha,beta,nu), [0 T_end], x0);
        yi = [yi ; interp1(t,x(:,1),ti)];
        zi = [zi ; interp1(t,x(:,2),ti)];
    end
end
yz_diff = yi-zi;
plot(ti,yz_diff,'b') 


alpha  = 200;
beta   = 15;
yi = [];
zi = [];
for r = 1:row
    for c = 1:column
        x0 = [Y0(r,c) Z0(r,c)];
        [t,x] = ode45(@(t,x) ODE_leaky(t,x,alpha,beta,nu), [0 T_end], x0);
        yi = [yi ; interp1(t,x(:,1),ti)];
        zi = [zi ; interp1(t,x(:,2),ti)];
    end
end
yz_diff = yi-zi;
plot(ti,yz_diff,'g')  

alpha  = 200;
beta   = 5;
yi = [];
zi = [];
for r = 1:row
    for c = 1:column
        x0 = [Y0(r,c) Z0(r,c)];
        [t,x] = ode45(@(t,x) ODE_leaky(t,x,alpha,beta,nu), [0 T_end], x0);
        yi = [yi ; interp1(t,x(:,1),ti)];
        zi = [zi ; interp1(t,x(:,2),ti)];
    end
end
yz_diff = yi-zi;
plot(ti,yz_diff,'r') 



axis square
grid on





%% ODE
function dx = ODE_leaky(t,x,alpha,beta,nu)

    dx = [0 ; 0];
    y = x(1);
    z = x(2);

    ey = nu+1/(1+z^2);
    ez = nu+1/(1+y^2);

    dx(1) = alpha*ey/(1 + beta*ey + beta*ez) - y;
    dx(2) = alpha*ez/(1 + beta*ey + beta*ez) - z;
end