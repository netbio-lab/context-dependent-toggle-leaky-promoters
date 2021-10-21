clear all
close all
clc

FontSize = 20;

a = 20;
b = 0;

N = 100;
T_end = 100;
N_time = 1000;


x_vector = 0:1.5*a/(N-1):1.5*a;
y_vector = x_vector;
[X Y] = meshgrid(x_vector, y_vector);
X = X+rand(N)*1e-4;
Y = Y+rand(N)*1e-4;

for i = 1:N
    100*i/N
    for j = 1:N
        v = Potential(a,b,[X(i,j) Y(i,j)],T_end,N_time);
        V(i,j) = v;
    end
end
[saddle_value saddle_index] = min(diag(V));
V = 1 + V - min(V(:));




ax = figure()
hAxes = axes;
%ax(1) = subplot(2,2,1)
surf(X,Y,log10(V),'FaceAlpha',0.75,'EdgeColor','none')
view([-45 20])
% colorbar
xlabel('y')
ylabel('z')
axis square
colormap(ax,othercolor('Greys9'))


hold on

N_time = 1e4;
dt = T_end/N_time;
ti = 0:dt:T_end;

[t,x] = ode45(@(t,x) ODE_toggle(t,x,a,b), [0 T_end], [30 20]);
y = x(:,1);
z = x(:,2);
yi = interp1(t,y,ti);
zi = interp1(t,z,ti);
ey = 1./(1+zi.^2);
ez = 1./(1+yi.^2);
dydt = a*ey./(1 + b*ey + b*ez) - yi;
dzdt = a*ez./(1 + b*ey + b*ez) - zi;
V = sum(dydt.^2+dzdt.^2)*dt;
dV = (dydt.^2+dzdt.^2)*dt;
Vi(1) = V;
for i = 2:length(dV)
    Vi(i) = Vi(i-1) - dV(i-1);
end
Vi = 1 + Vi - min(Vi);

plot3(yi,zi,log10(Vi)+0.05,'r','Color',[255 44 121]/255,'LineWidth',10)

[t,x] = ode45(@(t,x) ODE_toggle(t,x,a,b), [0 T_end], [20 30]);
y = x(:,1);
z = x(:,2);
yi = interp1(t,y,ti);
zi = interp1(t,z,ti);
ey = 1./(1+zi.^2);
ez = 1./(1+yi.^2);
dydt = a*ey./(1 + b*ey + b*ez) - yi;
dzdt = a*ez./(1 + b*ey + b*ez) - zi;
V = sum(dydt.^2+dzdt.^2)*dt;
dV = (dydt.^2+dzdt.^2)*dt;
Vi(1) = V;
for i = 2:length(dV)
    Vi(i) = Vi(i-1) - dV(i-1);
end
Vi = 1 + Vi - min(Vi);

plot3(yi,zi,log10(Vi)+0.05,'r','Color',[32 178 170]/255,'LineWidth',10)


set(gca,'ZLim',[0 3],'FontSize',FontSize,'XTickLabel',[],'YTickLabel',[],'XTick',[0 10 20 30],'YTick',[0 10 20 30])
Ax = gca;
Ax.ZAxis.Visible = 'off';
Ax.XGrid = 'off';
Ax.YGrid = 'off';
Ax.ZGrid = 'off';




% Plot Fig 4B

h1 = 0.545;
h2 = 2.039;
h = @(q,h1,h2) h1*(1./q-1).^h2;
q = 0:0.01:1;

figure()
semilogy(q,h(q,h1,h2),'b')
ylim([1e-4 1e4])



% Plot Fig 4C
T_end = 200;
N_sample = 100;
beta_vec = 0:0.1:2;


alpha_vec = [10 7.5 5];

for i = 1:length(alpha_vec)
    alpha = alpha_vec(i)
    for j = 1:length(beta_vec)
        j/length(beta_vec)
        beta = beta_vec(j);

        d = Depth(alpha,beta,T_end,N_sample);
        D1(i,j) = mean(d);
        D2(i,j) = h(2*(1+beta)/alpha,h1,h2)'; 
    end
end

for i = 1:length(D2(:))
    if ~isreal(D2(i))
        D2(i) = 0;
    end
end
D1(isnan(D1)) = 0;

figure()
plot(beta_vec,D1./D1(2,1))
hold on
plot(beta_vec,D2./D2(2,1),'ko')
ylim([0 2.5])










%% Depth
function d = Depth(a,b,T_end,N_sample)
    band = 0.01;
    
    [slope bias] = Separatrix(a,b,T_end,N_sample);

    if isnan(slope)
        d = [NaN NaN];
    else
        yi = 0:max(a)/N_sample:max(a);
        zi = slope*yi + bias;
        zi1 = zi - band;
        zi2 = zi + band;

        positive_samples = find((yi>0) .* (zi1>0));

        y  = yi(positive_samples);
        z1 = zi1(positive_samples);
        z2 = zi2(positive_samples);

        for i = 1:length(y)
            [V yf zf ypz] = Potential(a,b,[y(i) z1(i)],T_end,1000);
            V1(i)  = V;
            Yf1(i) = yf;
            Zf1(i) = zf;

            [V yf zf ypz] = Potential(a,b,[y(i) z2(i)],T_end,1000);
            V2(i)  = V;
            Yf2(i) = yf;
            Zf2(i) = zf;
        end
        d = [min(V1) min(V2)];

    end
end


%% Separatrix
function [slope bias] = Separatrix(a,b,T_end,N_sample)

    answer = Bistable(a,b,T_end);
    if answer
        N = N_sample;
        y_vector = 0:a/(N-1):a;
        z_vector = 0:a/(N-1):a;
        [Y Z] = meshgrid(y_vector, z_vector);
        Y = Y+rand(N)*1e-4;
        Z = Z+rand(N)*1e-4;

        for i = 1:N
            for j = 1:N
                [t,x] = ode45(@(time,x) ODE_toggle(time,x,a,b),[0 T_end], [Y(i,j) Z(i,j)]);
                yf = x(length(t),1);
                zf = x(length(t),2);
                Yf(i,j)  = yf;
                Zf(i,j)  = zf;
            end
        end
        predictor = [Y(:) Z(:)];
        label = kmeans([Yf(:) Zf(:)],2);


        Mdl = fitclinear(predictor,label);
        slope = -Mdl.Beta(1)/Mdl.Beta(2);
        bias = -Mdl.Bias/Mdl.Beta(2);

    else
        slope = NaN;
        bias  = NaN;
    end
end


%% Potential depth
function [V yf zf ypz] = Potential(a,b,x0,T_end,N)

    dt = T_end/N;

    ti = 0:dt:T_end;

    [t,x] = ode45(@(t,x) ODE_toggle(t,x,a,b), [0 T_end], x0);
    y = x(:,1);
    z = x(:,2);

    yf = y(length(y));
    zf = z(length(z));
    ypz = yf/zf;

    yi = interp1(t,y,ti);
    zi = interp1(t,z,ti);

    ey = 1./(1+zi.^2);
    ez = 1./(1+yi.^2);

    dydt = a*ey./(1 + b*ey + b*ez) - yi;
    dzdt = a*ez./(1 + b*ey + b*ez) - zi;

    V = sum(dydt.^2+dzdt.^2)*dt;

end

%% Bistable
function answer = Bistable(a,b,T_end)

    threshold = 0.001;
    [t,x] = ode45(@(time,x) ODE_toggle(time,x,a,b),[0 T_end], [a 0]);
    yf1 = x(length(t),1);
    zf1 = x(length(t),2);
    [t,x] = ode45(@(time,x) ODE_toggle(time,x,a,b),[0 T_end], [0 a]);
    yf2 = x(length(t),1);
    zf2 = x(length(t),2);
    
    answer = sqrt((yf1-yf2)^2+(zf1-zf2)^2)>threshold;
end

%% ODE
function dx = ODE_toggle(t,x,a,b)

    dx = [0 ; 0];
    y = x(1);
    z = x(2);

    ey = 1/(1+z^2);
    ez = 1/(1+y^2);

    dx(1) = a*ey/(1 + b*ey + b*ez) - y;
    dx(2) = a*ez/(1 + b*ey + b*ez) - z;
end