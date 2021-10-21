clear all
close all
clc

warning('off')



load border_fit.mat
T = [];
method = 'spline';
b1 = @(mu) interp1(border_data(:,1)',border_data(:,2)',mu,method,'extrap');
a1 = @(mu) interp1(border_data(:,1)',border_data(:,4)',mu,method,'extrap');

b2 = @(mu) interp1(border_data(:,1)',border_data(:,6)',mu,method,'extrap');
a2 = @(mu) interp1(border_data(:,1)',border_data(:,8)',mu,method,'extrap');

b3 = @(mu) interp1(border_data(:,1)',border_data(:,10)',mu,method,'extrap');
a3 = @(mu) interp1(border_data(:,1)',border_data(:,12)',mu,method,'extrap');


%% panel A
nu_vec = 0:0.001:0.125;

N = 100;
A = 100;
B = 50;
R = rand(N,2);

for n = 1:N
    n/N
    mA = A*R(n,1);
    mB = B*R(n,2);
    sA = mA/3;
    sB = mB/3;
    r  = 2*rand()-1;

    mQ = 2*(1+mB)/mA;
    sQ = sqrt(sA^2+4*sB^2-4*r*sA*sB)/mA;
    p(n,1) = normcdf((1-mQ)/sQ);

    for i = 2:length(nu_vec)
        nu = nu_vec(i);
        m2 = (mB - b2(nu))/(a2(nu)*mA);
        s2 = sqrt((a2(nu)*sA)^2 - 2*r*a2(nu)*sA*sB + sB^2)/(a2(nu)*mA);

        m1 = (mB - b1(nu))/(a1(nu)*mA);
        s1 = sqrt((a1(nu)*sA)^2 - 2*r*a1(nu)*sA*sB + sB^2)/(a1(nu)*mA);

        p(n,i) = normcdf((1-m1)/s1) - normcdf((1-m2)/s2);
    end
end

figure()
plot(nu_vec,p)
ylim([0 1])






%% panel B
nu = 0.05;
N = 1e5;
hist_edges = [-10:0.01:10];

figure('Position',[0 0 1000 1000])

for scenario = 1:7
    tic
    switch scenario
        % grey
        case 1 
            mu_alpha = 60;
            mu_beta  = 18;
            rho      = 0.75; 
        
        % light red
        case 2
            mu_alpha = 60;
            mu_beta  = 21;
            rho      = 0.75;
            
        % dark red
        case 3 
            mu_alpha = 60;
            mu_beta  = 24;
            rho      = 0.75;
            
        % light green
        case 4 
            mu_alpha = 60;
            mu_beta  = 3.5;
            rho      = 0.75;
        
        % dark green    
        case 5 
            mu_alpha = 60;
            mu_beta  = 8;
            rho      = 0.75;
            
        % light purple
        case 6 
            mu_alpha = 150;
            mu_beta  = 20;
            rho      = 0.75;
        
        % dark purple    
        case 7 
            mu_alpha = 150;
            mu_beta  = 30;
            rho      = 0.75;
    end
    sigma_alpha = mu_alpha/5;
    sigma_beta  = mu_beta/5;
    
    mu = [mu_alpha ; mu_beta];
    S  = [sigma_alpha^2 rho*sigma_alpha*sigma_beta ; rho*sigma_alpha*sigma_beta sigma_beta^2];
    
    R = mvnrnd(mu,S,N);
    
    for i = 1:N
        alpha = R(i,1);
        beta  = R(i,2);

        x0 = [alpha 0];
        eq= EQ_converge(alpha,beta,nu,x0);
        D{scenario}(i,1) = eq.diff(end);

        x0 = [0 alpha];
        eq= EQ_converge(alpha,beta,nu,x0);
        D{scenario}(i,2) = eq.diff(end);

        x0 = [alpha alpha] + (rand(1,2)-0.5)*10;
        eq = EQ_converge(alpha,beta,nu,x0);
        D{scenario}(i,3) = eq.diff(end);
    end
    
    
    subplot(3,3,scenario)
    h = histogram(D{scenario}(:),hist_edges,'normalization','probability','EdgeColor','none');
    [scenario toc/60]
end






%% CD
nu = 0.05;
N = 1e4;
T_end = 200;
N_sample = 10;
hist_edges = [0:2:200];

ai = 0:1:100;
figure()
hold on
plot(ai,a1(nu)*ai + b1(nu))
plot(ai,a2(nu)*ai + b2(nu))
plot(ai,a3(nu)*ai + b3(nu))

grid on
axis square

xlim([0 100])
ylim([0 45])


for scenario = 1:40
    tic
    mu_alpha = 75;
    mu_beta  = scenario;
    rho      = 0.75;

    sigma_alpha = mu_alpha/5;
    sigma_beta  = mu_beta/5;
    
    mu = [mu_alpha ; mu_beta];
    S  = [sigma_alpha^2 rho*sigma_alpha*sigma_beta ; rho*sigma_alpha*sigma_beta sigma_beta^2];
    
    R = mvnrnd(mu,S,N);
    alpha = R(:,1);
    beta  = R(:,2);
    
    Q1 = (beta-b1(nu))./(a1(nu)*alpha);
    Q2 = (beta-b2(nu))./(a2(nu)*alpha);
    Q3 = (beta-b3(nu))./(a3(nu)*alpha);

    for i = 1:N
        if Q1(i) > 1
            T(i) = 1;
        elseif Q2(i) > 1
            T(i) = 2;
        elseif Q3(i) > 1
            T(i) = 3;
        else
            T(i) = 1;
        end
    end
    percentage(scenario,:) = [100*sum(T==1)/N 100*sum(T==2)/N 100*sum(T==3)/N];
    
    
    for i = 1:N
        d = Depth(alpha(i),beta(i),nu,T_end,N_sample);
        D_depth(scenario,i) = mean(d);       
    end
    [scenario toc/60]
    D_median(scenario) = median(D_depth(scenario,:));
end
figure()
plot(D_median./max(D_median),scenario)
set(gca,'XScale','log')
axis square
grid on

%% equilibrium converge
function EQ = EQ_converge(alpha,beta,nu,x0)
    thr  = 1e-2; 
    stop = 0;
    T_end = 10;
    EQ.yz = [0 0];
    while ~stop
        [t,x] = ode45(@(t,x) ODE_leaky(t,x,alpha,beta,nu), [0 T_end], x0);
        xf = x(end,:);
        if sqrt(sum((EQ.yz - xf).^2)) > thr & T_end < 200
            stop = 0;
            T_end = T_end*2;
            EQ.yz = NaN;
            EQ.x  = NaN;
            EQ.t  = NaN;
            EQ.diff = NaN;
        else
            stop = 1;
            EQ.yz = xf;
            EQ.x  = x;
            EQ.t  = t;
            EQ.diff = EQ.x(:,1)-EQ.x(:,2);
        end
    end
end



%% Depth
function d = Depth(alpha,beta,nu,T_end,N_sample)
    band = 0.01;
    
    [slope bias] = Separatrix(alpha,beta,nu,T_end,N_sample);

    if isnan(slope)
        d = [0 0];
    else
        yi = 0:max(alpha)/N_sample:max(alpha);
        zi = slope*yi + bias;
        zi1 = zi - band;
        zi2 = zi + band;

        positive_samples = find((yi>0) .* (zi1>0));

        y  = yi(positive_samples);
        z1 = zi1(positive_samples);
        z2 = zi2(positive_samples);

        for i = 1:length(y)
            [V yf zf ypz] = Potential(alpha,beta,nu,[y(i) z1(i)],T_end,1000);
            V1(i)  = V;
            Yf1(i) = yf;
            Zf1(i) = zf;

            [V yf zf ypz] = Potential(alpha,beta,nu,[y(i) z2(i)],T_end,1000);
            V2(i)  = V;
            Yf2(i) = yf;
            Zf2(i) = zf;
        end
        if isempty(y)
            d = [0 0];
        else
            d = [min(V1) min(V2)];
        end

    end
end


%% Separatrix
function [slope bias] = Separatrix(alpha,beta,nu,T_end,N_sample)

    answer = Bistable(alpha,beta,nu,T_end);
    if answer
        N = N_sample;
        y_vector = 0:alpha/(N-1):alpha;
        z_vector = 0:alpha/(N-1):alpha;
        [Y Z] = meshgrid(y_vector, z_vector);
        Y = Y+rand(N)*1e-4;
        Z = Z+rand(N)*1e-4;

        for i = 1:N
            for j = 1:N
                [t,x] = ode45(@(time,x) ODE_leaky(time,x,alpha,beta,nu),[0 T_end], [Y(i,j) Z(i,j)]);
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
function [V yf zf ypz] = Potential(alpha,beta,nu,x0,T_end,N)

    dt = T_end/N;

    ti = 0:dt:T_end;

    [t,x] = ode45(@(t,x) ODE_leaky(t,x,alpha,beta,nu), [0 T_end], x0);
    y = x(:,1);
    z = x(:,2);

    yf = y(length(y));
    zf = z(length(z));
    ypz = yf/zf;

    yi = interp1(t,y,ti);
    zi = interp1(t,z,ti);

    ey = nu + 1./(1+zi.^2);
    ez = nu + 1./(1+yi.^2);

    dydt = alpha*ey./(1 + beta*ey + beta*ez) - yi;
    dzdt = alpha*ez./(1 + beta*ey + beta*ez) - zi;

    V = sum(dydt.^2+dzdt.^2)*dt;

end

%% Bistable
function answer = Bistable(alpha,beta,nu,T_end)

    threshold = 0.001;
    [t,x] = ode45(@(time,x) ODE_leaky(time,x,alpha,beta,nu),[0 T_end], [alpha 0]);
    yf1 = x(length(t),1);
    zf1 = x(length(t),2);
    [t,x] = ode45(@(time,x) ODE_leaky(time,x,alpha,beta,nu),[0 T_end], [0 alpha]);
    yf2 = x(length(t),1);
    zf2 = x(length(t),2);
    
    answer = sqrt((yf1-yf2)^2+(zf1-zf2)^2)>threshold;
end


%% ODE of toggle with leakiness nu
function dx = ODE_leaky(t,x,alpha,beta,nu)

    dx = [0 ; 0];
    y = x(1);
    z = x(2);

    ey = nu+1/(1+z^2);
    ez = nu+1/(1+y^2);

    dx(1) = alpha*ey/(1 + beta*ey + beta*ez) - y;
    dx(2) = alpha*ez/(1 + beta*ey + beta*ez) - z;
end