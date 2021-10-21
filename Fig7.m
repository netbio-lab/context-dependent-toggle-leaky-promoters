clear all
close all
clc



%% A
nu = 0.03;
N = 1e4;
hist_edges = [-10:0.01:10];

figure('Position',[0 0 1000 1000])

for scenario = 1:6
    tic
    switch scenario
        
        % light red
        case 1
            mu_alpha = 10;
            mu_beta  = 1;
            mu_c     = 0;
            
        % dark red
        case 2
            mu_alpha = 10;
            mu_beta  = 1;
            mu_c     = 3;
            
        % light green
        case 3 
            mu_alpha = 20;
            mu_beta  = 2;
            mu_c     = 0;
        
        % dark green    
        case 4 
            mu_alpha = 20;
            mu_beta  = 2;
            mu_c     = 3;
            
        % light purple
        case 5 
            mu_alpha = 150;
            mu_beta  = 5;
            mu_c     = 0;
        
        % dark purple    
        case 6 
            mu_alpha = 150;
            mu_beta  = 5;
            mu_c     = 15;
    end
    rho      = 0.75;
    rho1     = 0.75;
    rho2     = 0.75;
    sigma_alpha = mu_alpha/5;
    sigma_beta  = mu_beta/5;
    sigma_c     = mu_c/5;
    
    mu = [mu_alpha ; mu_beta ; mu_c];
    S  = [sigma_alpha^2 rho*sigma_alpha*sigma_beta rho1*sigma_alpha*sigma_c; rho*sigma_alpha*sigma_beta sigma_beta^2 rho2*sigma_beta*sigma_c ; rho1*sigma_alpha*sigma_c rho2*sigma_beta*sigma_c sigma_c^2];
    
    R = mvnrnd(mu,S,N);
    
    for i = 1:N
        alpha  = R(i,1);
        beta   = R(i,2);
        beta_c = R(i,3);

        x0 = [alpha/(1+beta_c) 0];
        eq = EQ_converge(alpha/(1+beta_c),beta/(1+beta_c),nu,x0);
        D{scenario}(i,1) = eq.diff(end);

        x0 = [0 alpha/(1+beta_c)];
        eq = EQ_converge(alpha/(1+beta_c),beta/(1+beta_c),nu,x0);
        D{scenario}(i,2) = eq.diff(end);

        x0 = [alpha/(1+beta_c) alpha/(1+beta_c)] + (rand(1,2)-0.5)*10;
        eq = EQ_converge(alpha/(1+beta_c),beta/(1+beta_c),nu,x0);
        D{scenario}(i,3) = eq.diff(end);
    end
    
    
    subplot(2,3,scenario)
    h = histogram(D{scenario}(:),hist_edges,'normalization','probability','EdgeColor','none');
    [scenario toc/60]
end








%% B
T_end = 200;
N_sample = 100;

nu = 0;
bc_vec = 0:0.1:2;

h1 = 0.545;
h2 = 2.039;
h = @(q) h1*(1./q-1).^h2;

b(1) = 2;
a(1) = 10;
q = 2*(1+b)/a;

b(2) = 1;
a(2) = 2*(1+b(2))/q;

b(3) = 0.2;
a(3) = 2*(1+b(3))/q;


for i = 1:length(bc_vec)
    i/length(bc_vec)
    
    Bc = bc_vec(i);
    for j = 1:length(a)
        
        A = a(j)/(1+Bc);
        B = b(j)/(1+Bc);
        d = Depth(A,B,nu,T_end,N_sample);
        D1(i,j) = mean(d);
        
        D2(i,j) = h(2*(1+B)/A)';
        
    end
    
end

for i = 1:length(D2(:))
    if ~isreal(D2(i))
        D2(i) = 0;
    end
end


figure()
plot(bc_vec,D1./D1(1,:))
hold on
plot(bc_vec,D2./D2(1,:),'ko')
set(gca,'YScale','log')
ylim([1e-2 1])




%% Depth
function d = Depth(alpha,beta,nu,T_end,N_sample)
    band = 0.01;
    [slope bias] = Separatrix(alpha,beta,nu,T_end,N_sample);

    if isnan(slope)
        d = [NaN NaN];
    else
        yi = 0:alpha/N_sample:alpha;
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
        d = [min(V1) min(V2)];

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

