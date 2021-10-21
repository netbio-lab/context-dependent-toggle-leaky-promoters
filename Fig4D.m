clear all
close all
clc


load border_fit.mat


b1 = @(nu) interp1(prism(:,1)',prism(:,2)',nu);
a1 = @(nu) interp1(prism(:,1)',prism(:,4)',nu);

b2 = @(nu) interp1(prism(:,1)',prism(:,6)',nu);
a2 = @(nu) interp1(prism(:,1)',prism(:,8)',nu);

b3 = @(nu) interp1(prism(:,1)',prism(:,10)',nu);
a3 = @(nu) interp1(prism(:,1)',prism(:,12)',nu);



nu_index = 3;
load h12.mat
nu = nu_vec(nu_index)



T_end = 1e5;
dt    = 1e-1;
zci = @(v) find(diff(sign(v)));


% Drift function
f = @(t,V,a,b,nu,sigma)[a*(nu+1/(1+V(2)^2))/(1 + b(1)*(nu+1/(1+V(2)^2)+nu+1/(1+V(1)^2))) - V(1); ...
                        a*(nu+1/(1+V(1)^2))/(1 + b(1)*(nu+1/(1+V(2)^2)+nu+1/(1+V(1)^2))) - V(2)];

% Diffusion function
g = @(t,sigma)[sigma; sigma];      

% Time vector
t = 0:dt:T_end;



figure()
subplot(1,2,1)
hold on

subplot(1,2,2)
hold on


for k = 1:3
    switch k
        case 1
            alpha = 25;
            beta_vec_fine{k} = 0:0.1:8;
            beta_vec{k} = 0:1:8;
            color_actual{k} = 'r';
            
        case 2
            alpha = 50;
            beta_vec_fine{k} = 0:0.1:20;
            beta_vec{k} = 0:1:20;
            color_actual{k} = 'g';
        
        case 3
            alpha = 75;
            beta_vec_fine{k} = 2:0.1:31;
            beta_vec{k} = 2:1:30;
            color_actual{k} = 'b';

    end
    
    beta_vec_fine{k} = 0:0.01:30;
    beta_vec{k} = 0:0.5:30;
    subplot(1,2,1)
    h_depth{k} = interp2(A,B,H1{nu_index},alpha*ones(1,length(beta_vec_fine{k})),beta_vec_fine{k});
    plot(beta_vec_fine{k},h_depth{k},color_actual{k})
    
    
    for s = 1:3
        switch s
            case 1
                sigma = 1;
            case 2
                sigma = 1.1;
            case 3
                sigma = 1.2;
        end
                
        ZC = [];
        for N_rep = 1:3
            for j = 1:length(beta_vec{k})

                beta = beta_vec{k}(j);
                b = beta;
                a = alpha;

                % initial condition
                A0 = [a*randn() ; 0];                
                % Create and initialize state vector
                lt = length(t);
                n = length(A0);
                L = zeros(n,lt);
                L(:,1) = A0;

                % Set seed and pre-calculate Wiener increments
                rng(1);
                r = sqrt(dt)*randn(lt-1,n).';

                % General Euler-Maruyama integration loop
                for i = 1:lt-1
                    L(:,i+1) = L(:,i)+f(t(i),L(:,i),a,b,nu,sigma)*dt+r(:,i).*g(t(i),sigma);
                end

                y = L(1,:);
                z = L(2,:);

                ZC(N_rep,j) = length(zci(y-z));
            end
        end
        tau{k}(s,:) = T_end./mean(ZC,1);
        subplot(1,2,2)
        plot(beta_vec{k},tau{k}(s,:),color_actual{k})
    end
    
end


for k = 1:3
    h_depth{k}(isnan(h_depth{k})) = 0;
    h_max(k) = max(h_depth{k});
    tau_max(k) = max(tau{k}(:));
end

prism1 = [beta_vec_fine{1}' 100*h_depth{1}'/max(h_max)];
prism2 = [beta_vec_fine{2}' 100*h_depth{2}'/max(h_max)];
prism3 = [beta_vec_fine{3}' 100*h_depth{3}'/max(h_max)];

prism4 = [beta_vec{1}' 100*tau{1}'/max(tau_max)];
prism5 = [beta_vec{2}' 100*tau{2}'/max(tau_max)];
prism6 = [beta_vec{3}' 100*tau{3}'/max(tau_max)];


