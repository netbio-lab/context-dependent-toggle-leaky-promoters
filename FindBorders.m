% Generate the qi = 1 constraints and save the results as they take time 
% and are reused in multiple figures

clear all
close all
clc


a_vec  = [2:10:1000];
b_vec  = [0:0.1:100];
nu_vec = 0.01:0.005:0.125;
Nw     = 1e4;


for i = 1:length(nu_vec)
    i/length(nu_vec)
    nu = nu_vec(i);
    [mdl1 mdl2 mdl3 A{i} B{i} T{i}] = BorderFit(a_vec,b_vec,nu,Nw);
    M1{i} = mdl1;
    M2{i} = mdl2;
    M3{i} = mdl3;
    
end

for i = 1:length(nu_vec)
    if ~isempty(M1{i})
        m1  = [table2array(M1{i}.Coefficients(1,1:2)) table2array(M1{i}.Coefficients(2,1:2))];
        m1R = M1{i}.Rsquared.Ordinary;
    else
        m1  = [NaN NaN NaN NaN];
        m1R = NaN;
    end
    
    if ~isempty(M2{i})
        m2  = [table2array(M2{i}.Coefficients(1,1:2)) table2array(M2{i}.Coefficients(2,1:2))];
        m2R = M2{i}.Rsquared.Ordinary;
    else
        m2  = [NaN NaN NaN NaN];
        m2R = NaN;
    end
    
    if ~isempty(M3{i})
        m3  = [table2array(M3{i}.Coefficients(1,1:2)) table2array(M3{i}.Coefficients(2,1:2))];
        m3R = M3{i}.Rsquared.Ordinary;
    else
        m3  = [NaN NaN NaN NaN];
        m3R = NaN;
    end
    
    border_data(i,:) = [nu_vec(i) m1 m2 m3];
    border_data_R(i,:) = [nu_vec(i) m1R m2R m3R];
end

clearvars -except border_data
save('border_fit.mat')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Determine stability profile in Fig 3

function [mdl1 mdl2 mdl3 A B T] = BorderFit(a_vec,b_vec,nu,Nw)


    [A,B] = meshgrid(a_vec, b_vec);
    [row, column] = size(A);


    for r = 1:row
        for c = 1:column
            a = A(r,c);
            b = B(r,c);
            w = logspace(-2,2,Nw);
            x_nu = Curve_w(nu,w);
            x_ab = Curve_abw(a,b,nu,w);

            w3 = Find_w3(a/(1+2*b*nu),b/(1+2*b*nu),[1:a/Nw:a]);
            w1 = sqrt(((1-2*nu) - sqrt(1-8*nu))/(2*nu));
            w2 = sqrt(((1-2*nu) + sqrt(1-8*nu))/(2*nu));


            wi = logspace(-2,log10(w3),Nw);
            nu_above_ab = sum(interp1(w,max(x_nu),wi) > interp1(w,max(x_ab),wi));


            roots = [w3 w1 w2];
            [value index] = sort(roots);


            type = [];
            switch find(index==1)
                case 1
                    type = 1;
                case 2
                    type = 2;
                case 3
                    if nu_above_ab
                        type = 3;
                    else
                        type = 1;
                    end
            end
            T(r,c) = type;

        end
    end
    T2 = T;
    T3 = T;



    % Bistable region
    T2(T==1) = 0;
    T2(T==2) = 1;
    T2(T==3) = 0;
    [row,column] = size(T2);
    [~,ind] = max(T2,[],2);
    first = ind;
    [~,ind] = max(fliplr(T2),[],2);
    second = column-ind+1;
    first(first == 1) = 1;
    second(second == column) = 1;


    x1 = a_vec(first);
    x1(x1==a_vec(1)) = NaN;
    y1 = b_vec;
    del_ind = find(isnan(x1));
    x1(del_ind) = [];
    y1(del_ind) = [];

    x2 = a_vec(second);
    x2(x2==a_vec(1)) = NaN;
    y2 = b_vec;
    del_ind = find(isnan(x2));
    x2(del_ind) = [];
    y2(del_ind) = [];

    if ~isempty(x1) & ~isempty(y1)
        mdl1 = fitlm(x1,y1);
    else
        mdl1 = [];
    end
    if ~isempty(x2) & ~isempty(y2)
        mdl2 = fitlm(x2,y2);
    else
        mdl2 = [];
    end



    % Tristable region
    T3(T==1) = 0;
    T3(T==2) = 0;
    T3(T==3) = 1;
    [row,column] = size(T3);

    [~,ind] = max(fliplr(T3),[],2);
    second = column-ind+1;
    second(second == column) = 1;


    x2 = a_vec(second);
    x2(x2==a_vec(1)) = NaN;
    y2 = b_vec;
    del_ind = find(isnan(x2));
    x2(del_ind) = [];
    y2(del_ind) = [];

    if ~isempty(x2) & ~isempty(y2)
        mdl3 = fitlm(x2,y2);
    else
        mdl3 = [];
    end





    figure()
    view([0 90])
    axis square
    hold on
    surf(A,B,T,'edgecolor','none')
    if ~isempty(mdl1)
        plot(a_vec,predict(mdl1,a_vec'),'m','LineWidth',5)
    end
    if ~isempty(mdl2)
        plot(a_vec,predict(mdl2,a_vec'),'m','LineWidth',5)
    end
    if ~isempty(mdl3)
        plot(a_vec,predict(mdl3,a_vec'),'c:','LineWidth',5)
    end

end



%% Find w3
function w3 = Find_w3(a,b,w)
    D = (1+w.^2).^2 + 2*b*(1+w.^2) - 2*a*w;
    [value, index] = min(abs(D));
    w3 = w(index);
end

%% Constraint 1
function x = Curve_w(w,z)
    yp = (z + (- 4*w^2*z.^4 - 8*w^2*z.^2 - 4*w^2 - 4*w*z.^2 - 4*w + z.^2).^(1/2))./(2*(w*z.^2 + w));
    ym = (z - (- 4*w^2*z.^4 - 8*w^2*z.^2 - 4*w^2 - 4*w*z.^2 - 4*w + z.^2).^(1/2))./(2*(w*z.^2 + w));

    index = find(- 4*w^2*z.^4 - 8*w^2*z.^2 - 4*w^2 - 4*w*z.^2 - 4*w + z.^2 < 0);

    yp(index) = NaN;
    ym(index) = NaN;

    x = [yp; ym];
end

%% Constraint 2
function x = Curve_abw(a,b,w,z)
    a = a/(1+2*b*w);
    b = b/(1+2*b*w);

    yp = (a + sqrt(a^2-4*(1+z.^2+b).*(1+z.^2+2*b+b*z.^2-a*z)))./(2*(1+z.^2+b));
    ym = (a - sqrt(a^2-4*(1+z.^2+b).*(1+z.^2+2*b+b*z.^2-a*z)))./(2*(1+z.^2+b));

    index = find(a^2-4*(1+z.^2+b).*(1+z.^2+2*b+b*z.^2-a*z) < 0 );
    yp(index) = NaN;
    ym(index) = NaN;

    x = [yp ; ym];
end


