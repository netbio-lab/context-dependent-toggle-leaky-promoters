clear all
close all
clc


figure()


subplot(1,2,1)
alpha = linspace(0,150,100);
beta = alpha/2-1;
plot(alpha,beta)
axis square
grid on
xlim([0 150])
ylim([0 60])



subplot(1,2,2)
w = linspace(1,10,100);
alpha_w = (1+w.^2).^2./(2*w);
nu_w    = (w.^2-1)./(1+w.^2).^2;
plot(alpha_w,nu_w)
axis square
grid on
xlim([0 150])
ylim([0 0.15])

