clear all
close all
clc

mu      = @(mu_alpha,mu_beta) [mu_alpha ; mu_beta];
S       = @(sigma_alpha, sigma_beta, rho) [sigma_alpha^2 rho*sigma_alpha*sigma_beta; rho*sigma_alpha*sigma_beta sigma_beta^2];
mu_Q    = @(mu_alpha,mu_beta) 2*(1+mu_beta)./mu_alpha;
sigma_Q = @(mu_alpha,sigma_alpha, sigma_beta, rho, q) sqrt(q.^2*sigma_alpha^2 + 4*sigma_beta^2 - 4*q*rho*sigma_alpha*sigma_beta)./mu_alpha;
F_Q     = @(mu_Q, sigma_Q, q) normcdf((q-mu_Q)./sigma_Q);

edges = linspace(0,1.5,1e3);
q = edges;
N = 1e5;


figure()

subplot(1,2,1)
hold on

mu_alpha = 10;
mu_beta  = 0;
sigma_alpha = mu_alpha/5;
sigma_beta  = mu_beta/5;
rho = 0;

R = mvnrnd(mu(mu_alpha,mu_beta),S(sigma_alpha, sigma_beta, rho),N);
alpha = R(:,1);
beta  = R(:,2);
Q = 2*(1+beta)./alpha;

cdf = F_Q(mu_Q(mu_alpha,mu_beta), sigma_Q(mu_alpha,sigma_alpha, sigma_beta, rho, q), q);
pdf = diff(cdf);

h1 = histogram(Q,edges,'EdgeColor','none','normalization','probability')
plot(q(1:end-1),pdf,'r')
pdf1 = pdf;



mu_alpha = 10;
mu_beta  = 2;
sigma_alpha = mu_alpha/5;
sigma_beta  = mu_beta/5;
rho = 0;

R = mvnrnd(mu(mu_alpha,mu_beta),S(sigma_alpha, sigma_beta, rho),N);
alpha = R(:,1);
beta  = R(:,2);
Q = 2*(1+beta)./alpha;

cdf = F_Q(mu_Q(mu_alpha,mu_beta), sigma_Q(mu_alpha,sigma_alpha, sigma_beta, rho, q), q);
pdf = diff(cdf);

h2 = histogram(Q,edges,'EdgeColor','none','normalization','probability')
plot(q(1:end-1),pdf,'r')
pdf2 = pdf;

mu_alpha = 10;
mu_beta  = 4;
sigma_alpha = mu_alpha/5;
sigma_beta  = mu_beta/5;
rho = 0;

R = mvnrnd(mu(mu_alpha,mu_beta),S(sigma_alpha, sigma_beta, rho),N);
alpha = R(:,1);
beta  = R(:,2);
Q = 2*(1+beta)./alpha;

cdf = F_Q(mu_Q(mu_alpha,mu_beta), sigma_Q(mu_alpha,sigma_alpha, sigma_beta, rho, q), q);
pdf = diff(cdf);

h3 = histogram(Q,edges,'EdgeColor','none','normalization','probability')
plot(q(1:end-1),pdf,'r')
pdf3 = pdf;

axis square 
grid on


subplot(1,2,2)
hold on

mu_alpha = 10;
mu_beta  = 2;
sigma_alpha = mu_alpha/5;
sigma_beta  = mu_beta/5;




rho = -1;
R = mvnrnd(mu(mu_alpha,mu_beta),S(sigma_alpha, sigma_beta, rho),N);
alpha = R(:,1);
beta  = R(:,2);
Q = 2*(1+beta)./alpha;

cdf = F_Q(mu_Q(mu_alpha,mu_beta), sigma_Q(mu_alpha,sigma_alpha, sigma_beta, rho, q), q);
pdf = diff(cdf);

h4 = histogram(Q,edges,'EdgeColor','none','normalization','probability')
plot(q(1:end-1),pdf,'r')
pdf4 = pdf;




rho = -1/3;
R = mvnrnd(mu(mu_alpha,mu_beta),S(sigma_alpha, sigma_beta, rho),N);
alpha = R(:,1);
beta  = R(:,2);
Q = 2*(1+beta)./alpha;

cdf = F_Q(mu_Q(mu_alpha,mu_beta), sigma_Q(mu_alpha,sigma_alpha, sigma_beta, rho, q), q);
pdf = diff(cdf);

h5 = histogram(Q,edges,'EdgeColor','none','normalization','probability')
plot(q(1:end-1),pdf,'r')
pdf5 = pdf;




rho = 1/3;
R = mvnrnd(mu(mu_alpha,mu_beta),S(sigma_alpha, sigma_beta, rho),N);
alpha = R(:,1);
beta  = R(:,2);
Q = 2*(1+beta)./alpha;

cdf = F_Q(mu_Q(mu_alpha,mu_beta), sigma_Q(mu_alpha,sigma_alpha, sigma_beta, rho, q), q);
pdf = diff(cdf);

h6 = histogram(Q,edges,'EdgeColor','none','normalization','probability')
plot(q(1:end-1),pdf,'r')
pdf6 = pdf;




rho = +1;
R = mvnrnd(mu(mu_alpha,mu_beta),S(sigma_alpha, sigma_beta, rho),N);
alpha = R(:,1);
beta  = R(:,2);
Q = 2*(1+beta)./alpha;

cdf = F_Q(mu_Q(mu_alpha,mu_beta), sigma_Q(mu_alpha,sigma_alpha, sigma_beta, rho, q), q);
pdf = diff(cdf);

h7 = histogram(Q,edges,'EdgeColor','none','normalization','probability')
plot(q(1:end-1),pdf,'r')
pdf7 = pdf; 

axis square 
grid on
