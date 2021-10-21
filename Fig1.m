clear all
close all
clc

a_tog = 4;
b_tog = 0.1;

a_rep = 10;
b_rep = 1;
nu = 0;


T_end = 1e4;
dt    = 1e-1; 
sigma = 0.1;


figure('Position',[100 100 1200 800])



subplot(2,2,1)
hold on

a = [a_tog a_rep];
b = [b_tog 0];
[X,Y,n_hist,C0,n_hist_red,C0_red,n_hist_green,C0_green,n_hist_purple1,C0_purple1,n_hist_purple2,C0_purple2] = Subpops(a,b,nu,T_end,dt,sigma);
surf(X,Y,n_hist_red,C0_red,'EdgeColor','none')
surf(X,Y,n_hist_green,C0_green,'EdgeColor','none')

a = [a_tog a_rep];
b = [b_tog b_rep];
[X,Y,n_hist,C0,n_hist_red,C0_red,n_hist_green,C0_green,n_hist_purple1,C0_purple1,n_hist_purple2,C0_purple2] = Subpops(a,b,nu,T_end,dt,sigma);
surf(X,Y,n_hist,C0,'EdgeColor','none')

shading interp
view(2)
axis square

a_tog = 8;
b_tog = 1.2;

subplot(2,2,2)
hold on

a = [a_tog a_rep];
b = [b_tog 0];
[X,Y,n_hist,C0,n_hist_red,C0_red,n_hist_green,C0_green,n_hist_purple1,C0_purple1,n_hist_purple2,C0_purple2] = Subpops(a,b,nu,T_end,dt,sigma);
surf(X,Y,n_hist_red,C0_red,'EdgeColor','none')
surf(X,Y,n_hist_green,C0_green,'EdgeColor','none')

a = [a_tog a_rep];
b = [b_tog b_rep];
[X,Y,n_hist,C0,n_hist_red,C0_red,n_hist_green,C0_green,n_hist_purple1,C0_purple1,n_hist_purple2,C0_purple2] = Subpops(a,b,nu,T_end,dt,sigma);
surf(X,Y,n_hist,C0,'EdgeColor','none')


shading interp
view(2)
axis square



% Drift function
f = @(t,V,a,b,nu,sigma)[a(1)*(nu+1/(1+V(2)^2))/(1 + b(1)*(nu+1/(1+V(2)^2)+nu+1/(1+V(1)^2)) + b(2)*(1/(1+V(5)^4)+1/(1+V(3)^4)+1/(1+V(4)^4))) - V(1); ...
                        a(1)*(nu+1/(1+V(1)^2))/(1 + b(1)*(nu+1/(1+V(2)^2)+nu+1/(1+V(1)^2)) + b(2)*(1/(1+V(5)^4)+1/(1+V(3)^4)+1/(1+V(4)^4))) - V(2); ...
                        a(2)*1/(1+V(5)^4)/(1 + b(1)*(nu+1/(1+V(2)^2)+nu+1/(1+V(1)^2)) + b(2)*(1/(1+V(5)^4)+1/(1+V(3)^4)+1/(1+V(4)^4)))      - V(3); ...
                        a(2)*1/(1+V(3)^4)/(1 + b(1)*(nu+1/(1+V(2)^2)+nu+1/(1+V(1)^2)) + b(2)*(1/(1+V(5)^4)+1/(1+V(3)^4)+1/(1+V(4)^4)))      - V(4); ...
                        a(2)*1/(1+V(4)^4)/(1 + b(1)*(nu+1/(1+V(2)^2)+nu+1/(1+V(1)^2)) + b(2)*(1/(1+V(5)^4)+1/(1+V(3)^4)+1/(1+V(4)^4)))      - V(5)];

% Diffusion function
g = @(t)[sigma; sigma ; 0 ; 0 ; 0];      





% Time vector
t = 0:dt:T_end; 

subplot(2,2,3)
b_rep = 0;
a = [a_tog a_rep];
b = [b_tog b_rep];

% initial condition
A0 = [1 ; 1 ; rand(3,1)];                
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
    L(:,i+1) = L(:,i)+f(t(i),L(:,i),a,b,nu,sigma)*dt+r(:,i).*g(t(i));
end

y = L(1,:);
z = L(2,:);

hold on
plot(t,y,'ro')
plot(t,z,'go')
axis square
grid on
ylim([0 4])


prism = [t' y' z'];


subplot(2,2,4)
b_rep = 1;
a = [a_tog a_rep];
b = [b_tog b_rep];

% initial condition
A0 = [1 ; 1 ; rand(3,1)];                
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
    L(:,i+1) = L(:,i)+f(t(i),L(:,i),a,b,nu,sigma)*dt+r(:,i).*g(t(i));
end

y = L(1,:);
z = L(2,:);

hold on
plot(t,y,'ro')
plot(t,z,'go')
axis square
grid on
ylim([0 4])

prism = [prism y' z'];


%% subpopulations
function [X,Y,n_hist,C0,n_hist_red,C0_red,n_hist_green,C0_green,n_hist_purple1,C0_purple1,n_hist_purple2,C0_purple2] = Subpops(a,b,nu,T_end,dt,sigma)

    % Drift function
    f = @(t,V,a,b,nu,sigma)[a(1)*(nu+1/(1+V(2)^2))/(1 + b(1)*(nu+1/(1+V(2)^2)+nu+1/(1+V(1)^2)) + b(2)*(1/(1+V(5)^4)+1/(1+V(3)^4)+1/(1+V(4)^4))) - V(1); ...
                            a(1)*(nu+1/(1+V(1)^2))/(1 + b(1)*(nu+1/(1+V(2)^2)+nu+1/(1+V(1)^2)) + b(2)*(1/(1+V(5)^4)+1/(1+V(3)^4)+1/(1+V(4)^4))) - V(2); ...
                            a(2)*1/(1+V(5)^4)/(1 + b(1)*(nu+1/(1+V(2)^2)+nu+1/(1+V(1)^2)) + b(2)*(1/(1+V(5)^4)+1/(1+V(3)^4)+1/(1+V(4)^4)))      - V(3); ...
                            a(2)*1/(1+V(3)^4)/(1 + b(1)*(nu+1/(1+V(2)^2)+nu+1/(1+V(1)^2)) + b(2)*(1/(1+V(5)^4)+1/(1+V(3)^4)+1/(1+V(4)^4)))      - V(4); ...
                            a(2)*1/(1+V(4)^4)/(1 + b(1)*(nu+1/(1+V(2)^2)+nu+1/(1+V(1)^2)) + b(2)*(1/(1+V(5)^4)+1/(1+V(3)^4)+1/(1+V(4)^4)))      - V(5)];
                            
    % Diffusion function
    g = @(t)[sigma; sigma ; 0 ; 0 ; 0];      
       
    % Time vector
    t = 0:dt:T_end; 
    
    
    
    
    % initial condition
    A0 = [a(1) ; 0 ; rand(3,1)];                
    % Create and initialize state vector (L here is transposed relative to sde_euler output)
    lt = length(t);
    n = length(A0);
    L = zeros(n,lt);
    L(:,1) = A0;

    % Set seed and pre-calculate Wiener increments with order matching sde_euler
    rng(1);
    r = sqrt(dt)*randn(lt-1,n).';

    % General Euler-Maruyama integration loop
    for i = 1:lt-1
        L(:,i+1) = L(:,i)+f(t(i),L(:,i),a,b,nu,sigma)*dt+r(:,i).*g(t(i));
    end

    y = L(1,:);
    z = L(2,:);




    % initial condition
    A0 = [0 ; a(1); rand(3,1)];                
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
        L(:,i+1) = L(:,i)+f(t(i),L(:,i),a,b,nu,sigma)*dt+r(:,i).*g(t(i));
    end

    y = [y L(1,:)];
    z = [z L(2,:)];

    load colorData

    Nq = 100;
    color1_q = othercolor('PuRd10',Nq+1);
    color2_q = othercolor('Greens10',Nq+1);
    color3_q = othercolor('Purples10',Nq+1);


    edges = 0:0.025:4;
    [N_hist,Xbins,Ybins] = hist2d(y,z,edges,edges);
    [X,Y] = meshgrid(Xbins,Ybins);
    n_hist = 1+round(Nq*(N_hist - min(N_hist(:)))/max(N_hist(:)));
    [row,column] = size(n_hist);


    n_hist_green = n_hist - triu(n_hist,1) + triu(ones(row,column),1);
    n_hist_red   = n_hist - tril(n_hist,-1) + tril(ones(row,column),-1);
    n_hist_purple1   = n_hist - tril(n_hist,-1) + tril(ones(row,column),-1);
    n_hist_purple2   = n_hist - triu(n_hist,1) + triu(ones(row,column),1);


    C0(:,:,1) = reshape(color3_q(n_hist(:),1),row,column);
    C0(:,:,2) = reshape(color3_q(n_hist(:),2),row,column);
    C0(:,:,3) = reshape(color3_q(n_hist(:),3),row,column);

    C0_purple1(:,:,1) = reshape(color3_q(n_hist_purple1(:),1),row,column);
    C0_purple1(:,:,2) = reshape(color3_q(n_hist_purple1(:),2),row,column);
    C0_purple1(:,:,3) = reshape(color3_q(n_hist_purple1(:),3),row,column);

    C0_purple2(:,:,1) = reshape(color3_q(n_hist_purple2(:),1),row,column);
    C0_purple2(:,:,2) = reshape(color3_q(n_hist_purple2(:),2),row,column);
    C0_purple2(:,:,3) = reshape(color3_q(n_hist_purple2(:),3),row,column);

    C0_green(:,:,1) = reshape(color2_q(n_hist_green(:),1),row,column);
    C0_green(:,:,2) = reshape(color2_q(n_hist_green(:),2),row,column);
    C0_green(:,:,3) = reshape(color2_q(n_hist_green(:),3),row,column);

    C0_red(:,:,1) = reshape(color1_q(n_hist_red(:),1),row,column);
    C0_red(:,:,2) = reshape(color1_q(n_hist_red(:),2),row,column);
    C0_red(:,:,3) = reshape(color1_q(n_hist_red(:),3),row,column);

end







% hist2d
function varargout=hist2d(x,y,varargin)
    %HIST2D Bivariable histogram plot
    %	HIST2D(X,Y) creates a bivariate histogram plot of vectors X and Y.
    %
    %	The function uses an automatic binning algorithm that returns bins with
    %	a uniform area, chosen to cover the range of elements in X and Y and 
    %	reveal the underlying shape of the distribution. HIST2D without output
    %	argument displays the bins as 3-D rectangular bars such that the height
    %	of each bar indicates the number of elements in the bin.
    %
    %	HIST2D(X,Y,NBINS) specifies the number of bins to use in each dimension
    %	of the histogram (default is 10).
    %
    %	HIST2D(X,Y,Xedges,Yedges) specifies the edges of the bins in each 
    %	dimension using the vectors Xedges and Yedges.
    %
    %	HIST2D(...,'tile') plots the result as a tiled 2-D image.
    %
    %	HIST2D(...,'bar3') plots the result as a 3-D bars. This is the default
    %	graph output but without this option, it will automatically switch to
    %	2-D if the number of total bins exceeds 2500 (e.g. 50x50).
    %
    %	N = HIST2D(...) returns the bin counts matrix N of size MxP where M is
    %	number of bins for Y and P the number of bins for X. No graph produced.
    %	Add 'bar3' or 'tile' option to force the graph.
    %
    %	[N,Xbins,Ybins] = HIST2D(...) returns also the two bins vectors.
    %
    %	It is also possible to normalize the bin counts matrix:
    %
    %	HIST2D(...,'probability') normalizes bin counts as a probability. The 
    %	height of each bar is the relative number of observations (number of 
    %	observations in bin / total number of observations). The sum of the bar
    %	heights is 1.
    %
    %	HIST2D(...,'countdensity') normalizes bin counts as count density. The 
    %	height of each bar is (number of observations in bin) / (area of bin). 
    %	The volume (height * area) of each bar is the number of observations 
    %	in the bin. The sum of the bar volumes is equal to numel(X) and numel(Y).
    %
    %	HIST2D(...,'pdf') normalizes bin counts as probability density function.
    %	The height of each bar is (number of observations in the bin) / (total 
    %	number of observations * area of bin). The volume of each bar is the 
    %	relative number of observations. The sum of the bar volumes is 1.
    %
    %	HIST2D(...,'cumcount') normalizes bin counts as cumulative counts. The 
    %	height of each bar is the cumulative number of observations in each bin
    %	and all previous bins in both the X and Y dimensions. The height of the
    %	last bar is equal to numel(X) and numel(Y).
    %
    %	HIST2D(...,'cdf') normalizes bin counts as cumulative density function. 
    %	The height of each bar is equal to the cumulative relative number of 
    %	observations in each bin and all previous bins in both the X and Y 
    %	dimensions. The height of the last bar is 1.
    %
    %
    %	Example:
    %		x = randn(1000,1);
    %		y = randn(1000,1);
    %		hist2d(x,y)
    %
    %
    %	Author: Francois Beauducel <beauducel@ipgp.fr>
    %	Created: 2018-03-24 in Yogyakarta, Indonesia
    %	Updated: 2018-03-27

    %	Copyright (c) 2018, François Beauducel, covered by BSD License.
    %	All rights reserved.
    %
    %	Redistribution and use in source and binary forms, with or without 
    %	modification, are permitted provided that the following conditions are 
    %	met:
    %
    %	   * Redistributions of source code must retain the above copyright 
    %	     notice, this list of conditions and the following disclaimer.
    %	   * Redistributions in binary form must reproduce the above copyright 
    %	     notice, this list of conditions and the following disclaimer in 
    %	     the documentation and/or other materials provided with the distribution
    %	                           
    %	THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" 
    %	AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE 
    %	IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE 
    %	ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER OR CONTRIBUTORS BE 
    %	LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR 
    %	CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF 
    %	SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS 
    %	INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN 
    %	CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) 
    %	ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE 
    %	POSSIBILITY OF SUCH DAMAGE.

    if nargin < 2
        error('Not enough input variable.')
    end

    if ~isvector(x) || ~isvector(y) || numel(x) ~= numel(y)
        error('X and Y must be vectors of the same length.')
    end

    if nargin > 2 && isscalar(varargin{1}) && round(varargin{1}) > 0
        nbins = round(varargin{1});
    else
        nbins = 10;
    end

    if nargin > 3 && isnumeric(varargin{1}) && isvector(varargin{1}) ...
                  && isnumeric(varargin{2}) && isvector(varargin{2})
        xedges = varargin{1};
        yedges = varargin{2};
    else
        xedges = linspace(min(x),max(x),nbins+1);
        yedges = linspace(min(y),max(y),nbins+1);
    end

    % plot options flag
    plot_bar3 = any(strcmpi(varargin,'bar3'));
    plot_tile = any(strcmpi(varargin,'tile'));

    % computes bins vectors (as middle of each edges couples)
    xbins = mean(cat(1,xedges(1:end-1),xedges(2:end)));
    ybins = mean(cat(1,yedges(1:end-1),yedges(2:end)));

    % computes bins width vectors and area matrix
    xbw = diff(xedges);
    ybw = diff(yedges);
    [xx,yy] = meshgrid(xbw,ybw);
    a = xx.*yy;

    % initiate the result matrix
    n = zeros(length(ybins),length(xbins));

    % main loop to fill the matrix with element counts
    for i = 1:size(n,1)
        k = find(y >= yedges(i) & y < yedges(i+1));
        for j = 1:size(n,2)
            n(i,j) = length(find(x(k) >= xedges(j) & x(k) < xedges(j+1)));
        end
    end

    % normalize options
    if any(strcmpi(varargin,'countdensity'))
        n = n./a;
    elseif any(strcmpi(varargin,'cumcount'))
        n = cumsum(cumsum(n,1),2);
    elseif any(strcmpi(varargin,'probability'))
        n = n/sum(n(:));
    elseif any(strcmpi(varargin,'pdf'))
        n = n./a/sum(n(:));
    elseif any(strcmpi(varargin,'cdf'))
        n = cumsum(cumsum(n,1),2)/sum(n(:));	
    end

    % plots a 3-D graph with indexed colors
    if nargout < 1 || plot_bar3 || plot_tile
        if plot_tile || (numel(n) > 2500 && ~plot_bar3)
            imagesc(xbins,ybins,n)
            hold on
            plot(x,y,'.k','MarkerSize',10)
            hold off
        else
            % unit cube XYZ coordinates for patch
            ux = [0 1 1 0 0 0;1 1 0 0 1 1;1 1 0 0 1 1;0 1 1 0 0 0];
            uy = [0 0 1 1 0 0;0 1 1 0 0 0;0 1 1 0 1 1;0 0 1 1 1 1];
            uz = [0 0 0 0 0 1;0 0 0 0 0 1;1 1 1 1 0 1;1 1 1 1 0 1];

            if ~ishold
                cla
            end
            for i = 1:size(n,1)
                for j = 1:size(n,2)
                    patch(ux*xbw(j) + xedges(j), ...
                          uy*ybw(i) + yedges(i), ...
                          uz*n(i,j),repmat(n(i,j)/max(n(:)),size(ux)))
                end
            end
            axis tight
            view(3)
            rotate3d on
        end
    end

    if nargout > 0
        varargout{1} = n;
    end

    if nargout > 2
        varargout{2} = xbins;
        varargout{3} = ybins;
    end

end