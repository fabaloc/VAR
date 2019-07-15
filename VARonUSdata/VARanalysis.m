%% VAR Analysis of US output growth, inflation and short-term interest rate

clear all;
close all;

%% 1) Load data 
load('VAR3.mat')

%% 2) Estimation of a VAR model with lags p={1,2...,7,8} with constant

% Firstly, I run VAR on all 8 lags: in the case of trivariate of lag 8, we would need 8 presample values
% so Z0 is for 1 to 8, in sample starts at position 9
L = length(gdp);
K = 3;
% Example with a precise number of lag, done with p=1:8 afterwards
p = 8;
T = L-p; %In-sample
Y = zeros(K,T);
Z = ones((K*p+1),T); %first row is already 1 to estimate the constant vector

for i = 1:T
    % Filling Y
    Y(1,i) = gdp(p+i);
    Y(2,i) = inf(p+i);
    Y(3,i) = str(p+i);
    % Filling Z
    for j = 0:(p-1)
    Z(2+K*j,i) = gdp(p+i-j-1);
    Z(3+K*j,i) = inf(p+i-j-1);
    Z(4+K*j,i) = str(p+i-j-1);
    end
end
Ahat = Y*Z'/(Z*Z');
Uhat = Y - Ahat*Z;

% c hat (the constant) is the first column of Ahat, A1hat would be the 3*3 matrix
% Ahat(2:4,1:3)

%%
% Now if I want to do it for every lag, it was coded so I can run the loop
% on the previous code:
pmax = 8;
Aglobal = zeros(K*pmax,(K*pmax+1));  % Global matrix where values are stored

for p = 1:8
    T = L-p; %In-sample
    Y = zeros(K,T);
    Z = ones((K*p+1),T); %first row is already 1 for estimating constant
    
    for i = 1:T
        % Filling Y
        Y(1,i) = gdp(p+i);
        Y(2,i) = inf(p+i);
        Y(3,i) = str(p+i);
        % Filling Z
        for j = 0:(p-1)
            Z(2+K*j,i) = gdp(p+i-j-1);
            Z(3+K*j,i) = inf(p+i-j-1);
            Z(4+K*j,i) = str(p+i-j-1);
        end
    end
    Ahat = Y*Z'/(Z*Z');
    Aglobal(1+K*(p-1):(K*p),1:length(Ahat))= Ahat;
end;

%% 3) Computing AIC and BIC

% Will use Aglobal to compute matrix Sigma hat for BIC
pmax = 8;
AIC = zeros(1,pmax);
BIC = zeros(1,pmax);
for p = 1:8
    T = L-p; %In-sample
    Y = zeros(K,T);
    Z = ones((K*p+1),T); %first row is already 1 for estimating constant
    
    for i = 1:T
        % Filling Y
        Y(1,i) = gdp(p+i);
        Y(2,i) = inf(p+i);
        Y(3,i) = str(p+i);
        % Filling Z
        for j = 0:(p-1)
            Z(2+K*j,i) = gdp(p+i-j-1);
            Z(3+K*j,i) = inf(p+i-j-1);
            Z(4+K*j,i) = str(p+i-j-1);
        end
    end
    Ahat = Y*Z'/(Z*Z');
    %Get fitted error term to have variance covariance matrix E
    Uhat = Y - Ahat*Z;
    % E = (Uhat*Uhat')/(length(Uhat)-1) % -1 to correct for bias but we can
    % get it with the function of Matlab
    E = cov(Uhat');               % it indeed gives th same result
    AIC(p) = log(det(E)) + 2/T*(K^2*p+K);
    BIC(p) = log(det(E)) + log(T)/T*(K^2*p+K);
end;

subplot(2,1,1)
plot(AIC)
title('AIC')
subplot(2,1,2)
plot(BIC)
title('BIC')
%% 4) Estimation of Var(2)

L = length(gdp);
K = 3;
p = 2;
T = L-p;
Y = zeros(K,T);
Z = ones((K*p+1),T); 

for i = 1:T
    Y(1,i) = gdp(p+i);
    Y(2,i) = inf(p+i);
    Y(3,i) = str(p+i);
    for j = 0:(p-1)
    Z(2+K*j,i) = gdp(p+i-j-1);
    Z(3+K*j,i) = inf(p+i-j-1);
    Z(4+K*j,i) = str(p+i-j-1);
    end
end
Ahat = Y*Z'/(Z*Z');

% I use the companion form to find stationarity. The companion matrix is
% given by:
% [A1, A2;
%   I,  0]
% in the case of a trivariate

companion = zeros(6,6);
A1 = Ahat(1:3,2:4);
companion(1:3,1:3) = A1;      
A2 = Ahat(1:3,5:7);
companion(1:3,4:6) = A2;      
companion(4:6,1:3) = eye(3);  % Identity matrix I(3)

% find the eigenvalues 
e = eig(companion);

% check eigenvalues setting a counter counting how many of them are above 1
counter = 0;
for i = 1:length(e)
    % disp(abs(e(i))); %check that they are all under 1 in absolute value
    if abs(e(i))> 1
        counter = counter + 1;
    end
end
disp(counter); % counter is 0, no eigenvalue is above 1

% I found that the 6 eigenvalues are below 1 => The VAR(2) with a constant
% is stationary !

%% 5)Plotting actual and fitted data

%gather all the data
for i = 1:length(gdp)
    Y(1,i) = gdp(i);
    Y(2,i) = inf(i);
    Y(3,i) = str(i);
end

% p=2 so we need two initial values and the first fitted value is for the 3
% value of Y
p = 2;
T = length(gdp)-p;
yhat = zeros(K,T);
cons = Ahat(1:3,1);
for i = 1:T
    % There's a two lag difference between yhat and Y so ith column if Y
    % is i'th second lag in yhat, similary i+1 in Y is first lag for yhat
    yhat(1:3,i) = cons + A1*Y(1:3,i+1) + A2*Y(1:3,i);
end


% Final plot
subplot(3,1,1)
plot(yhat(1,1:end))
hold on
plot(Y(1,3:end)) %starts at 3 because of 2 lags
hold off
title('GDP')

subplot(3,1,2)
plot(yhat(2,1:end))
hold on
plot(Y(2,3:end)) %starts at 3 because of 2 lags
hold off
title('Inflation')

subplot(3,1,3)
plot(yhat(3,1:end))
hold on
plot(Y(3,3:end)) %starts at 3 because of 2 lags
hold off
title('Short term rate')

%% 6) Autocorrelation of residuals up to order 20

eps = Y(1:3,3:end) - yhat;
subplot(3,1,1)
autocorr(eps(1,1:end),20);
title('GDP')
subplot(3,1,2)
autocorr(eps(2,1:end),20);
title('Inflation')
subplot(3,1,3)
autocorr(eps(3,1:end),20);
title('Short term rate')

% Three error terms seem to be uncorrelated in time which is what we look for