%Pricing Heston Utilizing Various Discretization Scheme%

s0 = 100;
K = 110;
T = 1;
eta = 0.1;
r = 0.025;
q = 0.0125;
mu0 = 0.05;
theta = 0.0625;
alpha = 1.5;
kappa = 2.75;
lambda = 0.0125;
rho = -0.65;


Euler_payOff1 = 0;
Milsten_payOff2 = 0;
RK_payOff3 = 0;
payOff4 = 0;

n = 10;
N = 2^n;


nSims = 100000;
m  = 52;
dt = T/m;
%%
for j = 1:nSims
    
    s1 = s0;
    s2 = s0;
    s3 = s0;
    mu1 = mu0;
    mu2 = mu0;
    mu3 = mu0;
    
    for i = 1:m
        z1 = randn;
        z2 = randn;
        % Euler
        %-------------------------------------------------------
        mu1 = mu1+ kappa*(theta-mu1)*dt + lambda*sqrt(mu1*dt)*z2;
        s1 = s1 + (r-q)*s1*dt + s1*sqrt(mu1*dt)*z1;    
        %--------------------------------------------------------
        % Milstein
        %---------------------------------------------------------
        mu2 = mu2 + kappa*(theta-mu2)*dt + lambda*sqrt(mu2*dt)*z2 + 0.25*dt*(z2^2-1);
        s2 = s2 + (r-q)*s2*dt + s2*sqrt(mu2*dt)*z1 + 0.5*mu2*s2*dt*(z1^2-1);
        %----------------------------------------------------------
        % Runge-Kutta
        %-----------------------------------------------------------
        mu3_tilde =  mu3 + kappa*(theta-mu3)*dt + lambda*sqrt(mu3*dt);
        mu3 = mu3 + kappa*(theta-mu3)*dt + lambda*sqrt(mu3*dt)*z2 + 0.5*(lambda*sqrt(mu3)*(mu3_tilde - mu3))*(z2^2-1)*sqrt(dt);
        s3_tilde = s3 + (r-q)*s3*dt + s3*sqrt(mu3*dt);
        s3 = s3 + (r-q)*s3*dt + s3*sqrt(mu3*dt)*z1 + 0.5*(sqrt(mu3)*s3_tilde - sqrt(mu3)*s3)*(z1^2-1)*sqrt(dt);
        %---------------------------------------------------------------
    end
      
    Euler_payOff1 = Euler_payOff1 + max(s1-K,0);
    Milsten_payOff2 = Milsten_payOff2 + max(s2-K,0);
    RK_payOff3 = RK_payOff3 + max(s3-K,0);
    
end

c1 = real(exp(-r*T)*Euler_payOff1/nSims);
c2 = real(exp(-r*T)*Milsten_payOff2/nSims);
c3 = real(exp(-r*T)*RK_payOff3/nSims);
%%
params = [kappa theta lambda rho mu0];
call = genericFFT('Heston', eta, alpha, N, s0, K, r, q, T, params);
disp("The actual price using Heston Model is: "+call);

disp("The price using Euler Method to simulate the Heston Model is: "+c1);
disp("The price using Milstein Method to simulate the Heston Model is: "+c2);
disp("The price using  Runge-Kutta Method to simulate the Heston Model is: "+c3);

% Output:
% The actual price using Heston Model is: 6.1881
% The price using Euler Method to simulate the Heston Model is: 6.24
% The price using Milstein Method to simulate the Heston Model is: 6.1848
% The price using  Runge-Kutta Method to simulate the Heston Model is: 6.2403