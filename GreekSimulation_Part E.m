


% ============================================================
%    Calculating Greeks via Simulation (Black-Merton-Scholes)
%
%           (1) w/o Common Random Numbers
%           (2) w/ Common Random Numbers
%           (3) via Pathwise Estimator
%           (4) via Likelihood Ratio
%          
% -------------------------------------------------------------


% Set up parameters
spot = 100;
strike = 115;
interest = 0.04;
dividend = 0.015;
volatility = 0.28;
time = 1;
barrier = 130;
monitoringTimes = 12;

delS = 0.1;
delSig = 0.001;

Simulations = 500000;

fprintf("Problem 2 \n")
fprintf("Digital Call \n")

% Simulation

z1 = randn(Simulations,1);
z2 = randn(Simulations,1);
z3 = randn(Simulations,1);

% Delta and Gamma

tmp1U = (spot+delS)*exp((interest-dividend-volatility^2/2)*time);
tmp1  =  spot      *exp((interest-dividend-volatility^2/2)*time);
tmp1D = (spot-delS)*exp((interest-dividend-volatility^2/2)*time);

tmp2 = volatility*sqrt(time);

% (1) W/ Common Random Numbers
s = tmp1*exp(tmp2*z1);
tmp = exp(-interest*time)*(s>strike);

s = tmp1U*exp(tmp2*z1);
tmpU = exp(-interest*time)*(s>strike);

s = tmp1D*exp(tmp2*z1);
tmpD = exp(-interest*time)*(s>strike);
%
delC_CRN1 = mean((tmpU-tmpD))/(2*delS);
%
gam_CRN  = mean((tmpU-2*tmp+tmpD))/(delS^2);

% (2) W/O Common Random Numbers
s =  tmp1*exp(tmp2*z1);
tmp = exp(-interest*time)*(s>strike);

s = tmp1U*exp(tmp2*z2);
tmpU = exp(-interest*time)*(s>strike);

s = tmp1D*exp(tmp2*z3);
tmpD = exp(-interest*time)*(s>strike);

delC_tilde1 = mean((tmpU-tmpD))/(2*delS);

gam_tilde  = mean((tmpU-2*tmp+tmpD))/(delS^2);

% (3) Pathwise Estimator

% Not continuous
delC_pathwise = NaN;
gam_pathwise = NaN;

% ----------------------------------
% (4) Likelihood ratio
s =  tmp1*exp(tmp2*z1);
scoreH = exp(-interest*time)*(s>strike).*z1/(spot*volatility*sqrt(time));
delC_likelihood = mean(scoreH);

scoreG = exp(-interest*time)*(s>strike).*((z1.^2-1)/(spot^2*volatility^2*time)-z1/(spot^2*volatility*sqrt(time)));
gam_likelihood = mean(scoreG);



%             Vega


tmp1U  =  spot *exp((interest-dividend-(volatility+delSig)^2/2)*time);
tmp1D  =  spot *exp((interest-dividend-(volatility-delSig)^2/2)*time);

tmp2U = (volatility+delSig)*sqrt(time);
tmp2D = (volatility-delSig)*sqrt(time);

% (1) W/ Common Random Numbers
s = tmp1U*exp(tmp2U*z1);
tmpU = exp(-interest*time)*(s>strike);

s = tmp1D*exp(tmp2D*z1);
tmpD = exp(-interest*time)*(s>strike);

vega_CRN = mean((tmpU-tmpD)/(2*delSig));

% (2) W/O Common Random Numbers
s = tmp1U*exp(tmp2U*z1);
tmpU = exp(-interest*time)*(s>strike);

s = tmp1D*exp(tmp2D*z2);
tmpD = exp(-interest*time)*(s>strike);

vega_tilde = mean((tmpU-tmpD)/(2*delSig));

% (3) Pathwise Estimator
% Not continuous 
vega_pathwise = NaN;

% (4) Likelihood Ratio

s =  tmp1*exp(tmp2*z1);
scoreH = exp(-interest*time)*(s>strike).*((z1.^2-1)/volatility - z1*sqrt(time));
vega_likelihood = mean(scoreH);



% Displaying Results

fprintf("The value of Delta is \n")
disp([ ' W/ CRN: ' num2str(delC_CRN1),  ' W/O CRN: ', num2str(delC_tilde1), ' Pathwise Estimator: ', num2str(delC_pathwise), ' Likelihood Ratio: ', num2str(delC_likelihood)]);

fprintf("The value of Gamma is \n")
disp([' W/ CRN: ' num2str(gam_CRN),  ' W/O CRN: ', num2str(gam_tilde), ' Pathwise Estimator: ', num2str(gam_pathwise), ' Likelihood Ratio: ', num2str(gam_likelihood)]);

fprintf("The value of Vega is \n")
disp([' W/ CRN: ' num2str(vega_CRN),  ' W/O CRN: ', num2str(vega_tilde), ' Pathwise Estimator: ', num2str(vega_pathwise), ' Likelihood Ratio: ', num2str(vega_likelihood)]);




fprintf("Up-and-Out Call \n")

% Simulation

%         Delta & Gamma

tmp1U = (spot+delS)*exp((interest-dividend-volatility^2/2)*time);
tmp1  =  spot      *exp((interest-dividend-volatility^2/2)*time);
tmp1D = (spot-delS)*exp((interest-dividend-volatility^2/2)*time);

tmp2 = volatility*sqrt(time);


% (1) W/ Common Random Numbers
s = tmp1*exp(tmp2*z1);
tmp = exp(-interest*time)*max(s-strike,0).*(s<barrier);

s = tmp1U*exp(tmp2*z1);
tmpU = exp(-interest*time)*max(s-strike,0).*(s<barrier);

s = tmp1D*exp(tmp2*z1);
tmpD = exp(-interest*time)*max(s-strike,0).*(s<barrier);

delC_CRN1 = mean((tmpU-tmpD))/(2*delS);

gam_CRN  = mean((tmpU-2*tmp+tmpD))/(delS^2);


% (2) W/O Common Random Numbers
s =  tmp1*exp(tmp2*z1);
tmp = exp(-interest*time)*max(s-strike,0).*(s<barrier);

s = tmp1U*exp(tmp2*z2);
tmpU = exp(-interest*time)*max(s-strike,0).*(s<barrier);

s = tmp1D*exp(tmp2*z3);
tmpD = exp(-interest*time)*max(s-strike,0).*(s<barrier);

delC_tilde1 = mean((tmpU-tmpD))/(2*delS);

gam_tilde  = mean((tmpU-2*tmp+tmpD))/(delS^2);


% (3) Pathwise Estimator

% Not Continuous
delC_pathwise = NaN;
gam_pathwise = NaN;

% (4) Likelihood ratio
s =  tmp1*exp(tmp2*z1);
scoreH = exp(-interest*time)*max(s-strike,0).*(s<barrier).*z1/(spot*volatility*sqrt(time));
delC_likelihood = mean(scoreH);

scoreG = exp(-interest*time)*max(s-strike,0).*(s<barrier).*((z1.^2-1)/(spot^2*volatility^2*time)-z1/(spot^2*volatility*sqrt(time)));
gam_likelihood = mean(scoreG);



%             Vega 

tmp1U  =  spot *exp((interest-dividend-(volatility+delSig)^2/2)*time);
tmp1D  =  spot *exp((interest-dividend-(volatility-delSig)^2/2)*time);

tmp2U = (volatility+delSig)*sqrt(time);
tmp2D = (volatility-delSig)*sqrt(time);

% (1) W/ Common Random Numbers
s = tmp1U*exp(tmp2U*z1);
tmpU = exp(-interest*time)*max(s-strike,0).*(s<barrier);

s = tmp1D*exp(tmp2D*z1);
tmpD = exp(-interest*time)*max(s-strike,0).*(s<barrier);

vega_CRN = mean((tmpU-tmpD)/(2*delSig));

% (2) W/O Common Random Numbers
s = tmp1U*exp(tmp2U*z1);
tmpU = exp(-interest*time)*max(s-strike,0).*(s<barrier);

s = tmp1D*exp(tmp2D*z2);
tmpD = exp(-interest*time)*max(s-strike,0).*(s<barrier);

vega_tilde = mean((tmpU-tmpD)/(2*delSig));

% (3) Pathwise Estimator
% Not continuous
vega_pathwise = NaN;

% (4) Likelihood Ratio

s =  tmp1*exp(tmp2*z1);
scoreH = exp(-interest*time)*max(s-strike,0).*(s<barrier).*((z1.^2-1)/volatility - z1*sqrt(time));
vega_likelihood = mean(scoreH);




fprintf("The value of Delta is \n")
disp([ ' W/ CRN: ' num2str(delC_CRN1),  ' W/O CRN: ', num2str(delC_tilde1), ' Pathwise Estimator: ', num2str(delC_pathwise), ' Likelihood Ratio: ', num2str(delC_likelihood)]);

fprintf("The value of Gamma is \n")
disp([' W/ CRN: ' num2str(gam_CRN),  ' W/O CRN: ', num2str(gam_tilde), ' Pathwise Estimator: ', num2str(gam_pathwise), ' Likelihood Ratio: ', num2str(gam_likelihood)]);

fprintf("The value of Vega is \n")
disp([' W/ CRN: ' num2str(vega_CRN),  ' W/O CRN: ', num2str(vega_tilde), ' Pathwise Estimator: ', num2str(vega_pathwise), ' Likelihood Ratio: ', num2str(vega_likelihood)]);





fprintf("Average Option Call \n")

% Simulation

z1 = randn(Simulations,monitoringTimes);
z2 = randn(Simulations,monitoringTimes);
z3 = randn(Simulations,monitoringTimes);
dt = time/monitoringTimes;
si = zeros(Simulations,monitoringTimes);
siD = si;
siU = si;

%         Delta & Gamma

tmp2 = volatility*sqrt(dt);
step = exp((interest-dividend-volatility^2/2)*dt);


% (1) W/ Common Random Numbers

si(:,1) = spot*step*exp(tmp2*z1(:,1));
siU(:,1) = (spot+delS)*step*exp(tmp2*z1(:,1));
siD(:,1) = (spot-delS)*step*exp(tmp2*z1(:,1));

for i=2:monitoringTimes
    si(:,i) = si(:,i-1)*step.*exp(tmp2*z1(:,i));
    siU(:,i) = siU(:,i-1)*step.*exp(tmp2*z1(:,i));
    siD(:,i) = siD(:,i-1)*step.*exp(tmp2*z1(:,i));   
end
means = mean(si,2);
meansU = mean(siU,2);
meansD = mean(siD,2);

temp1= exp(-interest*time)*max(means-strike,0);
temp1U = exp(-interest*time)*max(meansU-strike,0);
temp1D = exp(-interest*time)*max(meansD-strike,0);

delC_CRN1 = mean((temp1U-temp1D))/(2*delS);

gam_CRN  = mean((temp1U-2*temp1+temp1D))/(delS^2);


% (2) W/O Common Random Numbers
si(:,1) = spot*step*exp(tmp2*z1(:,1));
siU(:,1) = (spot+delS)*step*exp(tmp2*z2(:,1));
siD(:,1) = (spot-delS)*step*exp(tmp2*z3(:,1));

for i=2:monitoringTimes
    si(:,i) = si(:,i-1)*step.*exp(tmp2*z1(:,i));
    siU(:,i) = siU(:,i-1)*step.*exp(tmp2*z2(:,i));
    siD(:,i) = siD(:,i-1)*step.*exp(tmp2*z3(:,i));   
end
means = mean(si,2);
meansU = mean(siU,2);
meansD = mean(siD,2);

temp1 = exp(-interest*time)*max(means-strike,0);
temp1U = exp(-interest*time)*max(meansU-strike,0);
temp1D = exp(-interest*time)*max(meansD-strike,0);
delC_tilde1 = mean((temp1U-temp1D))/(2*delS);
%
gam_tilde  = mean((temp1U-2*temp1+temp1D))/(delS^2);


% (3) Pathwise Estimator
dhds0 = exp(-interest*time)*(means>strike).*means/spot;
delC_pathwise = mean(dhds0);

% Delta not continuous
gam_pathwise = NaN;

% (4) Likelihood ratio

scoreH = exp(-interest*time)*max(means-strike,0).*z1(:,1)/(spot*volatility*sqrt(dt));
delC_likelihood = mean(scoreH);

scoreG = ((z1(:,1).^2-1)/(spot^2*volatility^2*dt)-z1(:,1)/(spot^2*volatility*sqrt(dt)));
scoreG = exp(-interest*time)*max(means-strike,0).*scoreG;
gam_likelihood = mean(scoreG);


%             Vega Part C


tmp1U  =  exp((interest-dividend-(volatility+delSig)^2/2)*dt);
tmp1D  =  exp((interest-dividend-(volatility-delSig)^2/2)*dt);

tmp2U = (volatility+delSig)*sqrt(dt);
tmp2D = (volatility-delSig)*sqrt(dt);

% (1) W/ Common Random Numbers

siU(:,1) = spot*tmp1U*exp(tmp2U*z1(:,1));
siD(:,1) = spot*tmp1D*exp(tmp2D*z1(:,1));


for i=2:monitoringTimes
   
    siU(:,i) = siU(:,i-1)*tmp1U.*exp(tmp2U*z1(:,i));
    siD(:,i) = siD(:,i-1)*tmp1D.*exp(tmp2D*z1(:,i));   
end

meansU = mean(siU,2);
meansD = mean(siD,2);

tmpU = exp(-interest*time)*max(meansU-strike,0);
tmpD = exp(-interest*time)*max(meansD-strike,0);
%
vega_CRN = mean((tmpU-tmpD)/(2*delSig));

% (2) W/O Common Random Numbers
siU(:,1) = spot*tmp1U*exp(tmp2U*z2(:,1));
siD(:,1) = spot*tmp1D*exp(tmp2D*z3(:,1));

for i=2:monitoringTimes
   
    siU(:,i) = siU(:,i-1)*tmp1U.*exp(tmp2U*z2(:,i));
    siD(:,i) = siD(:,i-1)*tmp1D.*exp(tmp2D*z3(:,i));   
end

meansU = mean(siU,2);
meansD = mean(siD,2);

tmpU = exp(-interest*time)*max(meansU-strike,0);
tmpD = exp(-interest*time)*max(meansD-strike,0);
%
vega_tilde = mean((tmpU-tmpD)/(2*delSig));

% (3) Pathwise Estimator
dhdsig = 0;
for i=1:monitoringTimes
   dhdsig = dhdsig+si(:,i).*(log(si(:,i)/spot)-(interest-dividend+(1/2)*volatility^2)*i*dt);
end
dhdsig = exp(-interest*time)*(means>strike)*1/(monitoringTimes*volatility).*dhdsig;
vega_pathwise = mean(dhdsig);


% (4) Likelihood Ratio
scoreH = 0;
for i=1:monitoringTimes
    scoreH = scoreH+(z1(:,i).^2-1)/volatility-z1(:,i)*sqrt(dt);
end

scoreH = exp(-interest*time)*max(means-strike,0).*scoreH;
vega_likelihood = mean(scoreH);


fprintf("The value of Delta is \n")
disp([ ' W/ CRN: ' num2str(delC_CRN1),  ' W/O CRN: ', num2str(delC_tilde1), ' Pathwise Estimator: ', num2str(delC_pathwise), ' Likelihood Ratio: ', num2str(delC_likelihood)]);

fprintf("The value of Gamma is \n")
disp([' W/ CRN: ' num2str(gam_CRN),  ' W/O CRN: ', num2str(gam_tilde), ' Pathwise Estimator: ', num2str(gam_pathwise), ' Likelihood Ratio: ', num2str(gam_likelihood)]);

fprintf("The value of Vega is \n")
disp([' W/ CRN: ' num2str(vega_CRN),  ' W/O CRN: ', num2str(vega_tilde), ' Pathwise Estimator: ', num2str(vega_pathwise), ' Likelihood Ratio: ', num2str(vega_likelihood)]);



% Output
% Digital Call 
% The value of Delta is 
%  W/ CRN: 0.011991 W/O CRN: 0.015305 Pathwise Estimator: NaN Likelihood Ratio: 0.011689
% The value of Gamma is 
%  W/ CRN: 0.013451 W/O CRN: 0.33109 Pathwise Estimator: NaN Likelihood Ratio: 0.00011052
% The value of Vega is 
%  W/ CRN: 0.31994 W/O CRN: -0.67063 Pathwise Estimator: NaN Likelihood Ratio: 0.30945
% Up-and-Out Call 
% The value of Delta is 
%  W/ CRN: 0.029167 W/O CRN: 0.054532 Pathwise Estimator: NaN Likelihood Ratio: 0.025707
% The value of Gamma is 
%  W/ CRN: -0.054397 W/O CRN: -0.0086566 Pathwise Estimator: NaN Likelihood Ratio: -0.00057892
% The value of Vega is 
%  W/ CRN: -1.393 W/O CRN: -2.4876 Pathwise Estimator: NaN Likelihood Ratio: -1.621
% Average Option Call 
% The value of Delta is 
%  W/ CRN: 0.24958 W/O CRN: 0.29729 Pathwise Estimator: 0.2496 Likelihood Ratio: 0.2491
% The value of Gamma is 
%  W/ CRN: 0.018636 W/O CRN: 1.8552 Pathwise Estimator: NaN Likelihood Ratio: 0.018028
% The value of Vega is 
%  W/ CRN: 19.5517 W/O CRN: 24.2622 Pathwise Estimator: 19.5522 Likelihood Ratio: 19.6868