clc;
Spot = 1;
maturity = 35;
sig = 0.23;
mu = 0.07;
m = 35;
dt = maturity/m;

S = zeros(m,1);
t = zeros(m,1);

disp('simulating index values')
S = ones(10000,m);
%S(1) = Spot;
t(1) = 0;
T = maturity;
for n = 1:10000
    for i = 2:m
        
        T = T-dt;
        z = randn;
        
        S(n,i) = S(n,i-1) * exp((mu-sig*sig/2)*dt + sig*sqrt(dt)*z);
        t(i) = maturity - T;
        
    end
end

ret_index = ones(10000,m);
for n = 1:10000
    for i = 2:m
        
        ret_index(n,i) = S(n,i)/S(n,i-1)-1;
    end
end
disp('Sample plot of index for 20 simulations')
plot(1:35,S(1:20,:))
%%Portfolio S1
disp('--------------------------------------')

disp('Portfolio S1')
P1 = zeros(1,m);
P1(1,1) = 10000;

for j = 1:m-1
    P1(1,j+1) = P1(1,j)*1.03 + 10000;
end
P1(1,m) = P1(1,m);
disp('Value of Portfolio S1 at the end of 35 years period')
P1(1,35)
disp('--------------------------------------')



%Portfolio S2
disp('--------------------------------------')
disp('Portfolio S2')
P2 = zeros(n,m);
P2(:,1) = 10000;

for j=2:m
    for i = 1:10000
        P2(i,j) = (1+(ret_index(i,j)))*P2(i,j-1) + 10000;
    end
end


%Calculate probability of final portfolio value exceeding $1M
for j=1:10000
    if P2(j,m)>1000000
        P_2M(j,1)=1;
    end
end


Prob_avg = mean(P_2M(:,1))
Prob_std = std(P_2M(:,1))
P2_Std_Err_Prob = Prob_std/sqrt(10000)
P2_avg = mean(P2(:,m))
P2_med = median(P2(:,m))
P2_Std_Err_Avg = std(P2(:,m))/sqrt(10000)
P2_stddev = std(P2(:,m))
P2_perc10 = prctile(P2(:,m),10)
P2_perc50 = prctile(P2(:,m),50)
P2_perc90 = prctile(P2(:,m),90)
disp('--------------------------------------')

%Strategy 3
disp('--------------------------------------')
disp('Portfolio S3')
P3 = zeros(10000,m);
P3(:,1) = 10000;
P3_1M = zeros(n,1);

for j=2:m
    for i = 1:10000
        P3(i,j) = ((1+(ret_index(i,j)))*P3(i,j-1)*0.5) + (0.5*P3(i,j-1)*1.03) + 10000;
    end
end

%Calculate probability of final portfolio value exceeding $1M
for j=1:10000
    if P3(j,m)>1000000
        P3_1M(j,1)=1;
    end
end

Prob_avg3 = mean(P3_1M(:,1))
Prob_std3 = std(P3_1M(:,1))
P2_Std_Err_Prob = Prob_std3/sqrt(10000)
P3_avg = mean(P3(:,m))
P3_med = median(P3(:,m))
P3_stddev = std(P3(:,m))
P3_Std_Err_Avg = P3_stddev/sqrt(10000)
P3_perc10 = prctile(P3(:,m),10)
P3_perc50 = prctile(P3(:,m),50)
P3_perc90 = prctile(P3(:,m),90)

disp('--------------------------------------')

% Dynamic rebalancing strategy
disp('--------------------------------------')
disp('Portfolio S4: Dynamic rebalancing Strategy')
P4 = zeros(n,m);
P4(:,1) = 10000;
P4_1M = zeros(n,1);
sharpe = zeros(1,m);

for j=3:m
    %Growth optimal portfolio
    sharpe(1,j) = (mean(ret_index(:,j-1)-0.03))/var(ret_index(:,j-1));
end
sharpe(1,1)=0.5;
sharpe(1,2)=0.5;
for j=2:m
    for i = 1:n
        P4(i,j) = ((1+(ret_index(i,j)))*P4(i,j-1)*sharpe(1,j-1)) + ((1-sharpe(1,j-1))*P4(i,j-1)*1.03) + 10000;
    end
end


%Calculate probability of final portfolio value exceeding $1M
for j=1:n
    if P4(j,m)>1000000
        P4_1M(j,1)=1;
    end
end

disp('b* calculated every year and portfolio re-balanced on calculated b*')
Prob_avg4 = mean(P4_1M(:,1))
Prob_std4 = std(P4_1M(:,1))
%Strategy 4 - Portfolio Value
P4_Std_Err_Prob = Prob_std4/sqrt(10000)
P4_avg = mean(P4(:,m))
P4_med = median(P4(:,m))
P4_stddev = std(P4(:,m))
P4_Std_Err_Avg = P4_stddev/sqrt(10000)
disp('--------------------------------------')

disp('--------------------------------------')
disp('Constant ratio Rebalancing Strategy')
%Constant Ratio Rebalancing Strategy
P5 = zeros(n,m);
P5(:,1) = 10000;
P5_1M = zeros(n,1);
const_sharpe = zeros(1,m);


%Growth optimal portfolio

const_sharpe = (0.07-0.03)/(0.23^2)

sharpe(1,1)=0.5;
sharpe(1,2)=0.5;
for j=2:m
    for i = 1:n
        P5(i,j) = ((1+(ret_index(i,j)))*P5(i,j-1)*const_sharpe) + ((1-const_sharpe)*P5(i,j-1)*1.03) + 10000;
    end
end


%Calculate probability of final portfolio value exceeding $1M
for j=1:n
    if P5(j,m)>1000000
        P5_1M(j,1)=1;
    end
end

disp('b* calculated at initial time based on mu=7, rf=3, sigma=0.23')
disp('Rebalancing done every year on constant b*')
Prob_avg5 = mean(P5_1M(:,1))
Prob_std5 = std(P5_1M(:,1))

P5_Std_Err_Prob = Prob_std5/sqrt(10000)
P5_avg = mean(P5(:,m))
P5_med = median(P5(:,m))
P5_stddev = std(P5(:,m))
P5_Std_Err_Avg = P5_stddev/sqrt(10000)

disp('--------------------------------------')

%Part C
disp('--------------------------------------')
disp('Part C. Portfolio S2 Calculating Probability of hitting 20000 before 5000 for portfolio ')

P2new = zeros(10000,m);
P2new(:,1) = 10000;
P2new_1N = zeros(10000,1);
P3new = zeros(10000,m);
P3new(:,1) = 10000;
P3new_1N = zeros(10000,1);


for j=2:m
    for i = 1:10000
        P2new(i,j) = (1+(ret_index(i,j)))*P2new(i,j-1);
    end
end
counter=0;
for j=1:10000
    index_5000=min(find(P2new(j,:)<=5000));
    index_20000=min(find(P2new(j,:)>=20000));
    if isempty(index_20000) && isempty(index_5000)
        counter = counter+1;
        P2new_IN(j,1)=0;
    elseif isempty(index_20000)
        P2new_IN(j,1)=0;
    elseif isempty(index_5000)
        P2new_IN(j,1)=1;
    else
        if index_5000<index_20000
            P2new_IN(j,1)=0;
        else P2new_IN(j,1)=1;
        end
    end
    
end
sum1=sum(P2new_IN)
Prob2_new_avg = sum1/(n-counter)
Std_Prob2_new = std(P2new_IN)/sqrt(n-counter)
disp('--------------------------------------')


%Strategy 3 - new
disp('--------------------------------------')
disp('Part C. Portfolio S3 Calculating Probability of hitting 20000 before 5000 for portfolio ')
for j=2:m
    for i = 1:10000
        P3new(i,j) = ((1+(ret_index(i,j)))*P3new(i,j-1)*0.5) + (0.5*P3new(i,j-1)*1.03);
    end
end
counter=0;
for j=1:n
    index_5000_3=min(find(P3new(j,:)<=5000));
    index_20000_3=min(find(P3new(j,:)>=20000));
    if isempty(index_20000) && isempty(index_5000)
        counter = counter+1;
        P3new_IN(j,1)=0;
    elseif isempty(index_20000_3)
        P3new_IN(j,1)=0;
    elseif isempty(index_5000_3)
        P3new_IN(j,1)=1;
    else
        if index_5000_3<index_20000_3
            P3new_IN(j,1)=0;
        else P3new_IN(j,1)=1;
        end
    end
end

sum2=sum(P3new_IN);
Prob3_new_avg = sum2/(n-counter-1)
Std_Prob3_new = std(P3new_IN)/sqrt(n-counter-1)
disp('--------------------------------------')

%Output
% simulating index values
% Sample plot of index for 20 simulations
% --------------------------------------
% Portfolio S1
% Value of Portfolio S1 at the end of 35 years period
% 
% ans =
% 
%    6.0462e+05
% 
% --------------------------------------
% --------------------------------------
% Portfolio S2
% 
% Prob_avg =
% 
%     0.4411
% 
% 
% Prob_std =
% 
%     0.4965
% 
% 
% P2_Std_Err_Prob =
% 
%     0.0050
% 
% 
% P2_avg =
% 
%    1.4524e+06
% 
% 
% P2_med =
% 
%    8.6890e+05
% 
% 
% P2_Std_Err_Avg =
% 
%    1.9969e+04
% 
% 
% P2_stddev =
% 
%    1.9969e+06
% 
% 
% P2_perc10 =
% 
%    3.0132e+05
% 
% 
% P2_perc50 =
% 
%    8.6890e+05
% 
% 
% P2_perc90 =
% 
%    3.0153e+06
% 
% --------------------------------------
% --------------------------------------
% Portfolio S3
% 
% Prob_avg3 =
% 
%     0.3306
% 
% 
% Prob_std3 =
% 
%     0.4705
% 
% 
% P2_Std_Err_Prob =
% 
%     0.0047
% 
% 
% P3_avg =
% 
%    9.2544e+05
% 
% 
% P3_med =
% 
%    8.1087e+05
% 
% 
% P3_stddev =
% 
%    4.8153e+05
% 
% 
% P3_Std_Err_Avg =
% 
%    4.8153e+03
% 
% 
% P3_perc10 =
% 
%    4.6660e+05
% 
% 
% P3_perc50 =
% 
%    8.1087e+05
% 
% 
% P3_perc90 =
% 
%    1.5104e+06
% 
% --------------------------------------
% --------------------------------------
% Portfolio S4: Dynamic rebalancing Strategy
% b* calculated every year and portfolio re-balanced on calculated b*
% 
% Prob_avg4 =
% 
%     0.4027
% 
% 
% Prob_std4 =
% 
%     0.4905
% 
% 
% P4_Std_Err_Prob =
% 
%     0.0049
% 
% 
% P4_avg =
% 
%    1.0843e+06
% 
% 
% P4_med =
% 
%    8.4879e+05
% 
% 
% P4_stddev =
% 
%    8.3666e+05
% 
% 
% P4_Std_Err_Avg =
% 
%    8.3666e+03
% 
% --------------------------------------
% --------------------------------------
% Constant ratio Rebalancing Strategy
% 
% const_sharpe =
% 
%     0.7561
% 
% b* calculated at initial time based on mu=7, rf=3, sigma=0.23
% Rebalancing done every year on constant b*
% 
% Prob_avg5 =
% 
%     0.4191
% 
% 
% Prob_std5 =
% 
%     0.4934
% 
% 
% P5_Std_Err_Prob =
% 
%     0.0049
% 
% 
% P5_avg =
% 
%    1.1625e+06
% 
% 
% P5_med =
% 
%    8.6328e+05
% 
% 
% P5_stddev =
% 
%    1.0356e+06
% 
% 
% P5_Std_Err_Avg =
% 
%    1.0356e+04
% 
% --------------------------------------
% --------------------------------------
% Part C. Portfolio S2 Calculating Probability of hitting 20000 before 5000 for portfolio 
% 
% sum1 =
% 
%         7742
% 
% 
% Prob2_new_avg =
% 
%     0.7983
% 
% 
% Std_Prob2_new =
% 
%     0.0042
% 
% --------------------------------------
% --------------------------------------
% Part C. Portfolio S3 Calculating Probability of hitting 20000 before 5000 for portfolio 
% 
% Prob3_new_avg =
% 
%     0.9191
% 
% 
% Std_Prob3_new =
% 
%     0.0027
% 
% --------------------------------------