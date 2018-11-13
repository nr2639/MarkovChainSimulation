function exampleMarkovChain(initialState, T,P,Lamda )

%% simulation of CTMC 

P = [  0.3125,    0.1936 ,   0.0835,    0.4105;
    0.3156 ,   0.3319  ,  0.2433 ,   0.1093;
    0.2307,    0.2130   , 0.3121 ,   0.2443;
    0.4989   , 0.1075   , 0.2250  ,  0.1686;]

Lamda = [0.5,0.75,1,1.25];
T= 10;
initialState = 1;
nStates = length(P);





%P_unNormalized = rand(nStates, nStates);

%sum_rows = sum(P_unNormalized,2);
%sum_rows_extended = repmat(sum_rows,1,nStates);



% normlized to get transition probability matrix
%P = P_unNormalized./sum_rows_extended;

disp('Transition Probability Matrix');
disp(P);

x = zeros(nStates+1,1);


current_state = initialState;
x = [];
%initial time t =0
t=0;

while t < T
    
    %which state at current time 
    i = current_state;
    %disp(['current state:', num2str(i)]);
    
    %pick the corresponding row
    p = P(i,:);
    
    u = rand;
    %disp(u);
    hold_time = poissrnd(Lamda(i));
    t = t + hold_time;
    
    for i = 1:nStates
        if i == 1
            if u <= p(i)
                j = i;
                break;
            end
        elseif i == nStates
            j = i;
            break;
        else
            if ((u > sum(p(1:i-1))) && (u <= sum(p(1:i))))
                j = i;
                break;
            end    
        end
    end
    x = [x j];
    current_state = j;
    disp(['Switching time',num2str(t), ' transiting from state ', num2str(i), ' to state ',num2str(j)]);
   
    %disp(' ');
    title(['Markov Chain with ', num2str(nStates), ' states']); 
    plot(x, 'o-');
    pause(hold_time);
    
end

% Output:
% 
% P =
% 
%     0.3125    0.1936    0.0835    0.4105
%     0.3156    0.3319    0.2433    0.1093
%     0.2307    0.2130    0.3121    0.2443
%     0.4989    0.1075    0.2250    0.1686
% 
% Transition Probability Matrix
%     0.3125    0.1936    0.0835    0.4105
%     0.3156    0.3319    0.2433    0.1093
%     0.2307    0.2130    0.3121    0.2443
%     0.4989    0.1075    0.2250    0.1686
% 
% Switching time1 transiting from state 4 to state 4
% Switching time2 transiting from state 4 to state 4
% Switching time5 transiting from state 1 to state 1
% Switching time6 transiting from state 4 to state 4
% Switching time6 transiting from state 3 to state 3
% Switching time10 transiting from state 2 to state 2 


    
    
    
    

