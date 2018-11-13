
clear all; 
close all;
clc;
m = 1000;
qArray = zeros(m+1,2);
qArray(1,:) = [0.33 0.40];
e = 0.05;
L = 20;
rho = 0.00;%%Unused

for i = 2:m+1
qArray(i,:) = algorithm_HMC(qArray(i-1,:), e, L,rho);
end
q1Array = qArray(:,1);
q2Array = qArray(:,2);
figure(1);
hold on;
plot(q1Array, q2Array, 'ro');
%hold off;
title('Position Coordinates');
%hold on;

%% Traceplots

figure(2);
%Traceplots
subplot(2,1,1);
plot(q1Array, 'r');
xlabel('I_k');
ylabel('q_1')
title('q_1');
axis tight;
hold on;

subplot(2,1,2);
plot(q2Array, 'r');
xlabel('I_k');
ylabel('q_2')
title('q_2');
axis tight;
hold on;


%% Running Average
figure(3);

q1RunningAvg = zeros(m,1);
q2RunningAvg = zeros(m,1);

q1RunningAvg(1) = q1Array(1);
q2RunningAvg(1) = q2Array(1);
% recursive calculation of averages
for k = 2:m
    q1RunningAvg(k) = (q1RunningAvg(k-1)*(k-1)+q1Array(k))/k;
    q2RunningAvg(k) = (q2RunningAvg(k-1)*(k-1)+q2Array(k))/k;
end
subplot(2,1,1);
plot(q1RunningAvg, 'r');
xlabel('I_k');
ylabel('Running Average of q_1')
axis tight;
hold on;

subplot(2,1,2);
plot(q2RunningAvg, 'r');
xlabel('I_k');
ylabel('Running Average of q_2')
axis tight;
hold on;
