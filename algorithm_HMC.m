function qNew = algorithm_HMC(qCurrent, e, L,rho)
%
q = qCurrent;
d = length(q);
p = randn(1,d);
pCurrent = p;
% At the beginning, make a half step for momentum
gradU = 0;
p = p - (e)*gradU;
% alternate full steps for p & q
for i = 1:L
    q = q + e*p; % make a full step for the position
    if i < L
        % make a full step for the momentum, except at the end of the
        % trajectory
        gradU = 0;
        p = p - e*gradU;
    else
        % make a half step for the momentum at the end of the
        % trajectory
        p = p - (e/2)*gradU;
        
    end
end

% To make the proposal symmetric, negate momentum at the end of trajectory
p = -p;
% Evaluate potential & kinetic energies to start and end the trajectory
U = evaluateU(qCurrent);
K = sum(pCurrent.^2)/2;
% proposed U
U_tilde = evaluateU(q);
% proposed K
K_tilde = sum(p.^2)/2;
% accept/reject the state at the end of trajectory.
% i.e. returning either the position at the end of the trajectory or the
% initial position
if (U + K) <= (U_tilde + K_tilde)
    qNew = q;
else
    qNew = qCurrent;
end
end