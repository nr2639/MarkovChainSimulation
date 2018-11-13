function evalU = evaluateU(q, rho)
if 0<=q(1) && q(1)<=1&&0<=q(2)&&q(2)<=1
    evalU = 0.5;
elseif -1<=q(1)&&q(1)<=0 &&-1<=q(2)&&q(2)<=0
    evalU = 0.5;
else
    evalU = 0;
end
end