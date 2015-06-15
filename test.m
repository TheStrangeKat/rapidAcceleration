r0 = 2;
x0 = [pi/2-0.1;1.5];
LB = [0;r0/8];
UB = [pi/2;r0];
Aineq = []; Bineq = []; Aeq = []; Beq = [];
nonlcon = @heightConstraint;
%optimization options
options = optimset('display','Iter','MaxFunEvals',15000,'MaxIter',7000,'UseParallel','always','Algorithm','interior-point',...
    'TolCon',1e-6,'FinDiffRelStep',5e-4, 'TolFun', 1e-6, 'TolX', 1e-9, 'LargeScale', 'off','diffmaxchange',1,'diffminchange',5e-4);
optQ = fmincon(@funkyFunc,x0,Aineq,Bineq,Aeq,Beq,LB,UB,@heightConstraint,options);
disp(funkyFunc(optQ));
disp(heightConstraint(optQ));