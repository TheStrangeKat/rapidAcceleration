function [outConstraint,word] = heightConstraint(q)
%heightConstraint - Returns a value describing the height constraint of the
%function. This constraint is satisfied when the returned value is greater
%than 0.
%   Detailed explanation goes here
    r0 = 2;
    g = 9.81;
    vyi = sqrt( 2*g*(r0-q(2)*sin(q(1))) );
    dti = -vyi*cos(q(1));
    dri = -vyi*sin(q(1));
    qIn = [q(1); q(2); dti; dri];
%     disp(qIn);
    %kIn = [25;50;500;200;50;10];
    kIn = 50;
    options = simset('SrcWorkspace','current');
    simOut = sim('simSLIPModel',[],options);
    q = get(simOut, 'q');
    q = squeeze(q.signals.values);
    the1 = q(end,1);
    r1 = q(end,2);
    Dthe1 = q(end,3);
    Dr1 = q(end,4);
    vy = Dr1*sin(the1) + Dthe1*r1*cos(the1);
%     disp(vy);
%     disp([the1;r1;Dthe1;Dr1]);
    outConstraint = -(sin(the1) + sign(vy)*(vy^2)/(2*9.81*r0) - 1); %const should be 1, not 1.1 - represents multiple of height.
    word = 0;
end
