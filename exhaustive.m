function [vx] = exhaustive(q)
%exhaustive is used for an exhaustive search in the configuration space.
%   This function returns the final x velocity of the model, or -1 if the
%   function does not meet constraints.
    r0 = 2;
    g = 9.81;
    vyi = sqrt( 2*g*(r0-q(2)*sin(q(1))) );
    dti = -vyi*cos(q(1));
    dri = -vyi*sin(q(1));
    qIn = [q(1); q(2); dti; dri];
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
    outConstraint = (sin(the1) + sign(vy)*(vy^2)/(2*9.81*r0) - 1);
    if(outConstraint>0)    
       vx = (Dr1*cos(the1) - Dthe1*r1*sin(the1));
    else
       vx = -1;
    end
end

