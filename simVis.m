close all;
runSim = 1;
runs = 2 -1; %min of 2 (includes start and finish)
angles = 0:pi/(2*runs):pi/2;
r0 = 2;
g = 9.81;
vx = zeros(runs+1,1);
vy = zeros(runs+1,1);
val = zeros(runs+1,4);
    vyi = sqrt( 2*g*(r0-Thf0(2)*sin(Thf0(1))) );
    dti = -vyi*cos(Thf0(1));
    dri = -vyi*sin(Thf0(1));
for i = 1:1
    qIn = [Thf0(1); Thf0(2); dti; dri];
    disp(qIn);
    %kIn = [25;50;500;200;50;10];
    kIn = 50;
    %size(kIn,1)
    simOut = sim('simSLIPModel');
    q = get(simOut, 'q');
    q = squeeze(q.signals.values);
    signalSize = size(q);
    rows = max(signalSize);
    signalSizeS = min(signalSize)/2+1;
    the1 = q(end,1);
    r1 = q(end,2);
    Dthe1 = q(end,3);
    Dr1 = q(end,4);
    val(i,:) = [the1,r1,Dthe1,Dr1];
    vx(i) = Dr1*cos(the1) - Dthe1*r1*sin(the1);
    vy(i) = Dr1*sin(the1) + Dthe1*r1*cos(the1);
    model_sim;
end