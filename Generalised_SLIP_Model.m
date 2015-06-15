%% Generalised_SLIP_Model
%
%   Calculates the Equations of Motion for a N-Link Pendulum
%
%   NF Steenkamp
%   University of Cape Town
%   Electrical Engineering Department
%   22 May 2015 - Last Updated 25 May 2015
%% 
Passive = false;
Actuated = true;

%% Settings
joints = [Passive];
generateEulerLagrange = false;
generateEngergyFunction = true;
generateFlight = true;

%% Declare Variables-------------------------------------------------------
noLinks = max(size(joints));
disp('...');
disp('Declaring Variables and calculating positions/velocities...');

% General variables
q = sym('the', [noLinks 1]);    % Generalized Coordinates
sym(q,'real');
Dq = sym('Dthe', [noLinks 1]);  % First Derivative of q
D2q = sym('D2the', [noLinks 1]);% Second Derivative of q

r = sym('r', [noLinks 1]);      % Length of links
r0 = sym('r0', [noLinks 1]);    % Rest condition of spring
sym(r,'real');
Dr = sym('Dr', [noLinks 1]);    % First Derivative of r
D2r = sym('D2r', [noLinks 1]);  % Second Derivative of r
x = r;%/2;                        % COM of links (relative to frame)
% 
% Dr(:) = 0;
% D2r(:) = 0;

tau = sym('tau', [noLinks 1]);  % Torques (control inputs)
mass = sym('mass', [noLinks 1]);% Link Masses
sym(mass,'real');
%I = sym('I',[noLinks 1]);       % Mass Moment of Inertia for Links
I = 0;%mass.*(r.^2)/12;            % Mass Moment of Inertia for Links
R0 = cell(noLinks,1);           % Rotation Matrices
T = sym('T', [noLinks 1]);      % Kinectic Energy
U = sym('U', [noLinks 1]);      % Potential Energy
EL = sym('EL', [noLinks 1]);    % Euler-Lagrange
syms g k_spring Ttot Utot L;

disp('...Manipulator Equations...');
% Manipulator Equations M*D2q + C*Dq + G = Bu
M = sym('M', 2*noLinks); %Inertial Matrix
sym(M,'real');
C = sym('C', 2*noLinks); % Coriolis and Centripetal Terms
sym(C,'real');
G = sym('G', [2*noLinks 1]); % Gravitational Terms
sym(G,'real');
B = sym('B', 2*noLinks); % Torque Terms
sym(B,'real');

% Splitting variables into actuated and unactuated
noAct = sum(joints==Actuated);
noPas = sum(joints==Passive);
qp = q(joints==Passive);
qa = q(joints==Actuated);
Dqp = Dq(joints==Passive);
Dqa = Dq(joints==Actuated);

% Adjusting for actuated joints
tau(joints==Passive) = sym(0);

%--------------------------------------------------------------------------

%% Rotation Matrices-------------------------------------------------------
disp('...Rotation Matrices...');
% Rotation Matrices - Relative to frame
for i = 1:noLinks
    var = [cos(q(i)),-sin(q(i)),0;
           sin(q(i)), cos(q(i)),0;
           0        , 0        ,1];
         R0{i} = var;
end
% Rotation Matrices - Relative to inertial frame
for i = 2:noLinks
    R0{i} = simplify(R0{i-1}*R0{i});
end

%--------------------------------------------------------------------------

%% Positions---------------------------------------------------------------
% position mass relative to frame origin in inertial frame coords
pm = sym('pm', [3,noLinks]);
pm(1,:)=x.';
pm(2:3,:)=0;
for i = 1:noLinks
    pm(:,i) = R0{i}*pm(:,i);
end

% position of the next (n+1) frame origin relative to frame origin in the
% inertial frame coords
po = sym('po', [3,noLinks]);
po(1,:)=r.';
po(2:3,:)=0;
for i = 1:noLinks
    po(:,i) = R0{i}*po(:,i);
end

% position of mass relative to inertial origin
xm = sym('ym', [3,noLinks]);
xm(1,:)=r.';
xm(2:3,:)=0;
for i = 1:noLinks
    j = 1;
    xm(:,i) = pm(:,i);
    while (j<i)
        xm(:,i) = xm(:,i) + po(:,j);
        j = j+1;
    end
end

%--------------------------------------------------------------------------

%% Velocities--------------------------------------------------------------
% Angular Velocities 
omega = sym('omega', [3 noLinks]);
omega(1:2,:)=0;
omega(3,1) = Dq(1);
for i = 2:noLinks
    omega(3,i) = omega(3,i-1) + Dq(i);
end

% Rectilinear Velocities
v = sym('v', [3 noLinks]);
v(:,:) = 0;
for i = 1:noLinks
    j = 1;
    for k = 1:noLinks
        v(:,i) = v(:,i) + diff(pm(:,i),q(k))*Dq(k) + diff(pm(:,i),r(k))*Dr(k);
    end
    
    while (j<i)
        sym temp;
        temp = 0;
        for k = 1:noLinks
            temp = temp + diff(po(:,j),q(k))*Dq(k) + diff(po(:,j),r(k))*Dr(k);
        end        
        v(:,i) = v(:,i) + temp;
        j = j+1;
    end
end
v = simplify(v);

%--------------------------------------------------------------------------

%% Energy Equations--------------------------------------------------------
% Kinetic Energy
disp('Calculating Kinetic Energy...');
Ttot = 0;
for i = 1:noLinks
   T(i) = simplify(0.5*(v(:,i).')*mass(i)*v(:,i) + 0.5*(omega(:,i).')*I(i)*omega(:,i));
   Ttot = Ttot + T(i);
end
Ttot = simplify(Ttot);

% Potential Energy
disp('Calculating Potential Energy...');
Utot = 0;
for i = 1:noLinks
    j = 1;
    U(i) = pm(2,i);
    while (j<i)
        U(i) = U(i) + po(2,j);
        j = j+1;
    end
    U(i) = mass(i)*g*U(i) + 0.5*k_spring*(r(i)-r0(i)).^2;
    Utot = Utot + U(i);
end
Utot = simplify(Utot);

%% Full Coordinates
q = [q;r];
Dq = [Dq;Dr];
D2q = [D2q;D2r];

tau = [tau;0*tau];
joints = [joints;zeros(noLinks,1)];

%% Energy Function
if(generateEngergyFunction)
    energy = simplify(T+U);
    matlabFunction(energy, 'File', 'getEnergy','Vars', {[q;Dq],mass,r0,g,k_spring});   %Creating function
end

%% Lagrangian - (currently not used)
if(generateEulerLagrange)
    disp('Calculating Lagrangian Equations...');
    L = Ttot - Utot;
    EL = sym('EL', [size(q,1) 1]);    % Euler-Lagrange
    sym temp;
    for i = 1:size(q,1)
        EL(i) = diff(L,Dq(i));    % Euler-Lagrange
        temp = 0;
        for j = 1:size(q,1)
           temp = temp + diff(EL(i),q(j))*Dq(j) + diff(EL(i),Dq(j))*D2q(j);
        end
        EL(i) = temp;
        EL(i) = EL(i) - diff(L,q(i)) - tau(i);
        EL(i) = simplify(EL(i));
    end
end

%% Manipulator Equations---------------------------------------------------
disp('Calculating Manipulator Equations...');
% Inertial Matrix
for i = 1:2*noLinks
    for j = 1:2*noLinks
        M(i,j) = diff(diff(Ttot,Dq(i)),Dq(j));
    end
end
M = simplify(M);

% Coriolis Matrix
C(:,:) = 0;
for i = 1:2*noLinks
    for j = 1:2*noLinks
        for k = 1:2*noLinks
            %disp([int2str(i) ' ' int2str(j) ' ' int2str(k)]);
            C(i,j) = C(i,j)+ 0.5*(diff(M(i,j),q(k)) + diff(M(i,k),q(j)) - diff(M(j,k),q(i)))*Dq(k);
        end
    end
end
C = simplify(C);

% Gravitational Matrix
for i = 1:2*noLinks
    G(i) = diff(Utot,q(i));
end
G = simplify(G);
% Torque Matrix
B = eye(2*noLinks) - diag(~joints); % Torque Matrix

% Equations of motion
ddq = simplify(-M\(C*Dq + G - B*tau));      % Solving for D2q
tau = sym('tau', [2*noLinks 1]);              % Allowing for construction of function
matlabFunction(ddq, 'File', 'ddq_func','Vars', {[q;Dq],mass,r0,g,k_spring,tau});   %Creating function
%matlabFunction(ddq(noLinks/2+1:noLinks), 'File', 'ddr','Vars', {[q;Dq],mass,g,k_spring,r0,tau});   %Creating function
%tau(joints==Passive) = sym(0);              % Re-correcting for passive joints


%% Partial Feedback Linearization------------------------------------------
disp(' Calculation PFL Equations...');
% Splitting Matrices up into unactuated and actuaged equations
N11 = M(joints==Passive,joints==Passive);
N12 = M(joints==Passive,joints==Actuated);
N21 = M(joints==Actuated,joints==Passive);
N22 = M(joints==Actuated,joints==Actuated);
D1 = C(joints==Passive,:);
D2 = C(joints==Actuated,:);
L1 = G(joints==Passive,:);
L2 = G(joints==Actuated,:);


%Taskspace
% COM = (ym*mass)./sum(mass);
endEffector = sum(po(:,:),2);
length = sqrt(endEffector(1).^2 + endEffector(2).^2);
angle = atan2(endEffector(2),endEffector(1));
y = [length;q(2)];%something something taskspace...
y = y(y~=0);
clength = 2;% length(y);
matlabFunction(y, 'File', 'map2TaskSpace','Vars', {[q;Dq],mass,r,g});

if(~isempty(joints(joints==Actuated)))
    % Jacobian
    disp('Calculating Jacobians...');
    J1 = sym('J1', [clength noPas]);      % Jacobian of Passive Joints
    J2 = sym('J2', [clength noAct]);      % Jacobian of Actuated Joints
    for i = 1:noPas
       J1(:,i) = diff(y,qp(i)); 
    end
    J1 = simplify(J1);

    for i = 1:noAct
       J2(:,i) = diff(y,qa(i)); 
    end
    J2 = simplify(J2);
    J = [J1,J2];
    Jbar = simplify(J2 - J1/(N11)*N12);
    disp('Calculating Jbar+...');
    Jbar_plus = simplify(pinv(Jbar));%Jbar'/(Jbar*Jbar');

    DJ = sym('DJ', size(J));
    [a,b] = size(J);
    for i =1:a
        for j = 1:b
            temp = 0;
            for k = 1:noLinks
               temp = temp + diff(J(i,j),q(j))*Dq(j) + diff(J(i,j),Dq(j))*D2q(j);
            end
            DJ(i,j) = temp;
        end
    end

    disp('Calculating control variables...');
    control_v = sym ('control_v', [2 1]);   %TODO - generalise size of control_v
    qad = sym ('qad', size(qa));                    % desired actuated joints
    qpd = sym ('qpd', size(qp));                    % desired passive joints
    control_tau = N21*qpd + N22*qad + D2*Dq + L2;   % controller torques
    matlabFunction(control_tau, 'File', 'control_tau','Vars', {qpd,qad,[q;Dq],mass,r,g});
    qpd = N11\(N12*qad + D1*Dq +L1);               % find qpd
    matlabFunction(qpd, 'File', 'D2q1_desired','Vars', {qad,[q;Dq],mass,r,g});
    qad = Jbar_plus*(control_v - DJ*Dq + (J1/N11)*(D1*Dq + L1)); % find qad
    matlabFunction(qad, 'File', 'D2q2_desired','Vars', {control_v,[q;Dq],mass,r,g});
end


%% Other
disp('Complete.');