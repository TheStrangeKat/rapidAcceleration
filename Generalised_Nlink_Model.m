%% Generalised_Nlink_Model
%
%   Calculates the Equations of Motion for a N-Link Pendulum
%
%   NF Steenkamp
%   University of Cape Town
%   Electrical Engineering Department
%   11 May 2015 - Last Updated 20 May 2015
%% 
Passive = false;
Actuated = true;

%% Settings
noLinks = 3;
joints = [Passive;Actuated;Actuated];

%% Error Checking
if(length(joints)~=noLinks)
    % See settings above
    error('The Actuated Joints vector must match the number of links...');
end

% if(COMlink > noLinks)
%     error('The link containing the COM must be one of the links of the system.');
% end

%% Declare Variables-------------------------------------------------------
disp('...');
disp('Declaring Variables and calculating positions/velocities...');

% General variables
q = sym('the', [noLinks 1]);    % Generalized Coordinates
sym(q,'real');
Dq = sym('Dthe', [noLinks 1]);  % First Derivative of q
D2q = sym('D2the', [noLinks 1]);% Second Derivative of q
r = sym('r', [noLinks 1]);      % Length of links
sym(r,'real');
x = r/2;                        % COM of links (relative to frame)
%COM                            % COM of system
tau = sym('tau', [noLinks 1]);  % Torques (control inputs)
mass = sym('mass', [noLinks 1]);% Link Masses
sym(mass,'real');
%I = sym('I',[noLinks 1]);       % Mass Moment of Inertia for Links
I = mass.*(r.^2)/12;            % Mass Moment of Inertia for Links
R0 = cell(noLinks,1);           % Rotation Matrices
T = sym('T', [noLinks 1]);      % Kinectic Energy
U = sym('U', [noLinks 1]);      % Potential Energy
EL = sym('EL', [noLinks 1]);    % Euler-Lagrange
syms g Ttot Utot L;

disp('...Manipulator Equations...');
% Manipulator Equations M*D2q + C*Dq + G = Bu
M = sym('M', noLinks); %Inertial Matrix
sym(M,'real');
C = sym('C', noLinks); % Coriolis and Centripetal Terms
sym(C,'real');
G = sym('G', [noLinks 1]); % Gravitational Terms
sym(G,'real');
B = sym('B', noLinks); % Torque Terms
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
%position mass relative to frame origin
pm = sym('pm', [3,noLinks]);
pm(1,:)=x.';
pm(2:3,:)=0;
for i = 1:noLinks
    pm(:,i) = R0{i}*pm(:,i);
end

%position of frame origin relative to inertial frame
po = sym('po', [3,noLinks]);
po(1,:)=r.';
po(2:3,:)=0;
for i = 1:noLinks
    po(:,i) = R0{i}*po(:,i);
end

%position of mass relative to inertial frame
ym = sym('ym', [3,noLinks]);
ym(1,:)=r.';
ym(2:3,:)=0;
for i = 1:noLinks
    j = 1;
    ym(:,i) = pm(:,i);
    while (j<i)
        ym(:,i) = ym(:,i) + po(:,j);
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
for i = 1:noLinks
    j = 1;
    v(:,i) = cross(omega(:,i), pm(:,i));
    while (j<i)
        v(:,i) = v(:,i) + cross(omega(:,j), po(:,j));
        j = j+1;
    end
end

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
    U(i) = mass(i)*g*U(i);
    Utot = Utot + U(i);
end
Utot = simplify(Utot);

%% Lagrangian - (currently not used)
% disp('Calculating Lagrangian Equations...');
% L = Ttot - Utot;
% sym temp;
% for i = 1:noLinks
%     EL(i) = diff(L,Dq(i));    % Euler-Lagrange
%     temp = 0;
%     for j = 1:noLinks
%        temp = temp + diff(EL(i),q(j))*Dq(j) + diff(EL(i),Dq(j))*D2q(j);
%     end
%     EL(i) = temp;
%     EL(i) = EL(i) - diff(L,q(i)) - tau(i);
%     EL(i) = simplify(EL(i));
% end

%% Manipulator Equations---------------------------------------------------
disp('Calculating Manipulator Equations...');
% Inertial Matrix
for i = 1:noLinks
    for j = 1:noLinks
        M(i,j) = diff(diff(Ttot,Dq(i)),Dq(j));
    end
end
M = simplify(M);

% Coriolis Matrix
C(:,:) = 0;
for i = 1:noLinks
    for j = 1:noLinks
        for k = 1:noLinks
                C(i,j) = C(i,j)+ 0.5*(diff(M(i,j),q(k)) + diff(M(i,k),q(j)) - diff(M(j,k),q(i)))*Dq(k);
        end        
    end
end
C = simplify(C);

% Gravitational Matrix
for i = 1:noLinks
    G(i) = diff(Utot,q(i));
end
G = simplify(G);

% Torque Matrix
B = eye(noLinks) - diag(~joints); % Torque Matrix

% Equations of motion
ddq = simplify(-M\(C*Dq + G - B*tau));      % Solving for D2q
tau = sym('tau', [noLinks 1]);              % Allowing for construction of function
matlabFunction(ddq, 'File', 'ddq','Vars', {[q;Dq],mass,r,g,tau});   %Creating function
tau(joints==Passive) = sym(0);              % Re-correcting for passive joints

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

%% Other
disp('Complete.');