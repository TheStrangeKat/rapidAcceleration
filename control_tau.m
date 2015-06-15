function control_tau = control_tau(qpd1,in2,in3,in4,in5,g)
%CONTROL_TAU
%    CONTROL_TAU = CONTROL_TAU(QPD1,IN2,IN3,IN4,IN5,G)

%    This function was generated by the Symbolic Math Toolbox version 5.10.
%    21-May-2015 15:59:28

Dthe1 = in3(4,:);
Dthe2 = in3(5,:);
Dthe3 = in3(6,:);
mass2 = in4(2,:);
mass3 = in4(3,:);
qad1 = in2(1,:);
qad2 = in2(2,:);
r1 = in5(1,:);
r2 = in5(2,:);
r3 = in5(3,:);
the1 = in3(1,:);
the2 = in3(2,:);
the3 = in3(3,:);
t2 = r2.^2;
t3 = cos(the2);
t4 = the2+the3;
t5 = sin(the2);
t6 = mass2.*t2.*(1.0./3.0);
t7 = mass3.*t2;
t8 = r3.^2;
t9 = mass3.*t8.*(1.0./3.0);
t10 = cos(the3);
t11 = mass3.*r2.*r3.*t10;
t12 = the1+the2;
t13 = cos(t12);
t14 = sin(the3);
t15 = sin(t4);
t16 = mass3.*r1.*r3.*t15.*(1.0./2.0);
t17 = r3.*2.0;
t18 = cos(t4);
t19 = r2.*t10.*3.0;
t20 = t17+t19;
t21 = the1+the2+the3;
t22 = cos(t21);
control_tau = [qad1.*(t6+t7+t9+t11)+Dthe1.*(Dthe1.*(t16+mass2.*r1.*r2.*t5.*(1.0./2.0)+mass3.*r1.*r2.*t5)-Dthe3.*mass3.*r2.*r3.*t14.*(1.0./2.0))+qpd1.*(t6+t7+t9+t11+mass2.*r1.*r2.*t3.*(1.0./2.0)+mass3.*r1.*r2.*t3+mass3.*r1.*r3.*t18.*(1.0./2.0))+g.*mass3.*(r2.*t13+r3.*t22.*(1.0./2.0))+g.*mass2.*r2.*t13.*(1.0./2.0)+mass3.*qad2.*r3.*t20.*(1.0./6.0)-Dthe3.*mass3.*r2.*r3.*t14.*(Dthe1+Dthe2+Dthe3).*(1.0./2.0)-Dthe2.*Dthe3.*mass3.*r2.*r3.*t14.*(1.0./2.0);Dthe1.*(Dthe1.*(t16+mass3.*r2.*r3.*t14.*(1.0./2.0))+Dthe2.*mass3.*r2.*r3.*t14.*(1.0./2.0))+mass3.*qad2.*t8.*(1.0./3.0)+g.*mass3.*r3.*t22.*(1.0./2.0)+mass3.*qad1.*r3.*t20.*(1.0./6.0)+mass3.*qpd1.*r3.*(t17+t19+r1.*t18.*3.0).*(1.0./6.0)+Dthe2.*mass3.*r2.*r3.*t14.*(Dthe1+Dthe2).*(1.0./2.0)];
