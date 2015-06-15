function energy = getEnergy(in1,mass1,r01,g,k_spring)
%GETENERGY
%    ENERGY = GETENERGY(IN1,MASS1,R01,G,K_SPRING)

%    This function was generated by the Symbolic Math Toolbox version 6.2.
%    05-Jun-2015 17:21:03

Dr1 = in1(4,:);
Dthe1 = in1(3,:);
r1 = in1(2,:);
the1 = in1(1,:);
t2 = r01-r1;
energy = k_spring.*t2.^2.*(1.0./2.0)+mass1.*(Dr1.^2+Dthe1.^2.*r1.^2).*(1.0./2.0)+g.*mass1.*r1.*sin(the1);