function ddq = ddq_func(in1,mass1,r01,g,k_spring,in6)
%DDQ_FUNC
%    DDQ = DDQ_FUNC(IN1,MASS1,R01,G,K_SPRING,IN6)

%    This function was generated by the Symbolic Math Toolbox version 6.2.
%    05-Jun-2015 17:21:04

Dr1 = in1(4,:);
Dthe1 = in1(3,:);
r1 = in1(2,:);
the1 = in1(1,:);
ddq = [-(g.*cos(the1)+Dr1.*Dthe1.*2.0)./r1;(k_spring.*r01-k_spring.*r1-g.*mass1.*sin(the1)+Dthe1.^2.*mass1.*r1)./mass1];
