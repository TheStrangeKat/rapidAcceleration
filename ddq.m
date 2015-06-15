function ddq_out = ddq(in1,mass1,r01,g,k_spring,in6)
%DDQ
%    DDQ = DDQ(IN1,MASS1,R01,G,K_SPRING,IN6)

%    This function was generated by the Symbolic Math Toolbox version 5.10.
%    22-May-2015 14:44:14

    Dr1 = in1(4,:);
    Dthe1 = in1(3,:);
    r1 = in1(2,:);
    the1 = in1(1,:);
    ddq_out = [-(g.*cos(the1)+Dr1.*Dthe1.*2.0)./r1;-(k_spring.*r1-k_spring.*r01+g.*mass1.*sin(the1)-Dthe1.^2.*mass1.*r1)./mass1];
end