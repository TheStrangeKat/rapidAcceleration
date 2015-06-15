% Exhaustive Search

r0 = 2;
divA = 10;
divL = 10;
% angle = 1.1:(pi/2-1.1)/(divA-1):pi/2;
% length = r0/8:(r0 - r0/8)/(divL-1):r0;
angle = 1.1:(1.4-1.1)/(divA-1):1.4;
length = r0/8:(1.2 - r0/8)/(divL-1):1.2;

outSurf = zeros(divA,divL);

for i = 1:divA
    for j = 1:divL
        outSurf(i,j) = exhaustive([angle(i);length(j)]);
    end
end

surf(length,angle,outSurf);
xlabel('length');
ylabel('angles');
zlabel('horizontal velocity');

[A,Irow] = max(outSurf);
[val,Icol] = max(A);
Irow = Irow(Icol);
disp(outSurf(Irow,Icol));
%exhaustive([angle(Irow);length(Icol)]);