function y = wrap_angle(u)
%#codegen   

    y = mod(u,2*pi);
end