function x = joint(I1, I2)
x = zeros(256);
[r c] = size(I1);
    for u=1:r
       for v=1:c
            x(I1(u,v)+1, I2(u,v)+1) = x(I1(u,v)+1, I2(u,v)+1) +1 ; 
       end
    end
end