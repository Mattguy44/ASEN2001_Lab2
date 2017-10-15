
xy = randi(60,1);
xz = randi(60,1);
pick = randi(84,9);
beam = zeros(81,3);
for i = 1:9
    for j = 1:9
        for L = 1:81
            bar = input(pick(j,i))
            cos(xz)*cos(xy) = x;
            cos(xz)*sin(xy) = y;
            sin(xz) = z;
            bx = bar*x;
            by = bar*y;
            bz = bar*z;
            beam(L,1) = bx;
            beam(L,2) = by;
            beam(L,3) = bz;
        end
        







    end
end