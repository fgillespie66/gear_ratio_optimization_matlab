function I = rodInertia(mass, radius, length)
I = zeros(3,3);
I(1,1) = mass*radius^2/4 + mass*length^2/12;
I(2,2) = mass*radius^2/4 + mass*length^2/12;
I(3,3) = mass*radius^2/2;
end