function [p, u] = cfdGetPlanesIntersection(n1,a1,n2,a2)

p = [0 0 0];
u = cross(n1,n2);

%  test if the two planes are parallel
if norm(u) < 10^-7                % Plane 1 and Plane 2 are near parallel
    V = a1 - a2;
    if (dot(n1,V)==0)
        check = 1;                    % Plane 1 and Plane 2 coincide
        return
    else
        check = 0;                   %Plane 1 and Plane 2 are disjoint
        return
    end
end

check = 2;

% Plane 1 and Plane 2 intersect in a line
% first determine max abs coordinate of cross product
maxc = find(abs(u)==max(abs(u)));


%next, to get a point on the intersection line and
%zero the max coord, and solve for the other two

d1 = -dot(n1,a1);   %the constants in the Plane 1 equations
d2 = -dot(n2, a2);  %the constants in the Plane 2 equations

switch maxc
    case 1                   % intersect with x=0
        p(1) = 0;
        p(2) = (d2*n1(3) - d1*n2(3))/ u(1);
        p(3) = (d1*n2(2) - d2*n1(2))/ u(1);
    case 2                    %intersect with y=0
        p(1) = (d1*n2(3) - d2*n1(3))/ u(2);
        p(2) = 0;
        p(3) = (d2*n1(1) - d1*n2(1))/ u(2);
    case 3                    %intersect with z=0
        p(1) = (d2*n1(2) - d1*n2(2))/ u(3);
        p(2) = (d1*n2(1) - d2*n1(1))/ u(3);
        p(3) = 0;
end


