% get linear approximation on triangle
% PARAM_IN:
%   x, y  - [3, 1] arrays, coordinates of triangle's vertices,
%   z     - [3, 1] array, function values in vertices
% PARAM_OUT:
%   z_f   - anonymous function of two arguments, linear approximation,
%   dz_dx - DBL, gradient value, x-direction
%   dz_dy - DBL, gradient value, y-direction
function [z_f, dz_dx, dz_dy] = getLinearApprox(x, y, z)

    A = det([y(2)-y(1),z(2)-z(1);y(3)-y(1),z(3)-z(1)]);
    B = det([x(2)-x(1),z(2)-z(1);x(3)-x(1),z(3)-z(1)]);
    C = det([x(2)-x(1),y(2)-y(1);x(3)-x(1),y(3)-y(1)]);
        
    z_f = @(xx, yy) (-A*(xx - x(1)) + B*(yy - y(1)))/C + z(1);
    dz_dx = -A/C;
    dz_dy =  B/C;
    
end
