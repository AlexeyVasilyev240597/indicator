% integrate function of 2 arguments by triangle
% PARAM_IN:
%   f    - anonymous function of two arguments,
%   x, y - [3, 1] arrays, coordinates of triangle's vertices
%   mode - STR, method of integration:
%          'int2' - using embaded integral2 method,
%          'q2'   - second order quadrature (by default),
%          'q2s1' - 2nd order quadrature with 1 times splitting
% PARAM_OUT:
%   I   - integral of f(x, y) by triangle 
function I = intByTriangle(f, x, y, mode)
    arguments
        f (1, :)
        x (1, :) {mustBeNumeric,mustBeReal}
        y (1, :) {mustBeNumeric,mustBeReal}
        mode (1,:) char {mustBeMember(mode,{'int2','q2','q2s1'})} = 'q2'
    end

    if strcmp(mode, 'int2')
        [x, i] = sort(x);
        y = y(i);    

        l_1 = @(xx) (y(3)-y(2))/(x(3)-x(2))*(xx - x(2))+y(2);
        l_2 = @(xx) (y(1)-y(3))/(x(1)-x(3))*(xx - x(3))+y(3);
        l_3 = @(xx) (y(2)-y(1))/(x(2)-x(1))*(xx - x(1))+y(1);

        I = 0;
        if x(1) ~= x(2)
            if y(1) < y(2)
                I = I + integral2(f, x(1), x(2), l_2, l_3);
            else
                I = I + integral2(f, x(1), x(2), l_3, l_2);
            end
        end
        if x(2) ~= x(3)
            if y(2) < y(3)
                I = I + integral2(f, x(2), x(3), l_1, l_2);
            else
                I = I + integral2(f, x(2), x(3), l_2, l_1);
            end
        end
    elseif strcmp(mode, 'q2')
        J = det([x;y;ones(1, 3)]);
        A = [1, 1, 0; 0, 1, 1; 1, 0, 1]/2;
        I = J*sum(f(x*A, y*A))/6;
    elseif strcmp(mode, 'q2s1') 
        % x*A = xc - x-coordinates of middles of triangle edges
        A = [1, 1, 0; 1, 0, 1; 0, 1, 1]/2;
        X = [x, x*A]; Y = [y, y*A];         
        I = 0;
        I = I + intByTriangle(f, X([1,4,5]), Y([1,4,5]));
        I = I + intByTriangle(f, X([2,6,4]), Y([2,6,4]));
        I = I + intByTriangle(f, X([3,5,6]), Y([3,5,6]));
        I = I + intByTriangle(f, X([4,6,5]), Y([4,6,5]));        
    else
        error('Unknown integration method')
    end
end

