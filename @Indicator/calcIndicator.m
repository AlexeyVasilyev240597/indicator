% construct the indicator as a field on the mesh
% PARAM_IN:
%   p_h  - [N_p, 2], a projection of gradient of the numerical solution
% PARAM_OUT:
%   indr - [N_t, 1], a field of indicator defined on elements of the mesh
function indr = calcIndicator(obj, p_h)

    N_t = length(obj.t);
    indr = zeros(N_t, 1);
    z_p = cell(2, 1);
    for i = 1:N_t
        x = obj.p(1, obj.t(1:3, i));
        y = obj.p(2, obj.t(1:3, i));

        for k = 1:2
            z = p_h(obj.t(1:3, i), k);
            z_p{k} = Indicator.getLinearApprox(x, y, z);
        end

        f_f = @(xx, yy) (obj.du_h_dx(i) - z_p{1}(xx, yy)).^2 +...
                        (obj.du_h_dy(i) - z_p{2}(xx, yy)).^2;
        indr(i) = Indicator.intByTriangle(f_f, x, y);
    end   
    
end
