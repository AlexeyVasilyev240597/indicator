% construct a projection of the gradient of the numerical solution as a grid function
% (based on averaging the gradient)
% PARAM_OUT:
%   p_h - [N_p, 2], projection of gradient of the numerical solution defined in mesh nodes
function p_h = averageGrad(obj)
    
    N_p = length(obj.p);
    
    % grad_u_h = grad(u_h) on elements
    [du_h_dx, du_h_dy] = pdegrad(obj.p, obj.t, obj.u_h);
    grad_u_h = [du_h_dx; du_h_dy]';

    % p_h = grad(u_h) in nodes
    p_h = zeros(N_p, 2);

    % areas of all triangles
    ar = pdetrg(obj.p, obj.t)';
    for i = 1:N_p
        [~, patch_i] = find(obj.t(1:3, :) == i);
        T_j = ar(patch_i);
        w_i = sum(T_j);

        p_h(i, :) = sum(grad_u_h(patch_i, :).*(T_j*ones(1, 2)))/w_i;
    end
    
end