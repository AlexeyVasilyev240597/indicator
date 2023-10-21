% construct projection of the gradient of the numerical solution as a grid function
% (based on minimizing M_+ functional being 
% a majorant for an A-posteriori estimate of the error)
% PARAM_OUT:
%   p_h - [N_p, 2], projection of the gradient of the numerical solution defined in mesh nodes
function p_h = MPlus(obj)

    ar = pdetrg(obj.p, obj.t)';
    N_p = length(obj.p); 
    N_t = length(obj.t);
    
    % calc init value of p_h = grad(u) as average gradient
    p_h = obj.averageGrad();
    p_h = p_h(:);
    %  minimization of M_+
    max_iter = 1;
    options = optimoptions('fminunc', 'MaxIter', max_iter);
    [p_h, ~] = fminunc(@Majorant, p_h, options);

    p_h = reshape(p_h, N_p, 2);    

    % M_+ functional
    % PARAM_IN:    
    %   p_h    - projection of the gradient
    % PARAM_OUT:
    %   M_plus - DBL, value of the functional
    function M_plus = Majorant(p_h)   
        
        err = 0;
        rsd = obj.f_norm_L2_sq;
        
        p_h = reshape(p_h, N_p, 2);
                
        z_p = cell(2, 1);
        for i = 1:N_t
            x = obj.p(1, obj.t(1:3, i));
            y = obj.p(2, obj.t(1:3, i));
    
            z1 = p_h(obj.t(1:3, i), 1);
            [z_p{1}, dz1_dx, ~] = Indicator.getLinearApprox(x, y, z1);
            z2 = p_h(obj.t(1:3, i), 2);
            [z_p{2}, ~, dz2_dy] = Indicator.getLinearApprox(x, y, z2);
    
            trm1_f = @(xx, yy) (obj.du_h_dx(i) - z_p{1}(xx, yy)).^2 +...
                               (obj.du_h_dy(i) - z_p{2}(xx, yy)).^2;
            
            err = err + Indicator.intByTriangle(trm1_f, x, y);
            c_i = dz1_dx + dz2_dy;
            rsd = rsd + 2*c_i*Indicator.intByTriangle(obj.f, x, y) + c_i^2*ar(i);        
        end    
        
        rsd = obj.C_F_sq*rsd;        
        % optimal beta
        beta = sqrt(rsd/err);        
        M_plus = (1+beta)*err + (1+1/beta)*rsd;  
        
    end

end
