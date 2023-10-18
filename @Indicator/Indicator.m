
classdef Indicator < handle
    properties (Access = private)
        % mesh edges
        e = [];
        % function, PDE coefficient: -nabla*(c nabla u) + a u = f
        f = [];
        % mesh points
        p = [];
        % mesh triangles
        t = [];

        % geometry description matrix
        gd = [];

        % numerical solution
        u_h = [];
        % Friedrichs' constant for a rectangle in power of 2
        C_F_sq = [];
        % (||f||_L2)^2
        f_norm_L2_sq = [];
    end

    methods (Access = public)
        % PARAM_IN:
        %   gd -geometry description matrix
        %   a  - PDE coefficient: -nabla*(c nabla u) + a u = f
        %   b  - boundary conditions
        %   c  - PDE coefficient: -nabla*(c nabla u) + a u = f
        %   e  - mesh edges
        %   f  - PDE coefficient: -nabla*(c nabla u) + a u = f
        %   p  - mesh points
        %   t  - mesh triangles
        function obj = Indicator(gd, a, b, c, e, f, p, t)
            obj.e = e;
            obj.p = p;            
            obj.t = t;
            
            % a rectangular domain
            if gd(1) == 3 && gd(2) == 4
                obj.gd = gd;
            else
                error('Domain type is not a rectangle.')
            end

            try
                obj.u_h = assempde(b, p, e, t, c, a, f);
            catch exception
                error(['Wrong set of init arguments. Error while solving the system: ' exception])
            end
            obj.f = str2func(['@(x, y)' f]);            
        end

        % calculate the indicator as a field on the mesh
        % PARAM_IN:
        %   projection_type - STR, types of projections:
        %       'AG' - average gradient,
        %       'MP' - minimizing the M_+ functional
        % PARAM_OUT:
        %   indr - [N_t, 1], a field of indicator defined on elements of the mesh
        function indr = getIndicator(obj, projection_type)
            arguments
                obj
                projection_type (1,:) char {mustBeMember(projection_type,{'AG','MP'})}
            end
            if strcmp(projection_type, 'AG')
                p_h = obj.averageGrad();
            elseif strcmp(projection_type, 'MP')
                x = obj.gd(3:6);  xmin = min(x); xmax = max(x);
                y = obj.gd(7:10); ymin = min(y); ymax = max(y);
                l1 = xmax - xmin;
                l2 = ymax - ymin;
                
                obj.C_F_sq = 1/(pi^2*(1/l1^2 + 1/l2^2));
                f_sq = @(x, y) (obj.f(x, y)).^2;

                obj.f_norm_L2_sq = integral2(f_sq, xmin, xmax, ymin, ymax);

                p_h = obj.MPlus();
            else
                error('Unknown projection type')
            end
            indr = obj.calcIndicator(p_h);
        end

        % calculate the error in the natural energy norm
        % PARAM_IN:
        %   du_dx, du_dy - functions, components of gradient of exact solution
        % PARAM_OUT:
        %   err - field of the error
        function err = calcError(obj, du_dx, du_dy)
            % approximate gradient
            [du_h_dx, du_h_dy] = pdegrad(obj.p, obj.t, obj.u_h);
            N_t = length(obj.t);
            err = zeros(N_t, 1);
            for i = 1:N_t
                x = obj.p(1, obj.t(1:3, i));
                y = obj.p(2, obj.t(1:3, i));
                err_f = @(xx, yy) (du_dx(xx, yy) - du_h_dx(i)).^2 +...
                                  (du_dy(xx, yy) - du_h_dy(i)).^2;
                err(i) = Indicator.intByTriangle(err_f, x, y);
            end
        end
        
        % plot field on the mesh
        % PARAM_IN:
        %   fld       - [N_t, 1], some field on the mesh,
        %   title_str - CHAR, field name
        function plotFld(obj, fld, title_str)
            arguments
                obj
                fld (1, :) {mustBeNumeric}
                title_str (1,:) char = ''
            end
            figure
            pdeplot(obj.p, obj.e, obj.t, 'xydata', fld,...
            'mesh', 'on', 'xystyle', 'flat',...
            'colormap', 'cool', 'colorbar', 'off');
            title(title_str)
        end
    end

    methods (Access = private)
        p_h = averageGrad(obj)
        p_h = MPlus(obj)
        indr = calcIndicator(obj, p_h)
    end

    methods (Static)        
        [z_f, dz_dx, dz_dy] = getLinearApprox(x, y, z)        
        I = intByTriangle(f, x, y, mode)
        fld_bool = marker(fld, type, kappa)
    end
    
end
