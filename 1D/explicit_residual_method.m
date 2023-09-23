% ******************************
% *********** LAB #2 ***********
% ***Explicit Residual Method***
% ******************************

% problem:
%    ordinary differential equation
%       (alpha(x)*u(x)')' + f(x) = 0, x in [0, 1],
%    boundary conditions
%       u(0) = u0, u(1) = u1.
%
% resuidual r(v) := (alpha(x)*v(x)')' + f(x),
%
% u_h(x) is finite dimensional approximation of u(x) (e.g. Galerkin one),
%
% task:
%   compare
%       alpha_0*||u - u_h|| <> h/pi*||r(u_h)||, 
%       where ||.|| is L^2 norm.


% notation:
%   suffix 'd'   - derivative,
%   suffix 'ad'  - antiderivative,
%   suffix 'sym' - symbolic variable or symbolic expression,
%   suffix 'f'   - function dependent on real variable,
%   suffix 'p'   - polynomial function,
%   suffix 'val' - vector of real values,

close all
clear all

%% setting of research params
% number of intervals for discretisation
N = 128;

flag_Gal_approx = false;
% fraction of noise
frac_val = 0.015;


%% init params of problem
syms x_sym
alpha_sym   = (152*x_sym^3 - 234*x_sym^2 + 97*x_sym + 24)/24;
alpha_d_sym = diff(alpha_sym);
alpha_ad_sym = int(alpha_sym);

u_sym     = sin(8*pi*x_sym);
u_d_sym   = diff(u_sym);

f_sym     = -diff(alpha_sym*u_d_sym);
f_ad_sym  = -alpha_sym*u_d_sym;
f_mul_x_ad_sym = int(x_sym*f_sym);

%% functions for conversions
% convert symbolic function into real value function
sym2f   = @(f_sym) @(x) double(subs(f_sym, x_sym, x));
% claclulate symbolic function on real variable
sym2val = @(f_sym, x_val) double(subs(f_sym, x_sym, x_val));

%% discretisation of problem
h = 1/N;
x = 0:h:1;
alpha_0 = min(sym2val(alpha_sym, x));

u0 = sym2val(u_sym, x(1));
u1 = sym2val(u_sym, x(N+1));

% finite elements method with piecewise linear approximation
alpha_ad_val = sym2val(alpha_ad_sym, x);
D       =  1/h^2*(alpha_ad_val(3:N+1) - alpha_ad_val(1:N-1));
subD    = -1/h^2*(alpha_ad_val(3:N)   - alpha_ad_val(2:N-1));
superD  = -1/h^2*(alpha_ad_val(3:N)   - alpha_ad_val(2:N-1));  
A = diag(D) + diag(subD, -1) + diag(superD, 1);

F_f = @(a, b, x1, x2) a.*(sym2val(f_mul_x_ad_sym, x2) -...
                         sym2val(f_mul_x_ad_sym, x1)) +...
                      b.*(sym2val(f_ad_sym, x2) - sym2val(f_ad_sym, x1));
F = zeros(N-1, 1);
F(1)   = 1/h^2*(alpha_ad_val(2)  -alpha_ad_val(1))*u0;
F(N-1) = 1/h^2*(alpha_ad_val(N+1)-alpha_ad_val(N))*u1;

i_in = 2:N;
F = F + 1/h*(F_f(1, -x(i_in-1), x(i_in-1), x(i_in)) +...
        F_f(-1, x(i_in+1), x(i_in), x(i_in+1)))';

u_h_val = [u0; A\F; u1];

u_val = sym2val(u_sym, x)';
fprintf('relative error in calculated nodes = %.3e \n',...
        norm(u_val - u_h_val)/norm(u_val))
fprintf('error in calculated nodes = %.3e \n \n', norm(u_val - u_h_val))

if ~flag_Gal_approx
    % noise
    eps = norm(u_h_val, Inf)*frac_val*(2*rand(N+1, 1)-1);
    
    u_h_val = u_h_val + eps;    
end

u_h_p = mkpp(x, [(u_h_val(2:end) - u_h_val(1:end-1))/h, u_h_val(1:end-1)]); 
u_h_d_p  = fnder(u_h_p, 1);
u_h_dd_p = fnder(u_h_p, 2);

%% plot approximation results
NNN = N*10;
hhh = 1/NNN;
xxx = 0:hhh:1;

% line width
lw = 2;

figure
plot(x, u_val, 'o', 'LineWidth', lw)
hold on
plot(xxx, sym2val(u_sym, xxx), 'LineWidth', lw)
hold on
plot(xxx, ppval(u_h_p, xxx), 'r--', 'LineWidth', lw)
grid on
legend('u_i', 'u', 'u_h')

figure
plot(x, sym2val(u_d_sym, x), 'o', 'LineWidth', lw)
hold on
plot(xxx, sym2val(u_d_sym, xxx), 'LineWidth', lw)
hold on
plot(xxx, ppval(u_h_d_p, xxx), 'r--', 'LineWidth', lw)
grid on
legend('(du/dx)_i', 'du/dx', 'du_h/dx')

%% calculate parts of inequality
%
% error
err = zeros(N, 1);
% residual
res = zeros(N, 1);

for n = 1:N
    
    hh = h/1e3;
    xx = x(n):hh:x(n+1);
    
    er = (sym2val(u_d_sym, xx) - ppval(u_h_d_p, xx)).^2;
    err(n) = trapz(xx, er);
    
    % r(u_h) = f + alpha(x)'*u_h(x)' + alpha(x)*u_h(x)''
    rr = (sym2val(f_sym, xx) +...
          sym2val(alpha_d_sym, xx).*ppval(u_h_d_p, xx) +...
          sym2val(alpha_sym, xx).*ppval(u_h_dd_p, xx)).^2;
    res(n) = trapz(xx, rr);
    
    if ~rem(n, round(N/8))
            fprintf('%.1f %% is done \n', n/N*100)
    end
    
end
fprintf('\n\n')

norm_const_e = alpha_0^2;
norm_const_r = (h/pi)^2;

fprintf('error = %.3e \n\n', norm_const_e*sum(err))
fprintf('residual = %.3e \n\n', norm_const_r*sum(res))

% plot histograms
fig = figure('position', [200, 200, 600, 480]);

bar((x(1:N)+x(2:N+1))/2', norm_const_r*res, 'r');  
hold on
bar((x(1:N)+x(2:N+1))/2', norm_const_e*err, 0.4);      

grid on
xlim([x(1), x(end)])
legend('residual', 'error')
