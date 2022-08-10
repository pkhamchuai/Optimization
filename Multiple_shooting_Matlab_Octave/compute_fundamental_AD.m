% written by Pakpoom Khamchuai
% AD codes need a fix.

function U = compute_fundamental_AD(a, b, d, y)
% f = A(t).u(t) = Gradf_u * u
% d = dimension of fundamental matrix
% need all y solved by IVP to evaluate fundamental matrix
% with AD, the difference is how to evaluate GradF_func
% which is f_val in integrator_ee

h = 1/30;
U = diag(ones(d,1));
for i = 1:d     % evaluate each column of fundamental matrix
    U(:,i) = integrator_ee_AD(U(:,i), a, b, h, d, y);
end
end

function un = integrator_ee_AD(u0, t0, tn, h, d, y)
% Explicit Euler method integrator for Fundamental Matrix
% y is values of u(t,u), needed for evaluation of Fundamental Matrix
% f is the Jacobian/Gradient

u = u0;
t = t0;
epsilon = 10e-7;
i = 1;
while true
    f_val = autodiff(t, y(:,i), d)*u;
    if abs((t+h) - tn) < epsilon % stop computation when we reach tn
        un = u + h*f_val;
        break;
    else
        if t+h > tn    % reduce step size if we are clse to tn
            h = tn - t;
            un = u + h*f_val;
            break;
        else            % normal step size h
            un = u + h*f_val;
            t = t + h;
        end
    end
    u = un;
    i = i+1;
end
end

function f_val = autodiff(tm, y, d)
% each loop computes each column of Jacobian

for i = 1:d         
    % independence v
    vm2 = y(1);
    vm1 = y(2);
    v0  = tm;
    
    % v_dot: directional derivative of v
    if i == 1
        vm2d = 1;
        vm1d = 0;
        v0d  = 0;
    else
        vm2d = 0;
        vm1d = 1;
        v0d  = 0;
    end
    
    % intermediates
    v1  = v0*vm1;
    v1d = (v0*vm1d) + (v0d*vm1);
    v2  = 4*vm2.^(3/2);
    v2d = 6*sqrt(vm2);
    
    % dependent
    y1  = v1;
    y1d = v1d;
    y2  = v2;
    y2d = v2d;
    
    f_val(:,i) = [y1d; y2d];
end
end
