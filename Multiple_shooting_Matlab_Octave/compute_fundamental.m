% written by Pakpoom Khamchuai

function U = compute_fundamental(a, b, d, y, Gradf_u)
% d = dimension of fundamental matrix
% need all y solved by IVP to evaluate GradF

h = 1/30;
U = diag(ones(d,1));
for i = 1:d
    U(:,i) = integrator_ee(U(:,i), a, b, h, y, Gradf_u);
end
end

function un = integrator_ee(u0, t0, tn, h, y, f)
% Explicit Euler method integrator for Fundamental Matrix
% y is values of u(t,u), needed for evaluation of Fundamental Matrix
% f is the Jacobian/Gradient
u = u0;
t = t0;
epsilon = 10e-7;
i = 1;
while true
    f_val = f(t, y(1,i))*u ;
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
