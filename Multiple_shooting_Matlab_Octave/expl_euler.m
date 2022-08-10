% written by Pakpoom Khamchuai

function y = expl_euler(u0, t0, tn, h)
u = u0;
t = t0;
epsilon = 1e-7;
y = u0;

while true
    f_val = rhs_evaluate(t, u);   % evaluate values of rhs at t
    
    if abs((t+h) - tn) < epsilon    % stop computation when we reach tn
        un = u + h*f_val;           % un = new u
        y = [y,un];
        break;
    else
        if t+h > tn    % reduce step size if we are clse to tn
            h = tn - t;
            un = u + h*f_val;
            y = [y,un];
            break;
        else            % normal step size h
            un = u + h*f_val;
            t = t + h;
        end
    end
    y = [y,un];
    u = un;
end
end

function f_val = rhs_evaluate(t, u)
f = @(t, u) [t*u(2); 4*(u(1)^(3/2))];
f_val = f(t, u);
end
