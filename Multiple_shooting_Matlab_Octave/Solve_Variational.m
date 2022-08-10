% f = A(t).u(t) = Gradf_u * u
% d = dimension of grad_f matrix
function U = Solve_Variational(a, b, d, y, Gradf_u)
    h = 1/30;        
    
    % U is made as a diagonal matrix to extract initial values of vectors
    % U is used to store final fundamental matrix also
    U = diag(ones(d,1));
    for i = 1:d
        U(:,i) = integrator_ee(U(:,i), a, b, h, y, Gradf_u); 
    end
    U;
end

% Explicit Euler method integrator for Fundamental Matrix
% f is the Jacobian/Gradient
function un = integrator_ee(u0, t0, tn, h, y, f)
    um = u0;
    tm = t0;
    epsilon = 10e-7;
    i = 1;
    while true
        f_val = f(tm, y(1,i))*um ;   
        if abs((tm+h) - tn) < epsilon % we are sufficiently close
            un = um + h*f_val;        
            break;
        else
            if tm+h > tn
                % This should not hit if n is a correct multiple of h
                h = tn - tm;
                un = um + h*f_val;
                break;
            else
                un = um + h*f_val;
                tm = tm + h;
            end
        end
        um = un;
        i = i+1;
    end
end