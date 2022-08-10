% written by Pakpoom Khamchuai

function [deltav, ua] = MultipleShooting(v, t0, tn, ns)

h = 1/30;                   % step size for Explicit Euler
g = [1; 0];                 % boundary condition of original problem
Ba = [1 0; 0 0];
Bb = [0 0; 1 0];
tm = linspace(t0,tn,ns+1);  % time at shooting nodes/boundaries

u = zeros(2,ns+1);
GradF = zeros(2*(ns+1),2*(ns+1));

% analytical Jacobian (not used in AutoDiff)
GradF_func = @(t, u) [0 t; 6*u(1)^(1/2) 0];

% compute u(i) and GradF(v)
for i = 1:ns                % for every shooting nodes
    y = expl_euler(v(:,i), tm(i), tm(i+1), h);
    u(:,i) = y(:,end);      % store values of last y for uk(tk)
    
    if i == 1
        ua = y(:,1);        % get ua from 1st subinterval, needed for r
    end
    if i == ns
        ub = y(:,end);      % get ub from last subinterval, needed for r
    end
    
    % to compute GradF: put G, -I, Ba, Bb in their positions
    % G is computed by integration of U' over interval of shooting nodes
    % done in compute_Variational func
    
    %%%%%%%%%%% switch between analytical and with AD %%%%%%%%%%%%%%%
    GradF(2*i-1:2*i,2*i-1:2*i) = compute_fundamental(tm(i), tm(i+1), 2, y, GradF_func);
    %GradF(2*i-1:2*i,2*i-1:2*i) = compute_fundamental_AD(tm(i), tm(i+1), 2, y);
    
    GradF(2*i-1:2*i,2*i+1:2*i+2) = -eye(2);
    GradF(2*ns+1:2*ns+2,1:2) = Ba;
    GradF(2*ns+1:2*ns+2, 2*ns+1:2*ns+2) = Bb;
    
end

r = Ba*ua + Bb*ub - g;          % boundary condition function r
F = compute_F(u, v, r, ns);     % compute F

% disp(['u1(0) = ' num2str(ua(1))])
% disp(['u2(0) = ' num2str(ua(2))])

% solve linear system GradF*v = -F for v (update of shooting nodes)
deltav = linsolve(GradF,-F);

% --------------------
% Globalization with line search (with Armijo & Backtracking)
% reduce Newton step size with t, t = (0,1]
% t which minimizes F(u, v+t*deltav, r, ns)
% return deltav = t*deltav
deltav = backtracking(deltav, u, v, r, ns);
% --------------------

end

function F = compute_F(u, v, r, ns)
% to computes F, need u, v, r and ns

F = zeros(2,ns+1);
for i = 1:ns            % compute F(v)
    F(:,i) = u(:,i) - v(:,i+1);
end
F(:,end) = r;           % last 2 elements of F(v)
F = reshape(F,numel(F),1);
end

function deltav = backtracking(deltav, u, v, r, ns)
% compute a new deltav with backtrackng
% find t = (0,1] which minimizes F(u, v+t*deltav, r, ns)
% first check if t = 1 satisfies Armijo (mu = 0)
% if not, reduce t and check again
% returns deltav = t*deltav
% needs u, v, r, ns for compute_F func

t = 1;

for i = 1:18     % step from t = 1 back to some t that is not so small
    v = reshape(v,numel(v),1);
    v_new = v + t*deltav;
    v_new = reshape(v_new,2,ns+1);
    v = reshape(v,2,ns+1);
    
    % compute norm of F(x+td) and F(x)
    normF = norm(compute_F(u, v_new, r, ns));
    normF_old = norm(compute_F(u, v, r, ns));
    
    % if F is "better", return deltav, if not, reduce t and repeat
    % if there is no t satiefies the condition, use t = 1
    if normF < normF_old
        %disp(['t = ' num2str(t)])
        deltav = t*deltav;
        break
    else t = t - 0.05;
    end
end
end
