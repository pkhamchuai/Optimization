% written by Pakpoom Khamchuai

% The definition of u'(t,u) is in rhs_evaluate func in expl_euler.m
% The computation of F and GradF is done in MultipleShooting func
% MultipleShooting func returns 'deltav' with globalization
% Globalization is implemented in MultipleShooting func
% can switch GradF between analytical and AutoDiff in MultipleShooting func

clear all;

t0 = 0.0;
tn = 5.0;

for node_nr = 10:20 % compare iterations needed for each node nr.
    
    ns = node_nr;                    % number of shooting nodes
    tm = linspace(t0,tn,ns+1);  % time at shooting nodes/boundaries
    u0 = [1.0; 0.0];            % initial values of u(0)
    v = zeros(2,ns+1);          % initial values of shooting nodes vk
    v(:,1) = u0;
    
    % Newton method iteration
    iteration = 0;
    while true
        %disp(['Iteration ' num2str(iteration)])
        [deltav, ua] = MultipleShooting(v, t0, tn, ns);
        v = reshape(v,numel(v),1);  % reshape v corresponding to deltav's size
        v = v + deltav;
        v = reshape(v,2,ns+1);
        
        % iterate until deltav is small enough
        if norm(deltav) < 1e-10
            break
        end
        
        plot(tm,real(v(1,:)),'-',tm,real(v(2,:)),'-')
        legend('u_1(t)','u_2(t)','Location','southeast')
        xlabel('Time t')
        ylabel('u(t)')
        title(['Multiple Shooting: ' num2str(ns) ' nodes'])
        ylim([-4 1.5])
        xlim([0 5])
        drawnow
        
        iteration = iteration + 1;
    end
    
    % display some summary
    disp(['nr. nodes: ' num2str(ns) ' | nr. iterations: ' num2str(iteration) ' | u2(0) = ' num2str(ua(2))])
end
