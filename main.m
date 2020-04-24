%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% HOMEWORK #5
% Joshua Julian Damanik (20194701)
% AE551 - Introduction to Optimal Control
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear, clc, close all;
addpath('lib');

X0 = [1000, 0]';
N = 5;


eps = 10;
k = 15;

marker_list = {'o-', '+-', '^-', 'v-',  'x-', '.-'};
legends = {};

plot_data = zeros(N+1,3,1);

rho = 0.1;
delta_rho = 2;

iter = 0;


X = zeros(3*N,1);
X(1:N) = X0(1);
X(N+1:2*N) = X0(2);

[A, B, R] = cost_function_param(N, X0);
c = @(X) (A*X-B);
lambda = zeros(length(B),1);

for j = 1:5
    iter = iter + 1;
    
    La = @(X, lambda) (0.5*X'*R*X + lambda'*(A*X-B)+rho.*((A*X-B))'*(A*X-B));

    last_X_star = X;
    
    for ii = 1:10000
        
        
        f = @(X) La(X,lambda);
        X_star = minimize(X, eps, k, f);

        quasi = quasi_newton_class(length(X));
        for n = 1:1000
            p = quasi.bgfs(X, eps, f);
            [Xa, Xb] = unimodal_interval(X, eps, p, f);
            X_star = fibonacci_search(Xa, Xb, k, f);

            err = norm(X_star - X);
            X = X_star;
            if (err < eps * 0.01)
                break;
            end
        end
        
%         X_star = -(2*R+2*rho*(A'*A)) \ (lambda'*A-2*(B'*A))';
%         X = X_star;

        L_star = lambda + 2*rho*c(X_star);
        err = norm(L_star - lambda);
        fprintf('j=%d,ii=%d,err=%.4f\n',j,ii,err);
        if (err < 2*rho*0.01)
            break;
        end
        X = X_star;
        lambda = L_star;
    end

    plot_data(1,1,j) = X0(1);
    plot_data(1,2,j) = X0(2);
    plot_data(2:N+1,1,j) = X_star(1:N);
    plot_data(2:N+1,2,j) = X_star(N+1:2*N);
    
    plot_data(:,3,j) = [X_star(2*N+1:3*N); 0];
    
    legends = [legends, {sprintf('rho=%d', rho)}];
    
    if (A*X-B)'*(A*X-B) <= 0.01
        break;
    end
    
    rho = rho*delta_rho;
    
end

title = {'Vertical Position', 'Vertical Speed', 'Vertical Input'};
ylabels = {'Y (m)', 'Vy (m/s)', 'U (m/s2)'};
legend_position = {'NorthEast', 'SouthEast', 'SouthEast'};
for i = 1:3
    figure('Name', title{i}), hold on;
    set(gcf, 'Position', [100*i, 100*i, 1000, 400]);
    for j = 1:iter
        figure(i),plot(5000 / N * (0:N), plot_data(:,i,j), marker_list{mod(j-1,length(marker_list))+1});
    end
    grid on;
    xlabel('x (m)');
    ylabel(ylabels{i});
    legend(legends, 'Location', 'Best');
end
