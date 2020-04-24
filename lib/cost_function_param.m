function [A, B, R] = cost_function_param(N, X0)
    dt = 20 / N;
    
    A11 = eye(N) - diag(ones(N-1,1),-1);
    A12 = -dt*diag(ones(N-1,1),-1);
    A13 = zeros(N);
    A1 = [A11, A12, A13];
    
    A21 = A13;
    A22 = A11;
    A23 = -dt*eye(N);
    A2 = [A21, A22, A23];
    
    A3 = [zeros(1, N-1), 1, zeros(1, 2*N)];
    A4 = [zeros(1, 2*N-1), 1, zeros(1, N)];
    
    A = [A1; A2; A3];%; A4];
    
    B1 = [X0(1)-X0(2)*dt; zeros(N-1,1)];
    B2 = [X0(2); zeros(N-1,1)];
    B3 = 0;
    B4 = 0;
    
    B = [B1; B2; B3];%; B4];
    
    
    R = blkdiag(zeros(2*N),eye(N));
end