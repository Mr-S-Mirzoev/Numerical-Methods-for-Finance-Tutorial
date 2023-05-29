% Parameters
a = 0;
b = 1;
T = 1;
theta = 0.5; % Choose a value between 0 and 1

% Maximum values of N and M
N_max = 10;
M_max = 10;

% Arrays to store error values
errors_N = zeros(N_max, 1);
errors_M = zeros(M_max, 1);

% Loop over N values
for N = 1:N_max
    M = N^2; % Choose M = N^2 for convergence analysis
    
    % Compute error
    error = heateq_fdm(a, b, T, round(N), round(M), theta);
    
    % Store error value
    errors_N(N) = error;
end

% Loop over M values
for M = 1:M_max
    N = round(sqrt(M)); % Choose N = sqrt(M) for convergence analysis
    
    % Compute error
    error = heateq_fdm(a, b, T, round(N), round(M), theta);
    
    % Store error value
    errors_M(M) = error;
end

% Plot the convergence w.r.t N
figure;
plot(1:N_max, errors_N, '-o');
xlabel('N');
ylabel('Error');
title('Convergence of Error with respect to N');

% Plot the convergence w.r.t M
figure;
plot(1:M_max, errors_M, '-o');
xlabel('M');
ylabel('Error');
title('Convergence of Error with respect to M');
