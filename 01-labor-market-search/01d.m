clear all; close all; clc

% Initial values
c_0     = 0.5;
alpha_0 = 0.5;
beta_0  = 0.98;
lam_0   = 0.6;

% Parameter grid
cmin = 0.01; cmax = 0.99; c_i = 100;  
c_vals = linspace(cmin, cmax, c_i);

betamin = 0.01; betamax = 0.99; b_i = 100; 
beta_vals = linspace(betamin, betamax, b_i);

alphamin = 0.01; alphamax = 0.99; alpha_i = 100;
alpha_vals = linspace(alphamin, alphamax, alpha_i);

lammin = 0.01; lammax = 0.99; lambda_i = 100;
lambda_vals = linspace(lammin, lammax, lambda_i);

% Inicial value of w
w0 = 5;

%% Varying c
wb1 = zeros(1, c_i);
for j = 1:c_i
    G = @(w) w - c_vals(j) - ((alpha_0 * beta_0) / (1 - (beta_0 * (1 - lam_0)))) * (((log(10 / (w + 1))) + ((w +1 ) / 10)));
    wb1(j) = fsolve(G, w0, optimoptions('fsolve', 'Display', 'off'));
end

%% Varying beta
wb2 = zeros(1, b_i);
for j = 1:b_i
    G = @(w) w - c_0 - ((alpha_0 * beta_vals(j)) / (1 - (beta_vals(j) * (1 - lam_0)))) * (((log(10 / (w + 1))) + ((w +1 ) / 10)));
    wb2(j) = fsolve(G, w0, optimoptions('fsolve', 'Display', 'off'));
end

%% Varying alpha
wb3 = zeros(1, alpha_i);
for j = 1:alpha_i
    G = @(w) w - c_0 - ((alpha_vals(j) * beta_0) / (1 - (beta_0 * (1 - lam_0)))) * (((log(10 / (w + 1))) + ((w +1 ) / 10)));
    wb3(j) = fsolve(G, w0, optimoptions('fsolve', 'Display', 'off'));
end

%% Varying lambda
wb4 = zeros(1, lambda_i);
for j = 1:lambda_i
    G = @(w) w - c_0 - ((alpha_0 * beta_0) / (1 - (beta_0 * (1 - lambda_vals(j))))) * (((log(10 / (w + 1))) + ((w +1 ) / 10)));
    wb4(j) = fsolve(G, w0, optimoptions('fsolve', 'Display', 'off'));
end

%% Results

figure;
subplot(2, 2, 1);
plot(alpha_vals, wb3, 'LineWidth', 1.5);
title("W* Varying Alpha (all else being equal)");
xlabel("Alpha"); ylabel("W*");

subplot(2, 2, 2);
plot(beta_vals, wb2, 'LineWidth', 1.5);
title("W* Varying Beta (all else being equal)");
xlabel("Beta"); ylabel("W*");

subplot(2, 2, 3);
plot(lambda_vals, wb4, 'LineWidth', 1.5);
title("W* Varying Lambda (all else being equal)");
xlabel("Lambda"); ylabel("W*");

subplot(2, 2, 4);
plot(c_vals, wb1, 'LineWidth', 1.5);
title("W* Varying C (all else being equal)");
xlabel("C"); ylabel("W*");