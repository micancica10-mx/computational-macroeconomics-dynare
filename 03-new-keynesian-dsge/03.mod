// Endogenous Variables
var a c n y pi pi_star x1 x2 v_p w mc i r m v_pf n_f y_f x dm;

//Exogenous variables
varexo e_a e_m; 

//Parameters
parameters sigma beta psi eta theta epsilon phi rho_m rho_a pi_bar;
sigma=2;
beta=0.99;
psi=1;
eta=1;
theta=2;
epsilon=5;
phi=0.6;
rho_a=0.95;
rho_m=0;
pi_bar=0;

// Definition of model equations
model;
// 1. Consumption Euler Equation
    c^(-sigma) = beta * (c(+1)^(-sigma)) * (1 + i) / (1 + pi(+1));
    
    // 2. Labor Supply Equation
    psi * n^eta = c^(-sigma) * w;
    
    // 3. Money Demand Equation
    m = theta * ((1 + i) / i) * c^(sigma);
    
    // 4. Marginal Cost Equation
    w = a * mc; 
    
    // 5. Aggregate Output Equation
    c = y;
    
    // 6. Aggregate Output Production Function
    y = (a * n) / v_p;
    
    // 7. Price Dispersion Equation
    v_p = (1 - phi) * (1 + pi_star)^(-epsilon) * (1 + pi)^epsilon + (1 + pi)^epsilon * phi * v_p(-1);
    
    // 8. Recursive Form for Price Dispersion
    (1 + pi)^(1 - epsilon) = (1 - phi) * (1 + pi_star)^(1 - epsilon) + phi;
    
    // 9. Optimal Pricing Equation
    1 + pi_star = (epsilon/(epsilon - 1))* (1 + pi) * (x1/x2);
    
    // 10. Marginal Cost Components - x1
    x1 = c^(-sigma) * mc * y + phi * beta * (1 + pi(+1))^epsilon * x1(+1);
    
    // 11. Marginal Cost Components - x2
    x2 = c^(-sigma) * y + phi * beta * (1 + pi(+1))^(epsilon - 1) * x2(+1);
    
    // 12. Productivity Shock Process
    log(a) = rho_a * log(a(-1)) + e_a;

    // 14. Definition of Money Growth
    dm = log(m) - log(m(-1));
    
    // 13. Money Growth Equation
    dm = (1 - rho_m) * pi_bar - pi + rho_m * dm(-1) + rho_m * pi(-1) + e_m;
    
    // 15. Price Dispersion in Steady State
    v_pf = 1;
    
    // 16. Labor Supply in Steady State
    n_f = ((1 / psi) * ((epsilon - 1) / epsilon) * a^(1-sigma))^(1 / (sigma + eta));
    
    // 17. Output in the Steady State
    y_f = ((1 / psi) * ((epsilon - 1) / epsilon) * a^(1+eta))^(1/(sigma + eta));
    
    // 18. Output Gap Definition
    x   = log(y) - log(y_f);

    // 19. Fisher
    (1+r) = (1+i)/(1+pi(+1));  


end;

//Initial values
initval;
r=(1/beta)-1;
pi=pi_bar;
a=1;
mc = (( ((1 + pi_bar)^(1 - epsilon) - phi)/((1 - phi)*(1 + pi_bar)^(1-epsilon)))^(1/(1 - epsilon))) * ( (1 - phi * beta * (1 + pi_bar)^epsilon) / (1 - phi * beta * (1 + pi_bar)^(epsilon - 1))) * ((epsilon - 1) / epsilon);
pi_star=(( ((1 + pi_bar)^(1 - epsilon) - phi)/((1 - phi)))^(1/(1 - epsilon))) - 1;
w=mc;
v_p=((1-phi) * (1+pi_star)^(-epsilon) * (1 + pi_bar)^(epsilon))/(1-phi*(1+pi_bar)^(epsilon)); 
n=((w*v_p^(sigma))/psi)^(1/(eta+sigma)); 
i=((1+pi_bar)/beta)-1;
y=n/v_p; 
c=y; 
m=theta * ((1+i)/i)*c^(sigma);
x1=(c^(-sigma) * mc * y)/(1-phi*beta*(1+pi_bar)^(epsilon));
x2=(c^(-sigma) * y)/(1-phi*beta*(1+pi_bar)^(epsilon-1)); 
v_pf=1; 
n_f=((1/psi) * ((epsilon - 1)/epsilon))^(1/(eta+sigma)); 
y_f=((1/psi) * ((epsilon - 1)/epsilon))^(1/(eta+sigma));  
x=log(y)-log(y_f); 
dm=0;
end;

resid;
//SS
steady;

// stability
check;

//Shock
shocks;
var e_a = 0.01;
var e_m = 0.01;
end;

// Dinamic Process
stoch_simul(periods=1000, order = 1, irf = 20);

// Variables selected for plotting the IRFs
vars = {'a', 'y', 'pi', 'i', 'r', 'mc', 'n', 'm', 'pi_star', 'y_f', 'x'}; 

% Plot for shock e_a (shock to 'a')
figure;
num_vars = length(vars);  % number of variables

for j = 1:num_vars
    subplot(3, 4, j);  % 3x4 grid

    % Check if the variable IRF data exists
    if isfield(oo_.irfs, [vars{j} '_e_a'])
        plot(1:options_.irf, eval(['oo_.irfs.' vars{j} '_e_a']), 'LineWidth', 1.5);
        title(vars{j});
        xlabel('Periods');
        ylabel('IRF');
        grid on;

        % Customize the y-axis for the variable 'i'
        if strcmp(vars{j}, 'i')
            ylim([-0.1, 0.1]);  % Adjust y-axis range as needed
        end
    else
        % Display a message if data is missing
        title([vars{j} ' (No Data)']);
        axis off;
    end
end
sgtitle('IRF - PTF shock');



// impulse-response graph for shocks at i
figure;

for j = 1:num_vars
    subplot(3, 4, j);  % 3x4 grid

    % Check if the variable IRF data exists
    if isfield(oo_.irfs, [vars{j} '_e_m'])
        plot(1:options_.irf, eval(['oo_.irfs.' vars{j} '_e_m']), 'LineWidth', 1.5);
        title(vars{j});
        xlabel('Periods');
        ylabel('IRF');
        grid on;

        % Customize the y-axis for the variable 'i'
        if strcmp(vars{j}, 'i')
            ylim([-0.1, 0.1]);  % Adjust y-axis range as needed
        end
    else
        % Display a message if data is missing
        title([vars{j} ' (No Data)']);
        axis off;
    end
end
sgtitle('IRF - money supply shock');  
