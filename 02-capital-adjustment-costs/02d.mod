% Endogenous variables
var c a k n i y r w q;

%Definition of exogenous variables:
varexo e_a;

%Definition and calibration of model parameters:
parameters beta, psi, phi, delta, alpha, tau, rho; 

% Parameter values
beta = 0.93;
psi = 2;
phi = 0.5;
delta = 0.03;
alpha = 0.45;
tau = 0.5;
rho = 0.95;


model;
c=w/(psi*(n^(phi)));
r=alpha*(y/k(-1))-delta;
w=(1-alpha)*a*(k(-1)/n)^(alpha);
q=1/(1-tau*((i/k(-1))-delta));
q=beta*(c/c(+1))*(alpha*a(+1)*((k)^(alpha-1))*((n(+1))^(1-alpha))+q(+1)*(1-delta)-q(+1)*(tau/2)*(((i(+1)/k)-delta)^2)+q(+1)*tau*((i(+1)/k)-delta)*(i(+1)/k));
log(a)=rho*log(a(-1))+e_a;
y=c+i+(tau/2)*(((i/k(-1))-delta)^2)*k(-1);
y=a*((k(-1))^(alpha))*(n^(1-alpha));
k=i+(1-delta)*k(-1)-k(-1)*(tau/2)*(((i/k(-1))-delta)^2);
end;

%These values ​​were calculated in Excel

initval;


a=1;
w=(1-alpha)*a*(((1-alpha)*((alpha*beta/(1-beta*(1-delta)))^(alpha/(1-alpha))))^alpha);
n=(((1/psi)*((a*((alpha*beta*a)/(1-beta*(1-delta)))^(alpha/(1-alpha)))-(delta*((alpha*beta*a)/(1-beta*(1-delta)))^(1/(1-alpha))))^(-1))*w)^(1/(phi+1));
k=(((alpha*beta*a)/(1-beta*(1-delta)))^(1/(1-alpha)))*n;
y=a*(k^(alpha))*(n^(1-alpha));
i=delta*k;
c=y-i;
q=1;
r=(1/beta)-1;
end;

steady;
check;

shocks;
var e_a; stderr 0.01;
end;

%simulations depending on the number of periods. Sample moments are obtained.
stoch_simul(periods=1000, irf=50);
  
// Variables selected for plotting the IRFs
vars = {'a', 'c', 'k', 'n', 'i', 'y', 'r', 'w', 'q'}; 

% Plot for shock e_a (shock to 'a')
figure;
num_vars = length(vars); 

for j = 1:num_vars
    subplot(3, 3, j);  % 3x3 grid

    % Check if the variable IRF data exists
    if isfield(oo_.irfs, [vars{j} '_e_a'])
        plot(1:options_.irf, eval(['oo_.irfs.' vars{j} '_e_a']), 'LineWidth', 1.5);
        title(vars{j});
        xlabel('Periods');
        ylabel('IRF');
        grid on;

        % Customize the y-axis for the variable 'i'
        if strcmp(vars{j}, 'i')
            ylim([-0.1, 0.1]); 
        end
    else
        % Display a message if data is missing
        title([vars{j} ' (No Data)']);
        axis off;
    end
end
sgtitle('IRF -  PTF shock'); 
