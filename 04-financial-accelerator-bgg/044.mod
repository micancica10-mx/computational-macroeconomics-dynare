// BGG MODEL

//Endogenous Variables 
var y c i g c_e r n r_k q k x a h pi r_n spread;

//Exogenous Variables 
varexo e_a e_g e_r e_n;

//Parameters
parameters beta delta alpha omega eta G_Y I_Y C_Y C_e_Y fi kappa gamma RK_N dseta rho rho_a rho_g epsilon v;

beta=0.80; 
delta=0.040; 
alpha=0.45;
omega=1; 
eta=3; 
G_Y=0.2; 
I_Y=0.18; 
C_Y=0.5; 
C_e_Y=0.12;
fi=0.25; 
kappa=0.0858; 
gamma= 0.87;
RK_N=2; 
dseta=0.21; 
rho=0.90; 
rho_a=0.90; 
rho_g=0.90; 
epsilon=0.9; 
v=0.8;


model;

//Equations

//1
y= C_Y*c + I_Y*i+G_Y*g+C_e_Y*c_e;

//2
c=-r+c(+1);

c_e=n;

r_k(+1)-r=-v*(n-(q+k));

r_k=(1-epsilon)*(y-k(-1)-x) + epsilon*q-q(-1);

q=fi*(i-k(-1));

y=a+alpha*k(-1)+(1-alpha)*omega*h;

y-h-x-c= h/eta;

pi=-kappa*x+beta*pi(+1);

k=delta*i+(1-delta)*k(-1);

n=gamma*RK_N*(r_k-r(-1))+r(-1)+n(-1) + e_n;

r_n= rho*r_n(-1)+ dseta* pi(-1)+ e_r;

r_n=r+pi(+1);

a=rho_a*a(-1)+e_a;

g=rho_g*g(-1)+e_g;

spread=r_k-r;

end;

//stady state
steady;

//Guarantee of equation closure
resid;

//Blanchard-khan Conditions
check;

shocks;
var e_a=1;
var e_n=1;
end;

stoch_simul(order=1, irf=20);


// Variables selected for plotting the IRFs
vars = {'y', 'c', 'i', 'r_k', 'r_n', 'r', 'q', 'spread', 'n'}; 

% Plot for shock e_r (shock to 'r')
figure;
num_vars = length(vars);  

for j = 1:num_vars
    subplot(3, 3, j);  % 3x3 grid

    % Check if the variable IRF data exists
    if isfield(oo_.irfs, [vars{j} '_e_n'])
        plot(1:options_.irf, eval(['oo_.irfs.' vars{j} '_e_n']), 'LineWidth', 1.5);
        title(vars{j});
        xlabel('Periods');
        ylabel('IRF');
        grid on;

        % Customize the y-axis for the variable 'i'
        if strcmp(vars{j}, 'i')
            ylim([-0.1, 0.1]);  
        end
    else
        
        title([vars{j} ' (No Data)']);
        axis off;
    end
end
sgtitle('IRF - net worth shock'); 


// Impulse-response graph for shocks at a
figure;

for j = 1:num_vars
    subplot(3, 3, j);  % 3x3 grid


    if isfield(oo_.irfs, [vars{j} '_e_a'])
        plot(1:options_.irf, eval(['oo_.irfs.' vars{j} '_e_a']), 'LineWidth', 1.5);
        title(vars{j});
        xlabel('Periods');
        ylabel('IRF');
        grid on;


        if strcmp(vars{j}, 'i')
            ylim([-0.1, 0.1]);  
        end
    else
        
        title([vars{j} ' (No Data)']);
        axis off;
    end
end
sgtitle('IRF - PTF shock');  

//
// Shocks simulation
stoch_simul(order=1, irf=20);
irf_v_1_8 = oo_.irfs;  % Store the IRFs for phi_pi=1.2

// Set v = 0 and run the model again
v=0;
stoch_simul(order=1, irf=20);
irf_v_1_0 = oo_.irfs;  % Store the IRFs for phi_pi=1.9
//
// Impulse-response graph for shocks at n

for j = 1:num_vars
    subplot(3, 3, j);  

    % Check if the IRF data exists for both cases
    if isfield(irf_v_1_8, [vars{j} '_e_n']) && isfield(irf_v_1_0, [vars{j} '_e_n'])
        plot(1:options_.irf, eval(['irf_v_1_8.' vars{j} '_e_n']), 'LineWidth', 1.5, 'DisplayName', 'v = 0.8');
        hold on;
        plot(1:options_.irf, eval(['irf_v_1_0.' vars{j} '_e_n']), 'LineWidth', 1.5, 'DisplayName', 'v = 0');
        title(vars{j});
        xlabel('Periods');
        ylabel('IRF');
        legend;
        grid on;
        
        % Customize y-axis for 'i' if needed
        if strcmp(vars{j}, 'i')
            ylim([-0.1, 0.1]);  
        end
    else
        title([vars{j} ' (No Data)']);
        axis off;
    end
end
sgtitle('IRF - net-worth shock'); 

