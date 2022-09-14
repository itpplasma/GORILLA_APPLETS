close all

% Plotting Style
fsize = 28;
lwidth = 4;

% Green
color_green = [18/255,173/255,42/255];

% Path for RAW Data
path_data = '../RAW_Data/alpha_life_time_t_scan/t_1d0/';

x_70 = importdata([path_data,'alpha_life_time_gorilla_poly4_70.dat']);
x_100 = importdata([path_data,'alpha_life_time_gorilla_poly4_100.dat']);
x_200 = importdata([path_data,'alpha_life_time_gorilla_poly4_200.dat']);

x_70_rk45 = importdata([path_data,'alpha_life_time_gorilla_rk45_70.dat']);
x_100_rk45 = importdata([path_data,'alpha_life_time_gorilla_rk45_100.dat']);
x_200_rk45 = importdata([path_data,'alpha_life_time_gorilla_rk45_200.dat']);

x_odeint_1em8_part_1 = importdata([path_data,'alpha_life_time_gorilla_direct_1e-8_part_1.dat']);
x_odeint_1em8_part_2 = importdata([path_data,'alpha_life_time_gorilla_direct_1e-8_part_2.dat']);
x_odeint_1em8_part_3 = importdata([path_data,'alpha_life_time_gorilla_direct_1e-8_part_3.dat']);
x_odeint_1em8_part_4 = importdata([path_data,'alpha_life_time_gorilla_direct_1e-8_part_4.dat']);
x_odeint_1em8_part_5 = importdata([path_data,'alpha_life_time_gorilla_direct_1e-8_part_5.dat']);
x_odeint_1em8_part_6 = importdata([path_data,'alpha_life_time_gorilla_direct_1e-8_part_6.dat']);


t_confined_70 = x_70.data(:,6);
t_confined_100 = x_100.data(:,6);
t_confined_200 = x_200.data(:,6);

t_confined_70_rk45 = x_70_rk45.data(:,6);
t_confined_100_rk45 = x_100_rk45.data(:,6);
t_confined_200_rk45 = x_200_rk45.data(:,6);

t_confined_odeint_1em8 = [  x_odeint_1em8_part_1.data(:,6); ...
                            x_odeint_1em8_part_2.data(:,6);
                            x_odeint_1em8_part_3.data(:,6);
                            x_odeint_1em8_part_4.data(:,6);
                            x_odeint_1em8_part_5.data(:,6);
                            x_odeint_1em8_part_6.data(:,6) ];         

t_sample = logspace(-4,log10(max(t_confined_100)),100).';

confined_ratio_70 = zeros(size(t_sample));
confined_ratio_100 = zeros(size(t_sample));
confined_ratio_200 = zeros(size(t_sample));

confined_ratio_70_rk45 = zeros(size(t_sample));
confined_ratio_100_rk45 = zeros(size(t_sample));
confined_ratio_200_rk45 = zeros(size(t_sample));

confined_ratio_odeint_1em8 = zeros(size(t_sample));

n_sampling = numel(t_confined_100);

% Throw away extra samples
t_confined_odeint_1em8 = t_confined_odeint_1em8(1:n_sampling,1);

for i = 1:numel(t_sample)
   
    L_70 = t_confined_70 >= t_sample(i);
    L_100 = t_confined_100 >= t_sample(i);
    L_200 = t_confined_200 >= t_sample(i);
    
    L_70_rk45 = t_confined_70_rk45 >= t_sample(i);
    L_100_rk45 = t_confined_100_rk45 >= t_sample(i);
    L_200_rk45 = t_confined_200_rk45 >= t_sample(i);
    
    L_odeint_1em8 = t_confined_odeint_1em8 >= t_sample(i);
    
    confined_ratio_70(i) = sum(L_70)/n_sampling;
    confined_ratio_100(i) = sum(L_100)/n_sampling;
    confined_ratio_200(i) = sum(L_200)/n_sampling;
    
    confined_ratio_70_rk45(i) = sum(L_70_rk45)/n_sampling;
    confined_ratio_100_rk45(i) = sum(L_100_rk45)/n_sampling;
    confined_ratio_200_rk45(i) = sum(L_200_rk45)/n_sampling;
    
    confined_ratio_odeint_1em8(i) = sum(L_odeint_1em8)/n_sampling;
    
end


% standard deviation
confined_sigma_70 = sqrt(confined_ratio_70.*(1-confined_ratio_70) / n_sampling);
confined_sigma_100 = sqrt(confined_ratio_100.*(1-confined_ratio_100) / n_sampling);
confined_sigma_200 = sqrt(confined_ratio_200.*(1-confined_ratio_200) / n_sampling);

confined_sigma_70_rk45 = sqrt(confined_ratio_70_rk45.*(1-confined_ratio_70_rk45) / n_sampling);
confined_sigma_100_rk45 = sqrt(confined_ratio_100_rk45.*(1-confined_ratio_100_rk45) / n_sampling);
confined_sigma_200_rk45 = sqrt(confined_ratio_200_rk45.*(1-confined_ratio_200_rk45) / n_sampling);

confined_sigma_odeint_1em8 = sqrt(confined_ratio_odeint_1em8.*(1-confined_ratio_odeint_1em8) / n_sampling);







%% Plot 1: Confined ratio vs. time
fh = figure;
fh.Units = 'normalized';
fh.Position = [0.2995 0 0.5292 0.9037];%[0.3542 0.0787 0.4833 0.8250];

hold on
% Shading

%N_grid = 70
% f0 = fill([t_sample;flipud(t_sample)],[confined_ratio_70-1.96.*confined_sigma_70;flipud(confined_ratio_70+1.96.*confined_sigma_70)],'r','linestyle','none');
% 
% %N_grid = 100
% f1 = fill([t_sample;flipud(t_sample)],[confined_ratio_100-1.96.*confined_sigma_100;flipud(confined_ratio_100+1.96.*confined_sigma_100)],color_green,'linestyle','none');
% 
% %N_grid = 200
% f2 = fill([t_sample;flipud(t_sample)],[confined_ratio_200-1.96.*confined_sigma_200;flipud(confined_ratio_200+1.96.*confined_sigma_200)],'b','linestyle','none');
% 
% 
% %N_grid = 70, RK45
% f3 = fill([t_sample;flipud(t_sample)],[confined_ratio_70_rk45-1.96.*confined_sigma_70_rk45;flipud(confined_ratio_70_rk45+1.96.*confined_sigma_70_rk45)],'m','linestyle','none');
% 
% %N_grid = 100, RK45
% f4 = fill([t_sample;flipud(t_sample)],[confined_ratio_100_rk45-1.96.*confined_sigma_100_rk45;flipud(confined_ratio_100_rk45+1.96.*confined_sigma_100_rk45)],'c','linestyle','none');
% 
% %N_grid = 200, RK45
% f5 = fill([t_sample;flipud(t_sample)],[confined_ratio_200_rk45-1.96.*confined_sigma_200_rk45;flipud(confined_ratio_200_rk45+1.96.*confined_sigma_200_rk45)],'y','linestyle','none');
% 
% 
% %ODEINT, eps = 1E-8
 f6 = fill([t_sample;flipud(t_sample)],[confined_ratio_odeint_1em8-1.96.*confined_sigma_odeint_1em8;flipud(confined_ratio_odeint_1em8+1.96.*confined_sigma_odeint_1em8)],'k','linestyle','none');


 set([f6],'FaceAlpha',0.2);

h0 = plot(t_sample,confined_ratio_70);
h1 = plot(t_sample,confined_ratio_100);
h2 = plot(t_sample,confined_ratio_200);

h3 = plot(t_sample,confined_ratio_70_rk45);
h4 = plot(t_sample,confined_ratio_100_rk45);
h5 = plot(t_sample,confined_ratio_200_rk45);

h6 = plot(t_sample,confined_ratio_odeint_1em8);

hold off

set([h0,h1,h2,h3,h4,h5,h6],'LineWidth',lwidth)


h0.Color = 'r';
h0.LineStyle = '-';
h3.Color = 'r';
h3.LineStyle = ':';

h1.Color = color_green;
h1.LineStyle = '-';
h4.Color = color_green;
h4.LineStyle = ':';

h2.Color = 'b';
h2.LineStyle = '-';
h5.Color = 'b';
h5.LineStyle = ':';

h6.Color = 'k';
h6.LineStyle = '-.';





xlabel('$t_\mathrm{trace}$ [s]','FontSize',fsize+2,'Interpreter','latex');
ylabel('$f_c$','FontSize',fsize+2,'Interpreter','latex');

a1 = gca;
a1.FontSize = fsize;
a1.XScale = 'log';
grid on
box on

l = legend([h6,h0,h3,h1,h4,h2,h5],{ ...
    'Exact Orbits', ...
    'GORILLA Poly4, $N_{s}=N_{\vartheta}=N_{\varphi} = 70$', ...
    'GORILLA RK4/5, $N_{s}=N_{\vartheta}=N_{\varphi} = 70$', ...
    'GORILLA Poly4, $N_{s}=N_{\vartheta}=N_{\varphi} = 100$', ...
    'GORILLA RK4/5, $N_{s}=N_{\vartheta}=N_{\varphi} = 100$', ...
    'GORILLA Poly4, $N_{s}=N_{\vartheta}=N_{\varphi} = 200$', ...
    'GORILLA RK4/5, $N_{s}=N_{\vartheta}=N_{\varphi} = 200$' ...
    });
l.Interpreter = 'latex';
l.Location = 'southwest';
l.FontSize = fsize+3;

t1 = title('$\mathbf{a}$: $t_\mathrm{trace} = 1.0$ s, $N_\mathrm{particles} = 10^3$');
t1.Interpreter ='latex';
t1.FontSize = fsize+2;

ylim([0.75,1])

