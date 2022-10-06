close all

% Load Poincaré invariant over time steps
poincare_invariant_data = load(['./results/poincare_invariance_gorilla_n_1E7_e_var_jperp_var_dtdtau_ham.dat']);



f1 = figure;
f1.Units = 'normalized';
f1.Position = [0.0797 0.1657 0.6271 0.6565];  
    
%
p_1 = plot(poincare_invariant_data(:,1),poincare_invariant_data(:,2),'-');
p_1.DisplayName = '$a_\perp = 0.1$, $a_H = 0.1$: Correct (Hamiltonian) gyro-phase $\phi$ \& time dynamics $\mathrm{d} t / \mathrm{d} \tau$';

%
l_1 = legend;
set([p_1],'LineWidth',3);
set(gca,'FontSize',18);
l_1.FontSize = 22;
l_1.Location = 'northwest';
l_1.Interpreter = 'latex';
%
ylabel('$J$','interp','latex')
xlabel('$\tau_\mathrm{step}$','interp','latex')
%
box on
xlim([0,200]);





f2 = figure;
f2.Units = 'normalized';
f2.Position = [0.0797 0.1657 0.6271 0.6565];   
%
p_2 = semilogy(normalize_tau(poincare_invariant_data(:,1)),normalize_pi(poincare_invariant_data));
p_2.DisplayName = p_1.DisplayName;
p_2.Color = p_1.Color;
p_2.Marker = p_1.Marker;
p_2.LineStyle = p_1.LineStyle;

%
l_2 = legend;
set([p_2],'LineWidth',3);
set(gca,'FontSize',18);
l_2.Location = 'southwest';
l_2.FontSize = 22;
l_2.Interpreter = 'latex';
%
ylabel('$\Delta J ~/~ J (\tau=0)$','interp','Latex')
xlabel('$\tau~/~\tilde{\tau}$','interp','latex')
%
box on
grid on
xlim([0,1]);
ylim([1E-20,1E-2]);









function [inv_out] = normalize_pi(inv_in)
    inv_out(:) = abs((inv_in(:,2)-inv_in(1,2))/inv_in(1,2));
end

function [tau_out] = normalize_tau(tau_in)
    tau_out(:) = tau_in(:)/tau_in(end);
end
