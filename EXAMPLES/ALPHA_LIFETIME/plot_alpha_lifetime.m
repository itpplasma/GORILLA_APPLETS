close all

% Plotting Style
fsize = 28;
lwidth = 4;

% Green
color_green = [18/255,173/255,42/255];

t_confined = load('alpha_life_time_gorilla.dat');
t_sample = logspace(-4,log10(max(t_confined)),100).';

confined_ratio = zeros(size(t_sample));
n_sampling = numel(t_confined);


for i = 1:numel(t_sample)
    L = t_confined >= t_sample(i);
    confined_ratio(i) = sum(L)/n_sampling;
end


% standard deviation
confined_sigma = sqrt(confined_ratio.*(1-confined_ratio) / n_sampling);



%% Plot 1: Confined ratio vs. time
fh = figure;
fh.Units = 'normalized';
fh.Position = [0.2995 0 0.5292 0.9037];%[0.3542 0.0787 0.4833 0.8250];

hold on
% Shading

 f6 = fill([t_sample;flipud(t_sample)],[confined_ratio-1.96.*confined_sigma;flipud(confined_ratio+1.96.*confined_sigma)],'k','linestyle','none');


 set([f6],'FaceAlpha',0.2);

h0 = plot(t_sample,confined_ratio);

hold off

set([h0],'LineWidth',lwidth)


h0.Color = 'r';
h0.LineStyle = '-';


xlabel('$t_\mathrm{trace}$ [s]','FontSize',fsize+2,'Interpreter','latex');
ylabel('$f_c$','FontSize',fsize+2,'Interpreter','latex');

a1 = gca;
a1.FontSize = fsize;
a1.XScale = 'log';
grid on
box on

l = legend(h0,'GORILLA Poly4, $N_{s}=N_{\vartheta}=N_{\varphi} = 101$');
l.Interpreter = 'latex';
l.Location = 'southwest';
l.FontSize = fsize+3;

t1 = title('$\mathbf{a}$: $t_\mathrm{trace} = 0.01$ s, $N_\mathrm{particles} = 10^4$');
t1.Interpreter ='latex';
t1.FontSize = fsize+2;

ylim([0.9,1])

