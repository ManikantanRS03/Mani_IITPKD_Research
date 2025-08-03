% Initialization
alpha = -1e4;
beta = 1e4;

% Function definitions
F = @(P, E) (alpha*P.^2)/2 + (beta*P.^4)/4 - beta*P.*E;

% Plot Double Well with Switching Dynamics
P = linspace(-2, 2, 100);
E = 0.2;
figure;
plot(P,F(P,E));
xlabel('P');
ylabel('F(P, E)');
title('Double Well with Switching Dynamics');

% Solve for PA, PB, PC using fsolve
solFP = @(P, E) alpha*P + beta*P.^3 - beta*E;

function [PA, PB, PC] = findPAPC(E, solFP)
    PA = fsolve(@(P) solFP(P, E), -2);
    PB = fsolve(@(P) solFP(P, E), 2);
    PC = fsolve(@(P) solFP(P, E), 0);
end

% Function to find Kramers time
function rk = findrk(Vpluse, D, alpha, beta)
    E = Vpluse;
    [PA, ~, PC] = findPAPC(E, @(P, E) alpha*P + beta*P.^3 - beta*E);
    F2derPA = alpha + 3*beta*PA^2;
    F2derPC = alpha + 3*beta*PC^2;
    deltaF = abs((alpha*PC^2)/2 + (beta*PC^4)/4 - beta*PC*E - ((alpha*PA^2)/2 + (beta*PA^4)/4 - beta*PA*E));
    
    % Kramers time calculation
    rk = sqrt(abs(F2derPA * F2derPC)) / (2*pi) * exp(-deltaF / D);
end

% Function to find dU
function [du,PA,PC] = finddU(Vpluse, alpha, beta)
    E = Vpluse;
    [PA, ~, PC] = findPAPC(E, @(P, E) alpha*P + beta*P.^3 - beta*E);
    du = abs((alpha*PC^2)/2 + (beta*PC^4)/4 - beta*PC*E - ((alpha*PA^2)/2 + (beta*PA^4)/4 - beta*PA*E));
end

% Example usage
Vps = 0.25;
D = 1000;
rk = findrk(Vps, D, alpha, beta);

% Time and Rate calculation for different voltages
Vp = 0.3;
tp = logspace(-7, -1, 100);
D = 10000;

figure;
semilogx(tp, 1 - exp(-tp * findrk(Vp, D, alpha, beta)));
hold on;
semilogx(tp, 1 - exp(-tp * findrk(0.2, D, alpha, beta)));
semilogx(tp, 1 - exp(-tp * findrk(0.1, D, alpha, beta)));
xlabel('Time');
ylabel('Switching Probability');
title('Switching Dynamics over Time');

% Bias voltage vs. switching probability at two noise levels
Vps = linspace(0.1, 0.3, 100);
D1 = 10000;
D2 = 20000;
tps = 5e-4;

P1 = zeros(1, length(Vps));
P2 = zeros(1, length(Vps));

for i = 1:length(Vps)
    rk1 = findrk(Vps(i), D1, alpha, beta);
    P1(i) = 1 - exp(-tps * rk1);
    rk2 = findrk(Vps(i), D2, alpha, beta);
    P2(i) = 1 - exp(-tps * rk2);
end

figure;
plot(Vps, P1, Vps, P2);
xlabel('Bias Voltage');
ylabel('P(V_s)');
title('Bias Voltage vs Switching Probability');
legend('D1 = 10000', 'D2 = 20000');

% Bias voltage vs. Delta U
DU = zeros(1, length(Vps));
PA = zeros(1, length(Vps));
PC = zeros(1, length(Vps));

for i = 1:length(Vps)
    [DU(i),PA(i),PC(i)] = finddU(Vps(i), alpha, beta);
end

figure;
plot(Vps, DU);
xlabel('Bias Voltage');
ylabel('\Delta U');
title('Bias Voltage vs Potential Difference (\Delta U)');

figure;
plot(Vps, PA);
xlabel('Bias Voltage');
ylabel('PA');
title('min point');

figure;
plot(Vps, PC);
xlabel('Bias Voltage');
ylabel('PC');
title('max point');

%% Simulink output ploting (set P track, fixed noise)
d = out;
time = out.tout;
Pset = squeeze(out.p_set.signals.values);
Ptrk = squeeze(out.P_track.signals.values);
Dext = squeeze(out.Dext.signals.values);

figure(1);
clf;  % clear previous content

% Set figure size (in inches): [left bottom width height]
set(gcf, 'Units', 'inches', 'Position', [1, 1, 15, 5]);  % 6x4 inch figure

% Plot curves
plot(time, Pset, 'LineWidth', 1.5, 'DisplayName', 'P_{set}');
hold on;
plot(time, Ptrk, 'LineWidth', 1.5, 'DisplayName', 'P_{trk}');
ylim([0,1]);

% Labels and legend
xlabel('Time (s)', 'FontSize', 12);
ylabel('Probability (P)', 'FontSize', 12);
title('Probability Setpoint vs Tracked', 'FontSize', 14);
legend('Location', 'best');
grid on;

% Improve aesthetics
set(gca, 'FontSize', 11);
box on;

% Save as PDF in current folder
print(gcf, 'Probability_track', '-dpdf');

%% Simulink output ploting (set P track, variable noise)

d = out;
time1 = out.tout;
Ptrk1 = squeeze(out.P_track1.signals.values);
Pset1 = 0.5*ones(length(Ptrk1),1);
Dext = squeeze(out.Dext.signals.values);

figure(1);
clf;  % clear previous content

% Set figure size (in inches): [left bottom width height]
set(gcf, 'Units', 'inches', 'Position', [1, 1, 15, 5]);  % 6x4 inch figure

% Plot curves
plot(time1, Pset1, 'LineWidth', 1.5, 'DisplayName', 'P_{set}');
hold on;
plot(time1, Ptrk1, 'LineWidth', 1.5, 'DisplayName', 'P_{trk}');
ylim([0,1]);
hold off;

% Labels and legend
xlabel('Time (s)', 'FontSize', 12);
ylabel('Probability (P)', 'FontSize', 12);
title('Probability Setpoint vs Tracked', 'FontSize', 14);
legend('Location', 'best');
grid on;

% Improve aesthetics
set(gca, 'FontSize', 11);
box on;

% Save as PDF in current folder
print(gcf, 'Probability_track1', '-dpdf');


figure(2);
clf;  % clear previous content

% Set figure size (in inches): [left bottom width height]
set(gcf, 'Units', 'inches', 'Position', [1, 1, 15, 5]);  % 6x4 inch figure
plot(Dext);
% Improve aesthetics
set(gca, 'FontSize', 11);
box on;

% Save as PDF in current folder
print(gcf, 'dext_variation', '-dpdf');