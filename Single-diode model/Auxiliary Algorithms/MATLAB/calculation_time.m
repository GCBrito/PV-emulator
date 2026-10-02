clear; clc; close all;
% Evaluates the calculation time at different output operating points.
% All data are embedded directly in this script.

%% 1) Data

useCapturedData = true;
I_raw = [
    3.02;
    2.99;
    2.96;
    2.88;
    2.86;
    2.68;
    2.58;
    2.21;
    2.02;
    1.76;
    1.50;
    1.29;
    1.11;
    0.98;
    0.73;
    0.52;
    0.38;
    0.27;
];
V_raw = [
    14.78;
    16.27;
    16.65;
    17.20;
    17.55;
    18.38;
    18.75;
    19.62;
    20.03;
    20.33;
    20.50;
    20.66;
    20.86;
    21.03;
    21.31;
    21.49;
    21.74;
    21.69;
];
I_capt = [
    3.000;
    2.955;
    2.933;
    2.888;
    2.840;
    2.715;
    2.613;
    2.290;
    2.130;
    1.910;
    1.727;
    1.620;
    1.456;
    1.340;
    1.140;
    0.995;
    0.830;
    0.830;
];
V_capt = [
    14.8100;
    16.2700;
    16.6700;
    17.1860;
    17.5700;
    18.3000;
    18.7016;
    19.5300;
    19.8750;
    20.1500;
    20.3500;
    20.4800;
    20.6600;
    20.7850;
    21.0130;
    21.1800;
    21.3600;
    21.3600;
];
t_ms = [
    3.11418;
    3.74697;
    4.15935;
    4.79214;
    5.21163;
    6.05061;
    6.67629;
    6.88248;
    6.88248;
    7.08867;
    7.09578;
    7.08156;
    7.09578;
    7.09578;
    7.09578;
    7.08156;
    7.08867;
    7.08867;
];
if useCapturedData
    I_data = I_capt;
    V_data = V_capt;
else
    I_data = I_raw;
    V_data = V_raw;
end

%% 2) Plot settings

xLabelCurrent = '$I_{out}^{OT}\,[\mathrm{A}]$';
xLabelVoltage = '$V_{out}^{OT}\,[\mathrm{V}]$';
fontSize = 30;
lineWidth = 2;

%% 3) Calculation time versus current

figure('Color','w','Position',[100 100 900 500]);
plot(I_data, t_ms, '-ko', ...
    'LineWidth', lineWidth, ...
    'MarkerSize', 6, ...
    'MarkerFaceColor', 'k');
grid on;
xlabel(xLabelCurrent, 'Interpreter','latex','FontSize',fontSize);
ylabel('$t_{calc}\,[\mathrm{ms}]$', 'Interpreter','latex','FontSize',fontSize);
set(gca, 'FontSize', fontSize);
set(findall(gcf, '-property','FontName'), 'FontName', 'Times New Roman');

%% 4) Calculation time versus voltage

figure('Color','w','Position',[100 100 900 500]);
plot(V_data, t_ms, '-ko', ...
    'LineWidth', lineWidth, ...
    'MarkerSize', 6, ...
    'MarkerFaceColor', 'k');
grid on;
xlabel(xLabelVoltage, 'Interpreter','latex','FontSize',fontSize);
ylabel('$t_{calc}\,[\mathrm{ms}]$', 'Interpreter','latex','FontSize',fontSize);
set(gca, 'FontSize', fontSize);
set(findall(gcf, '-property','FontName'), 'FontName', 'Times New Roman');
