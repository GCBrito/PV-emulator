
clear; clc; close all;

%% 1. Data Input 

% Data format for all matrices:
% [Vout_OT [V], Iout_OT [A], Voltage error [%], Current error [%]]
%
% Each row corresponds to one operating point.

% --- Canadian Solar CS6P-250P ---
data_CS = [
    22.642405, 6.816042, 1.809375000, 0.058035191;
    23.414888, 6.793228, 1.981219512, 4.723309958;
    24.056635, 6.770828, 1.805480322, 1.299883382;
    24.833183, 6.729959, 1.775340164, 1.752423358;
    25.833637, 6.648664, 1.627210858, 1.937109145;
    26.700890, 6.539245, 1.485708856, 1.813138138;
    27.416254, 6.408448, 1.353988909, 2.310243902;
    28.057871, 6.251021, 1.218870851, 1.867802198;
    28.952494, 5.944597, 1.055825480, 2.065947282;
    29.632019, 5.614259, 0.857791014, 1.676725044;
    30.294250, 5.214557, 0.678796943, 1.982011278;
    30.687252, 4.894783, 0.548007864, 0.923360825;
    30.963587, 4.659785, 0.400736057, 0.644243070;
    31.429310, 4.263731, 0.316980530, 0.843465116;
    31.895897, 3.582153, 0.081258237, 0.060139665;
    32.502819, 2.653505, 0.175617322, 0.893726236;
    33.319305, 1.404204, 0.569068935, 7.191145038;
    34.069290, 0.256652, 0.817205240, 220.8150000
];

% --- Kyocera KB260-6BPA ---
data_KB260 = [
    23.384464, 6.631636, 1.760069626, 1.753540741;
    24.569763, 6.605496, 1.444108175, 0.669233083;
    25.593412, 6.571879, 1.440396354, 2.058435171;
    26.668896, 6.513233, 1.364104903, 2.496511976;
    27.661718, 6.433278, 1.213750457, 2.378179059;
    28.555620, 6.326908, 1.153453773, 1.450031153;
    29.455330, 6.173790, 0.943557231, 3.232131661;
    30.391579, 5.944713, 0.767834881, 3.338000000;
    31.405592, 5.583774, 0.562254243, 2.039052632;
    31.972301, 5.297785, 0.415518216, 2.254889299;
    32.411503, 5.064976, 0.283115718, 1.841550388;
    32.913673, 4.688922, 0.163338405, 1.699748428;
    33.195103, 4.473986, 0.105859469, 1.018008850;
    33.651192, 4.120906, 0.055859816, 0.701060241;
    34.301086, 3.119885, 0.258546089, 0.317845659;
    35.331486, 1.532781, 0.642615298, 5.709034483;
    36.170616, 0.240290, 0.983805092, 200.3625000
];

% --- Kyocera KC200GT ---
data_KC200 = [
    12.533917, 4.225679, 2.568878887, 3.065341463;
    14.705232, 4.217923, 2.261696801, 2.588383372;
    17.256174, 4.190397, 1.566650971, 1.217318841;
    18.053396, 4.175287, 1.537660292, 2.085256724;
    18.770275, 4.161449, 1.296681058, 1.498756098;
    19.421268, 4.134554, 1.257914494, 0.131545894;
    20.231714, 4.094382, 1.007059411, 1.095851852;
    21.174362, 4.022306, 0.926415634, 0.683802469;
    21.890760, 3.941084, 0.739806719, 0.537857143;
    22.502890, 3.847137, 0.684071588, 0.185859375;
    23.111118, 3.723294, 0.570574413, 0.629567568;
    23.679567, 3.574845, 0.550178344, 0.699859155;
    24.440063, 3.309085, 0.369868583, 1.195259939;
    24.939966, 3.073253, 0.160506024, 0.433104575;
    25.292833, 2.895774, 0.169635644, 1.963873239;
    25.808134, 2.549612, 0.007229756, 2.806935484;
    26.232281, 2.259961, 0.143582033, 2.725500000;
    26.693003, 1.693082, 0.362064203, 5.160372671;
    27.100122, 1.176412, 0.550011009, 15.33450980;
    27.427097, 0.761451, 0.698417813, 33.58789474;
    27.814808, 0.269414, 0.838474153, 349.0233333
];

% --- Kyocera KC85TS (R Load) ---
data_KC85_R = [
    4.290837, 2.135176, 4.146529126, 10.63088083;
    10.134062, 2.133931, 1.952334004, 4.094195122;
    12.669167, 2.131800, 1.515761218, 2.985507246;
    13.826797, 2.129210, 1.147015362, 1.390952381;
    15.001449, 2.120053, 0.883987895, 2.915194175;
    15.612347, 2.109718, 0.595019330, 3.417549020;
    16.172581, 2.092691, 0.763744548, 3.088226601;
    16.606222, 2.071157, 0.704802911, 1.032048780;
    17.095009, 2.033670, 0.618063567, 1.683500000;
    17.511776, 1.984805, 0.584583573, 2.309536082;
    18.038937, 1.889737, 0.383622705, 4.405359116;
    18.432518, 1.777257, 0.394978214, 4.544529412;
    18.831503, 1.632836, 0.274243876, 6.721307190;
    19.143658, 1.463260, 0.123734310, 8.389629630;
    19.473244, 1.284213, 0.034681725, 12.65026316;
    19.819204, 0.873231, 0.205417925, 28.41632353;
    20.202765, 0.406569, 0.380843195, 125.8716667;
    20.301470, 0.286478, 0.482990196, 472.9560000
];

% --- Kyocera KC85TS (RL Load) ---
data_KC85_RL = [
    4.212106, 2.135191, 5.302650000, 6.759550000;
    10.206792, 2.133915, 2.170090090, 5.639356436;
    12.650108, 2.131843, 1.607293173, 4.502107843;
    13.811143, 2.129245, 1.329002201, 1.877751196;
    14.974236, 2.120420, 1.040728745, 4.454187192;
    15.629344, 2.109370, 0.834477419, 3.400490196;
    16.181227, 2.092373, 0.817613707, 2.567303922;
    16.588015, 2.072266, 0.838996960, 2.082068966;
    17.113245, 2.031851, 0.606966490, 2.618737374;
    17.581903, 1.974965, 0.640543789, 3.945526316;
    18.077917, 1.881181, 0.488699277, 4.510055556;
    18.393583, 1.790574, 0.292164667, 3.501387283;
    18.813585, 1.642570, 0.178833866, 6.660389610;
    19.156639, 1.456207, 0.139252483, 8.672164179;
    19.463346, 1.289590, 0.017194245, 13.12192982;
    19.813057, 0.880713, 0.186110831, 29.51661765;
    20.199619, 0.410394, 0.445446033, 127.9966667;
    20.309837, 0.276297, 0.490754532, 452.5940000
];

%% 2. Figure Settings
lineW = 1.5; fLabelSize = 30; fLegendSize = 22;

figure(1);
semilogy(data_CS(:,2), data_CS(:,4), '-o', 'LineWidth', lineW); hold on;
semilogy(data_KC200(:,2), data_KC200(:,4), '-s', 'LineWidth', lineW);
semilogy(data_KB260(:,2), data_KB260(:,4), '-^', 'LineWidth', lineW);

xlabel('$I_{out}^{OT}\,[\mathrm{A}]$', 'Interpreter','latex','FontSize',fLabelSize);
ylabel('$E_{I_{out}}\,[\%]$', 'Interpreter','latex','FontSize',fLabelSize);

lgd1 = legend({'$\mathrm{Canadian\ Solar\ CS6P\!-\!250P}$', '$\mathrm{Kyocera\ KC200GT}$', '$\mathrm{Kyocera\ KB260\!-\!6BPA}$'}, ...
               'Interpreter','latex','Location','northeast');
set(lgd1,'FontSize',fLegendSize, 'FontName', 'Times New Roman');
set(gca,'FontSize',fLabelSize, 'TickLabelInterpreter','latex');

ax = gca;
ax.YScale = 'log';
ax.YLim = [1e-3 1e3];
ax.YTick = 10.^(-3:3);
ax.YAxis.MinorTickValues = sort([ ...
    (2:9)*1e-3, ...
    (2:9)*1e-2, ...
    (2:9)*1e-1, ...
    (2:9)*1e0, ...
    (2:9)*1e1, ...
    (2:9)*1e2 ]);

ax.YGrid = 'on';
ax.YMinorGrid = 'on';
ax.XGrid = 'on';
ax.XMinorGrid = 'off';

ax.GridLineStyle = '-';
ax.MinorGridLineStyle = ':';
ax.GridAlpha = 0.25;
ax.MinorGridAlpha = 0.35;
%% Figure 2: Voltage Error (Original 3 Panels)
figure(2);

semilogy(data_CS(:,1), data_CS(:,3), '-o', 'LineWidth', lineW); hold on;
semilogy(data_KC200(:,1), data_KC200(:,3), '-s', 'LineWidth', lineW);
semilogy(data_KB260(:,1), data_KB260(:,3), '-^', 'LineWidth', lineW);

xlabel('$V_{out}^{OT}\,[\mathrm{V}]$', 'Interpreter','latex','FontSize',fLabelSize);
ylabel('$E_{V_{out}}\,[\%]$', 'Interpreter','latex','FontSize',fLabelSize);

lgd2 = legend({ ...
    '$\mathrm{Canadian\ Solar\ CS6P\!-\!250P}$', ...
    '$\mathrm{Kyocera\ KC200GT}$', ...
    '$\mathrm{Kyocera\ KB260\!-\!6BPA}$'}, ...
    'Interpreter','latex', 'Location','northeast');

set(lgd2, 'FontSize', fLegendSize, 'FontName', 'Times New Roman');
set(gca, 'FontSize', fLabelSize, 'TickLabelInterpreter', 'latex');

%% 4. Configure log grid on Y axis
ax = gca;

% Log scale and limits
ax.YScale = 'log';
ax.YLim = [1e-3 1e3];

% Major ticks: powers of 10
ax.YTick = 10.^(-3:3);

% Minor ticks: 2 to 9 within each decade
ax.YAxis.MinorTickValues = sort([ ...
    (2:9)*1e-3, ...
    (2:9)*1e-2, ...
    (2:9)*1e-1, ...
    (2:9)*1e0, ...
    (2:9)*1e1, ...
    (2:9)*1e2 ]);

% Grid only where wanted
ax.YGrid = 'on';
ax.YMinorGrid = 'on';
ax.XGrid = 'on';
ax.XMinorGrid = 'off';

% Grid appearance
ax.GridLineStyle = '-';
ax.MinorGridLineStyle = ':';
ax.GridAlpha = 0.25;
ax.MinorGridAlpha = 0.35;

%% Figure 3: Current Error (KC85TS: R vs RL)

figure(3);

semilogy(data_KC85_R(:,2),  data_KC85_R(:,4),  '-d', 'LineWidth', lineW); hold on;
semilogy(data_KC85_RL(:,2), data_KC85_RL(:,4), '-x', 'LineWidth', lineW);

xlabel('$I_{out}^{OT}\,[\mathrm{A}]$', 'Interpreter','latex','FontSize',fLabelSize);
ylabel('$E_{I_{out}}\,[\%]$', 'Interpreter','latex','FontSize',fLabelSize);

lgd3 = legend({ ...
    '$\mathrm{Kyocera\ KC85TS\ (R\ load)}$', ...
    '$\mathrm{Kyocera\ KC85TS\ (RL\ load)}$'}, ...
    'Interpreter','latex', 'Location','northeast');

set(lgd3, 'FontSize', fLegendSize, 'FontName', 'Times New Roman');
set(gca, 'FontSize', fLabelSize, 'TickLabelInterpreter', 'latex');

%% --- Log grid configuration (10^-2 to 10^3) ---
ax = gca;

ax.YScale = 'log';
ax.YLim = [1e-2 1e3];

% Major ticks (decades)
ax.YTick = 10.^(-2:3);

% Minor ticks (2..9 per decade)
ax.YAxis.MinorTickValues = sort([ ...
    (2:9)*1e-2, ...
    (2:9)*1e-1, ...
    (2:9)*1e0, ...
    (2:9)*1e1, ...
    (2:9)*1e2 ]);

% Grid
ax.YGrid = 'on';
ax.YMinorGrid = 'on';
ax.XGrid = 'on';
ax.XMinorGrid = 'off';

% Appearance
ax.GridLineStyle = '-';
ax.MinorGridLineStyle = ':';
ax.GridAlpha = 0.25;
ax.MinorGridAlpha = 0.35;

%% Figure 4: Voltage Error (KC85TS: R vs RL)

figure(4);

semilogy(data_KC85_R(:,1),  data_KC85_R(:,3),  '-d', 'LineWidth', lineW); hold on;
semilogy(data_KC85_RL(:,1), data_KC85_RL(:,3), '-x', 'LineWidth', lineW);

%% --- Limites do eixo X (margem desejada) ---
xlim([0 25]);   % <<< AQUI está a margem

xlabel('$V_{out}^{OT}\,[\mathrm{V}]$', 'Interpreter','latex','FontSize',fLabelSize);
ylabel('$E_{V_{out}}\,[\%]$', 'Interpreter','latex','FontSize',fLabelSize);

lgd4 = legend({ ...
    '$\mathrm{Kyocera\ KC85TS\ (R\ load)}$', ...
    '$\mathrm{Kyocera\ KC85TS\ (RL\ load)}$'}, ...
    'Interpreter','latex', 'Location','northeast');

set(lgd4, 'FontSize', fLegendSize, 'FontName', 'Times New Roman');
set(gca, 'FontSize', fLabelSize, 'TickLabelInterpreter', 'latex');

%% --- Configuração log no eixo Y (10^-2 a 10^3) ---
ax = gca;

ax.YScale = 'log';
ax.YLim = [1e-2 1e3];

% Ticks principais
ax.YTick = 10.^(-2:3);

% Ticks menores (grade log verdadeira)
ax.YAxis.MinorTickValues = sort([ ...
    (2:9)*1e-2, ...
    (2:9)*1e-1, ...
    (2:9)*1e0, ...
    (2:9)*1e1, ...
    (2:9)*1e2 ]);

% Grade
ax.YGrid = 'on';
ax.YMinorGrid = 'on';
ax.XGrid = 'on';
ax.XMinorGrid = 'off';

% Estilo
ax.GridLineStyle = '-';
ax.MinorGridLineStyle = ':';
ax.GridAlpha = 0.25;
ax.MinorGridAlpha = 0.35;