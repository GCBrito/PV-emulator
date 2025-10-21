clear, clc, close all;

%% 1) Display of the three curves + experimental points
figure('Color','w','Position',[100 100 900 500]);
hold on; grid on;

% Single diode emulator - KC85TS CS6P-250P (S = 765 W/m^2 and T = 44.5°C)
V_emulator = [
    0
    22.15;
    23.53;
    25.08;
    27.04;
    27.92;
    29.02;
    30.23;
    32.15;
    33.59;
    33.81;
    34;
];

I_emulator = [
    6.82;
    6.81;
    6.79;
    6.66;
    6.43;
    6.07;
    5.59;
    4.95;
    2.53;
    0.71;
    0.27;
    0;
];

h1 = plot(V_emulator,  I_emulator,   '-b',  'LineWidth',2);     % plot single diode emulator

% Real panel data - KC85TS CS6P-250P (S = 765 W/m^2 and T = 44.5°C)
V_pannel = [
    0.047996217771366645;
    2.5473136313978317;
    4.686489351002983;
    7.605990366562553;
    10.181947580363017;
    12.610309040298777;
    14.79342630980413;
    17.20541669131795;
    19.415489834424953;
    21.816632013769198;
    24.03214755291164;
    26.237334370932288;
    28.26322957001055;
    29.749962170030145;
    31.046437882063465;
    31.92801141798128;
    32.689215781035756;
    33.25368190774834;
    33.59447586154349;
    34.38443600476654
];

I_pannel = [
    7.007823910286788;
    7.0079635369891475;
    7.001379133258391;
    7.003776871281191;
    6.93017776227956;
    6.9347826983712295;
    6.841049913320265;
    6.845653933494394;
    6.8524813120762555;
    6.83697299652661;
    6.848269952562965;
    6.68079538185797;
    6.2719700120743145;
    5.5323882664885256;
    4.567097566839892;
    3.6531803358379182;
    2.83534576727134;
    2.1046510438953376;
    1.3538320938536454;
    0.08683712275429656
];

h2 = plot(V_pannel, I_pannel, 'ko', ...      % plot real panel data
          'MarkerSize',6, 'MarkerFaceColor','k');

xlabel('Voltage (V)', 'FontSize', 30);
ylabel('Current (A)', 'FontSize', 30); 
legend([h1 h2], ...
       {'PV Emulator – Single Diode Model', 'Experimental Results – Real PV Panel'}, ...
       'Location','best');

% Axis ticks font configuration
set(gca, 'FontSize', 30); % Adjusted font size based on your image

% Apply Times New Roman font to all text elements
set(findall(gcf, '-property','FontName'),'FontName','Times New Roman');
set(findall(gca, '-property','FontName'),'FontName','Times New Roman');
