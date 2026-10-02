clear; clc; close all;

% Compares three PV emulators with measured I-V curves from three PV panels.
% RMSE, nRMSE and MAPE are evaluated on a common voltage grid.
% Relative errors are also calculated for Isc, Voc, Vmpp, Impp and Pmpp.
% MAPE excludes reference values below 10% of Isc or Pmpp.
% Results are displayed in MATLAB only; no files or folders are created.

%% 1) Settings
fontAxisLabel = 30;
fontTicks     = 26;
fontLegend    = 15;
fontCallout   = 15;
lineWidthEmu  = 2.2;
markerSizeRef = 6;
% Number of interpolation points used for common-grid metrics.
nCommonGrid = 1000;

%% 2) Input data
% Real-panel measurements are stored in V_panel and I_panel.
% makeEmu inputs: full name, short name, RGB color, line style, V vector, I vector.

panels = struct([]);

%% Panel 1: ME Solar MESM-50W
panels(1).name = 'ME Solar MESM-50W';
panels(1).condition = 'G = 1000 W/m^2, T = 25 degC';
panels(1).V_panel = [
    0.1345090408816909;
    3.9418622020801637;
    7.1102974135229875;
    10.316096247185497;
    13.219249737750683;
    14.493349274149194;
    15.386339859898216;
    16.114930504317343;
    16.794948437178412;
    17.542220892707384;
    18.26333881312798;
    18.943356748216104;
    19.39172022153349;
    19.65700194241112;
    19.90733821529171;
    20.228665369498877;
    20.51636526682619;
    20.882528769478625;
    21.1702286645789;
    21.446719472567864;
    21.71200119344551
];
panels(1).I_panel = [
    3.0107092341738317;
    2.9969690844678114;
    2.9832289352435426;
    2.970297029807906;
    2.9614063447606873;
    2.9533239039838524;
    2.9533239035021017;
    2.9525156597134687;
    2.932309557771381;
    2.880581935547086;
    2.7868256213795988;
    2.643766417317217;
    2.503940189854619;
    2.3099616083200765;
    2.122448979503352;
    1.8508789653549917;
    1.5663770457709956;
    1.1840775911493435;
    0.8163265305040968;
    0.4598908867537168;
    0.030713275530074746
];
panels(1).emu(1) = makeEmu('Open-source PV Emulator - Single Diode Model', 'Single diode', [0 0 1], '-', [
    0;
    5.922083;
    12.551384;
    13.841137;
    15.125275;
    15.540461;
    15.992638;
    16.356771;
    16.715565;
    17.007088;
    17.418226;
    17.789234;
    18.138689;
    18.464882;
    18.834526;
    19.264559;
    19.65144;
    19.971132;
    20.552086;
    21.464041;
    21.969181;
    22.3
], [
    3.03;
    3.026813;
    3.019652;
    3.0118;
    2.995558;
    2.984576;
    2.96906;
    2.951598;
    2.929032;
    2.905429;
    2.862687;
    2.811599;
    2.749611;
    2.676852;
    2.569205;
    2.419445;
    2.236233;
    2.070245;
    1.553807;
    0.743124;
    0.29408;
    0
]);
panels(1).emu(2) = makeEmu('Open-source PV Emulator - Simplified Exp. Model', 'Simplified exp.', [1 0 0], '-', [
    0;
    8.08962;
    11.45164;
    14.526568;
    16.046844;
    16.94203;
    17.753181;
    18.2904;
    19.387537;
    20.072838;
    20.594084;
    21.106194;
    21.724693;
    22.3
], [
    3.03;
    3.028738;
    3.021295;
    2.99427;
    2.949897;
    2.878464;
    2.813736;
    2.679999;
    2.389157;
    2.202046;
    1.848619;
    1.501389;
    0.7471;
    0
]);
panels(1).emu(3) = makeEmu('Commercial PV emulator', 'Commercial', [1 0 1], '-', [
    1.1;
    2.7;
    4.6;
    6.5;
    8.5702;
    9.4;
    13.06;
    13.8;
    16.2;
    17.6774;
    18.83;
    19.4;
    20.1742;
    20.939;
    21.12;
    21.485;
    21.7584;
    21.93;
    22.1151
], [
    3.0261;
    3.0203;
    3.0129;
    3.0024;
    2.991;
    2.985;
    2.95;
    2.9382;
    2.88;
    2.818;
    2.59;
    2.205;
    1.61;
    1.033;
    0.89;
    0.6146;
    0.4062;
    0.2774;
    0.1343
]);

%% Panel 2: Renogy RNG-50DB-H
panels(2).name = 'Renogy RNG-50DB-H';
panels(2).condition = 'G = 1000 W/m^2, T = 25 degC';
panels(2).V_panel = [
    0.061238999661659355;
    3.0378744826113886;
    6.06593628879359;
    9.019207316748718;
    11.551928051656516;
    12.757528558290698;
    13.855704699997231;
    14.771648303222042;
    15.865224726830606;
    16.729902361354295;
    17.664819810269254;
    18.450378783493765;
    19.100584026063657;
    19.605999152099663;
    20.04147065780234;
    20.416281613412224;
    20.80544927623029;
    21.20426937792104;
    21.462738787922408;
    21.777746419471736;
    22.02257696007481;
    22.272296443566514;
    22.475133634590215;
    22.579352135807298
];
panels(2).I_panel = [
    2.9194630872527965;
    2.9194630872527965;
    2.915548098563706;
    2.915548098563706;
    2.913870246173136;
    2.9149888143223963;
    2.9077181208521656;
    2.8987695749913653;
    2.879753914787145;
    2.848993288182294;
    2.794742729109542;
    2.7125279639717688;
    2.6045861294008548;
    2.4854586130039786;
    2.342281879231173;
    2.1851230434259357;
    1.9737136465477885;
    1.713087248101959;
    1.4787472036639384;
    1.170022371466326;
    0.8747203580599141;
    0.5447427290262028;
    0.23937360194308832;
    0.012304250975298636
];
panels(2).emu(1) = makeEmu('Open-source PV Emulator - Single Diode Model', 'Single diode', [0 0 1], '-', [
    0;
    5.894531;
    12.125871;
    13.816704;
    16.145973;
    16.715784;
    17.231913;
    17.664825;
    18.180428;
    18.727219;
    19.188766;
    19.602093;
    19.970377;
    20.420042;
    20.759869;
    21.242069;
    21.980999;
    22.318256;
    22.6
], [
    2.92;
    2.915734;
    2.909352;
    2.902917;
    2.874516;
    2.856222;
    2.831469;
    2.801653;
    2.751547;
    2.673501;
    2.579515;
    2.459445;
    2.335244;
    2.107949;
    1.91376;
    1.412265;
    0.643768;
    0.293014;
    0
]);
panels(2).emu(2) = makeEmu('Open-source PV Emulator - Simplified Exp. Model', 'Simplified exp.', [1 0 0], '-', [
    0;
    6.817072;
    11.080575;
    15.057276;
    16.803661;
    17.610468;
    18.552885;
    19.459518;
    20.528112;
    21.02478;
    21.487143;
    22.002245;
    22.6
], [
    2.92;
    2.919855;
    2.918193;
    2.894286;
    2.843843;
    2.780185;
    2.695215;
    2.44177;
    2.114295;
    1.793764;
    1.483582;
    0.796883;
    0
]);
panels(2).emu(3) = makeEmu('Commercial PV emulator', 'Commercial', [1 0 1], '-', [
    0;
    1;
    1.2329;
    3.09;
    5.12;
    8.5471;
    10.24;
    12.12;
    14.44;
    15.2639;
    16.61;
    17.83;
    18.35;
    19.41;
    19.69;
    20.02;
    20.6976;
    21.1701;
    21.5176;
    21.8835;
    22.0688;
    22.24;
    22.5
], [
    2.92;
    2.917;
    2.9162;
    2.9106;
    2.9031;
    2.8875;
    2.8768;
    2.8623;
    2.8358;
    2.8236;
    2.795;
    2.756;
    2.71;
    2.51;
    2.289;
    2.0328;
    1.5;
    1.127;
    0.852;
    0.5613;
    0.4157;
    0.2815;
    0.0755
]);

%% Panel 3: Shell Solar SQ150-PC
panels(3).name = 'Shell Solar SQ150-PC';
panels(3).condition = 'G = 1000 W/m^2, T = 25 degC';
panels(3).V_panel = [
    0.07780604683747905;
    4.2189436898142905;
    7.660223341663279;
    10.863995324097095;
    14.176868703729975;
    17.155893827390162;
    20.141363211697026;
    23.036948208895787;
    24.995167557029617;
    26.42690527335035;
    27.81375668458678;
    29.001602344647385;
    30.36284434236471;
    31.512299834342564;
    32.59122368036738;
    33.76659988606946;
    34.743098058886865;
    35.65550884640268;
    36.446045441351224;
    37.16601846173125;
    37.82826853208371;
    38.47780170128618;
    39.10175339833374;
    39.68073097512157;
    40.29205756197019;
    40.91626468796069;
    41.54063145321372;
    42.171478402636254;
    42.8474591077566;
    43.38200934569101
];
panels(3).I_panel = [
    4.804197165221235;
    4.795540075268057;
    4.798166796403276;
    4.789358636164255;
    4.791964661586494;
    4.79244477341078;
    4.786709596059194;
    4.780959931956091;
    4.775059198186533;
    4.773217832998629;
    4.756864468979553;
    4.734262703388705;
    4.699256226304568;
    4.642458468761035;
    4.56182008609613;
    4.442863222890885;
    4.284504230031027;
    4.096089321638853;
    3.878645240109937;
    3.6456489653261293;
    3.397102567352327;
    3.1164364242871767;
    2.809864790583181;
    2.5115743462965194;
    2.1490537237230516;
    1.7761746301469223;
    1.3618533741500927;
    0.9319923422019132;
    0.4524079577975444;
    0.02045929702601157
];
panels(3).emu(1) = makeEmu('Open-source PV Emulator - Single Diode Model', 'Single diode', [0 0 1], '-', [
    0;
    23.428795;
    25.957676;
    27.723593;
    28.600252;
    29.371901;
    30.186039;
    31.42337;
    32.285755;
    33.184868;
    33.870052;
    34.74757;
    35.65852;
    36.230217;
    36.874702;
    37.490471;
    38.099937;
    39.200233;
    40.635674;
    42.10854;
    43.159355;
    43.4
], [
    4.8;
    4.771824;
    4.754732;
    4.736563;
    4.721648;
    4.70373;
    4.679128;
    4.624589;
    4.569901;
    4.492322;
    4.415938;
    4.290498;
    4.119548;
    3.97847;
    3.803591;
    3.594777;
    3.382138;
    2.680002;
    1.764004;
    0.824125;
    0.153565;
    0
]);
panels(3).emu(2) = makeEmu('Open-source PV Emulator - Simplified Exp. Model', 'Simplified exp.', [1 0 0], '-', [
    0;
    25.567003;
    27.200571;
    28.811605;
    30.372927;
    32.341831;
    34.1399;
    36.391037;
    38.724979;
    40.008648;
    41.242683;
    42.378609;
    42.994045;
    43.4
], [
    4.8;
    4.747;
    4.733705;
    4.687959;
    4.643624;
    4.51567;
    4.370661;
    3.898574;
    3.404761;
    2.759869;
    2.083102;
    0.986255;
    0.391991;
    0
]);
panels(3).emu(3) = makeEmu('Commercial PV emulator', 'Commercial', [1 0 1], '-', [
    1.55;
    3.4;
    8.22;
    11.75;
    14.85;
    17.9;
    19.52;
    21.8;
    24.55;
    27.47;
    30.13;
    33.09;
    35.9;
    36.8505;
    37.8;
    38.5;
    38.8934;
    39.3103;
    39.5975;
    39.99;
    40.2553;
    40.7371;
    41.228;
    41.7376;
    42.0248;
    42.4047;
    42.6456
], [
    4.7949;
    4.7884;
    4.7695;
    4.7529;
    4.735;
    4.715;
    4.7027;
    4.68;
    4.65;
    4.6095;
    4.555;
    4.4518;
    4.1523;
    3.7759;
    3.23;
    2.8289;
    2.6045;
    2.3639;
    2.1974;
    1.9679;
    1.8173;
    1.5363;
    1.252;
    0.9905;
    0.7912;
    0.571;
    0.4338
]);

%% 3) Initialize result arrays
nPanels = numel(panels);
nEmu = numel(panels(1).emu);
nRows = nPanels * nEmu;
PanelName = strings(nRows,1);
Condition = strings(nRows,1);
Emulator = strings(nRows,1);
Npoints = zeros(nRows,1);
CommonGrid_Vmin_V = zeros(nRows,1);
CommonGrid_Vmax_V = zeros(nRows,1);
Isc_realPV_A = zeros(nRows,1);
Isc_emu_A = zeros(nRows,1);
Isc_error_percent = zeros(nRows,1);
Voc_realPV_V = zeros(nRows,1);
Voc_emu_V = zeros(nRows,1);
Voc_error_percent = zeros(nRows,1);
Vmpp_realPV_V = zeros(nRows,1);
Vmpp_emu_V = zeros(nRows,1);
Vmpp_error_percent = zeros(nRows,1);
Impp_realPV_A = zeros(nRows,1);
Impp_emu_A = zeros(nRows,1);
Impp_error_percent = zeros(nRows,1);
Pmpp_realPV_W = zeros(nRows,1);
Pmpp_emu_W = zeros(nRows,1);
Pmpp_error_percent = zeros(nRows,1);
RMSE_I_A = zeros(nRows,1);
nRMSE_I_percent_Isc = zeros(nRows,1);
MAPE_I_percent = zeros(nRows,1);
RMSE_P_W = zeros(nRows,1);
nRMSE_P_percent_Pmpp = zeros(nRows,1);
MAPE_P_percent = zeros(nRows,1);
Npoints_MAPE_I = zeros(nRows,1);
Npoints_MAPE_P = zeros(nRows,1);

%% 4) Calculate metrics and plot I-V curves
row = 0;
for p = 1:nPanels
    panel = panels(p);
    [Vref, Iref] = prepareCurve(panel.V_panel, panel.I_panel);
    Isc_ref = estimateIscFromCurve(Vref, Iref);
    Voc_ref = estimateVocFromCurve(Vref, Iref);
    [Vmpp_ref, Impp_ref, Pmpp_ref] = findMPPFromCurve(Vref, Iref, 3000);
    % Build one voltage range shared by the real panel and all emulators.
    VminCandidates = zeros(nEmu + 1, 1);
    VmaxCandidates = zeros(nEmu + 1, 1);
    [Vref_metric, Iref_metric] = augmentCurveWithIscVoc(Vref, Iref);
    VminCandidates(1) = min(Vref_metric);
    VmaxCandidates(1) = max(Vref_metric);
    for kk = 1:nEmu
        [Vtmp, Itmp] = prepareCurve(panel.emu(kk).V, panel.emu(kk).I);
        [Vtmp_metric, ~] = augmentCurveWithIscVoc(Vtmp, Itmp);
        VminCandidates(kk + 1) = min(Vtmp_metric);
        VmaxCandidates(kk + 1) = max(Vtmp_metric);
    end
    Vgrid_min = max(VminCandidates);
    Vgrid_max = min(VmaxCandidates);
    if Vgrid_max <= Vgrid_min
        error('Invalid common voltage grid for panel %s.', panel.name);
    end
    Vgrid_common = linspace(Vgrid_min, Vgrid_max, nCommonGrid).';
    Iref_grid = interp1(Vref_metric, Iref_metric, Vgrid_common, 'pchip');
    Iref_grid(Iref_grid < 0) = 0;
    Pref_grid = Vgrid_common .* Iref_grid;
    panelMetrics = struct([]);
    for k = 1:nEmu
        row = row + 1;
        [Vemu, Iemu] = prepareCurve(panel.emu(k).V, panel.emu(k).I);
        [Vemu_metric, Iemu_metric] = augmentCurveWithIscVoc(Vemu, Iemu);
        % Interpolate each emulator on the same voltage grid as the reference.
        Iemu_grid = interp1(Vemu_metric, Iemu_metric, Vgrid_common, 'pchip');
        Iemu_grid(Iemu_grid < 0) = 0;
        % Current RMSE, nRMSE and filtered MAPE.
        eI = Iemu_grid - Iref_grid;
        RMSE_I = sqrt(mean(eI.^2));
        nRMSE_I = 100 * RMSE_I / abs(Isc_ref);
        I_threshold = 0.10 * abs(Isc_ref);
        validMAPE_I = Iref_grid > I_threshold;
        MAPE_I = mean(abs((Iemu_grid(validMAPE_I) - Iref_grid(validMAPE_I)) ./ Iref_grid(validMAPE_I))) * 100;
        % Power RMSE, nRMSE and filtered MAPE.
        Pemu_grid = Vgrid_common .* Iemu_grid;
        eP = Pemu_grid - Pref_grid;
        RMSE_P = sqrt(mean(eP.^2));
        nRMSE_P = 100 * RMSE_P / abs(Pmpp_ref);
        P_threshold = 0.10 * abs(Pmpp_ref);
        validMAPE_P = Pref_grid > P_threshold;
        MAPE_P = mean(abs((Pemu_grid(validMAPE_P) - Pref_grid(validMAPE_P)) ./ Pref_grid(validMAPE_P))) * 100;
        % Characteristic-point values and relative errors.
        Isc_emu = estimateIscFromCurve(Vemu, Iemu);
        Voc_emu = estimateVocFromCurve(Vemu, Iemu);
        [Vmpp_emu, Impp_emu, Pmpp_emu] = findMPPFromCurve(Vemu, Iemu, 3000);
        Isc_err = relativeErrorPercent(Isc_emu, Isc_ref);
        Voc_err = relativeErrorPercent(Voc_emu, Voc_ref);
        Vmpp_err = relativeErrorPercent(Vmpp_emu, Vmpp_ref);
        Impp_err = relativeErrorPercent(Impp_emu, Impp_ref);
        Pmpp_err = relativeErrorPercent(Pmpp_emu, Pmpp_ref);
        % Store values used in the results table.
        PanelName(row) = string(panel.name);
        Condition(row) = string(panel.condition);
        Emulator(row) = string(panel.emu(k).shortName);
        Npoints(row) = numel(Vgrid_common);
        CommonGrid_Vmin_V(row) = Vgrid_min;
        CommonGrid_Vmax_V(row) = Vgrid_max;
        Isc_realPV_A(row) = Isc_ref;
        Isc_emu_A(row) = Isc_emu;
        Isc_error_percent(row) = Isc_err;
        Voc_realPV_V(row) = Voc_ref;
        Voc_emu_V(row) = Voc_emu;
        Voc_error_percent(row) = Voc_err;
        Vmpp_realPV_V(row) = Vmpp_ref;
        Vmpp_emu_V(row) = Vmpp_emu;
        Vmpp_error_percent(row) = Vmpp_err;
        Impp_realPV_A(row) = Impp_ref;
        Impp_emu_A(row) = Impp_emu;
        Impp_error_percent(row) = Impp_err;
        Pmpp_realPV_W(row) = Pmpp_ref;
        Pmpp_emu_W(row) = Pmpp_emu;
        Pmpp_error_percent(row) = Pmpp_err;
        RMSE_I_A(row) = RMSE_I;
        nRMSE_I_percent_Isc(row) = nRMSE_I;
        MAPE_I_percent(row) = MAPE_I;
        RMSE_P_W(row) = RMSE_P;
        nRMSE_P_percent_Pmpp(row) = nRMSE_P;
        MAPE_P_percent(row) = MAPE_P;
        Npoints_MAPE_I(row) = sum(validMAPE_I);
        Npoints_MAPE_P(row) = sum(validMAPE_P);
        % Store values used in the I-V callouts.
        panelMetrics(k).shortName = panel.emu(k).shortName;
        panelMetrics(k).Isc_emu = Isc_emu;
        panelMetrics(k).Voc_emu = Voc_emu;
        panelMetrics(k).Vmpp_emu = Vmpp_emu;
        panelMetrics(k).Impp_emu = Impp_emu;
        panelMetrics(k).Pmpp_emu = Pmpp_emu;
        panelMetrics(k).Isc_err = Isc_err;
        panelMetrics(k).Voc_err = Voc_err;
        panelMetrics(k).Vmpp_err = Vmpp_err;
        panelMetrics(k).Impp_err = Impp_err;
        panelMetrics(k).Pmpp_err = Pmpp_err;
    end
    plotPanelIV(panel, Vref, Iref, Isc_ref, Voc_ref, Vmpp_ref, Impp_ref, panelMetrics, ...
        fontAxisLabel, fontTicks, fontLegend, fontCallout, lineWidthEmu, markerSizeRef);
end

%% 5) Command Window summary
fprintf('\n================ VALIDATION SUMMARY ================\n');

for p = 1:nPanels
    rows = find(PanelName == string(panels(p).name));

    fprintf('\n%s\n', panels(p).name);

    fprintf('  Curve errors [%%]\n');
    fprintf('  %-16s %9s %9s %9s %9s\n', ...
        'Emulator', 'nRMSE_I', 'MAPE_I', 'nRMSE_P', 'MAPE_P');
    for r = rows.'
        fprintf('  %-16s %9.2f %9.2f %9.2f %9.2f\n', ...
            char(Emulator(r)), ...
            nRMSE_I_percent_Isc(r), MAPE_I_percent(r), ...
            nRMSE_P_percent_Pmpp(r), MAPE_P_percent(r));
    end

    fprintf('\n  Characteristic-point errors [%%]\n');
    fprintf('  %-16s %7s %7s %7s %7s %7s\n', ...
        'Emulator', 'Isc', 'Voc', 'Vmpp', 'Impp', 'Pmpp');
    for r = rows.'
        fprintf('  %-16s %7.2f %7.2f %7.2f %7.2f %7.2f\n', ...
            char(Emulator(r)), ...
            Isc_error_percent(r), Voc_error_percent(r), ...
            Vmpp_error_percent(r), Impp_error_percent(r), ...
            Pmpp_error_percent(r));
    end
end

fprintf('\n====================================================\n');

%% 6) Summary plots
summaryLabels = strcat(PanelName, " - ", Emulator);
figS1 = figure('Color','w','Position',[100 100 1300 560]);
bar([nRMSE_I_percent_Isc, MAPE_I_percent]);
grid off; box on;
xticks(1:nRows);
xticklabels(summaryLabels);
xtickangle(35);
ylabel('Current error metrics [%]', 'FontSize', fontAxisLabel);
legend({'$nRMSE_I$', '$MAPE_I$'}, 'Interpreter','latex', 'Location','best', 'FontSize', fontLegend);
set(gca, 'FontSize', 16);
set(findall(gcf, '-property','FontName'), 'FontName', 'Times New Roman');
figS2 = figure('Color','w','Position',[100 100 1300 560]);
bar([nRMSE_P_percent_Pmpp, MAPE_P_percent]);
grid off; box on;
xticks(1:nRows);
xticklabels(summaryLabels);
xtickangle(35);
ylabel('Power error metrics [%]', 'FontSize', fontAxisLabel);
legend({'$nRMSE_P$', '$MAPE_P$'}, 'Interpreter','latex', 'Location','best', 'FontSize', fontLegend);
set(gca, 'FontSize', 16);
set(findall(gcf, '-property','FontName'), 'FontName', 'Times New Roman');
fprintf('\nValidation completed for 3 panels and 3 emulators using common voltage grids.\n');

%% Helper functions

function emu = makeEmu(name, shortName, color, lineStyle, V, I)
    % Store one emulator dataset in a consistent structure.
    emu.name = name;
    emu.shortName = shortName;
    emu.color = color;
    emu.lineStyle = lineStyle;
    emu.V = V(:);
    emu.I = I(:);
end

function [Vout, Iout] = prepareCurve(V, I)
    % Sort by voltage and remove repeated voltage samples.
    Vout = V(:);
    Iout = I(:);
    [Vout, idx] = sort(Vout);
    Iout = Iout(idx);
    [Vout, uniqueIdx] = unique(Vout, 'stable');
    Iout = Iout(uniqueIdx);
end

function [Vaug, Iaug] = augmentCurveWithIscVoc(V, I)
    % Add the estimated (0, Isc) and (Voc, 0) endpoints.
    [V, I] = prepareCurve(V, I);
    Isc = estimateIscFromCurve(V, I);
    Voc = estimateVocFromCurve(V, I);
    Vaug = [0; V(:); Voc];
    Iaug = [Isc; I(:); 0];
    [Vaug, Iaug] = prepareCurve(Vaug, Iaug);
end

function err = relativeErrorPercent(value, reference)
    % Absolute relative error in percent.
    err = 100 * abs(value - reference) / abs(reference);
end

function Isc = estimateIscFromCurve(V, I)
    % Estimate short-circuit current at V = 0.
    [V, I] = prepareCurve(V, I);
    idxZero = find(abs(V) < 1e-9, 1, 'first');
    if ~isempty(idxZero)
        Isc = I(idxZero);
        return;
    end
    Isc = interp1(V, I, 0, 'linear', 'extrap');
end

function Voc = estimateVocFromCurve(V, I)
    % Estimate open-circuit voltage at I = 0.
    [V, I] = prepareCurve(V, I);
    idxZero = find(abs(I) < 1e-9, 1, 'last');
    if ~isempty(idxZero)
        Voc = V(idxZero);
        return;
    end
    IscApprox = max(I);
    lowCurrentRegion = I <= 0.20 * IscApprox;
    if sum(lowCurrentRegion) < 2
        lowCurrentRegion = false(size(I));
        lowCurrentRegion(max(1, numel(I)-3):numel(I)) = true;
    end
    Ilow = I(lowCurrentRegion);
    Vlow = V(lowCurrentRegion);
    [Ilow, idx] = sort(Ilow);
    Vlow = Vlow(idx);
    [Ilow, uniqueIdx] = unique(Ilow, 'stable');
    Vlow = Vlow(uniqueIdx);
    Voc = interp1(Ilow, Vlow, 0, 'linear', 'extrap');
end

function [Vmpp, Impp, Pmpp] = findMPPFromCurve(V, I, nPoints)
    % Estimate the MPP from a dense PCHIP interpolation of the I-V curve.
    if nargin < 3
        nPoints = 3000;
    end
    [V, I] = prepareCurve(V, I);
    Vgrid = linspace(min(V), max(V), nPoints);
    Igrid = interp1(V, I, Vgrid, 'pchip');
    Igrid(Igrid < 0) = 0;
    Pgrid = Vgrid .* Igrid;
    [Pmpp, idxMPP] = max(Pgrid);
    Vmpp = Vgrid(idxMPP);
    Impp = Igrid(idxMPP);
end

function plotPanelIV(panel, Vref, Iref, Isc_ref, Voc_ref, Vmpp_ref, Impp_ref, metrics, ...
    fontAxisLabel, fontTicks, fontLegend, fontCallout, lineWidthEmu, markerSizeRef)
    % Plot the I-V curves and characteristic-point error callouts.
    fig = figure('Color','w','Position',[100 100 1100 650]);
    hold on; grid on; box on;
    nEmu = numel(panel.emu);
    hEmu = gobjects(nEmu,1);
    for k = 1:nEmu
        hEmu(k) = plot(panel.emu(k).V, panel.emu(k).I, panel.emu(k).lineStyle, ...
            'Color', panel.emu(k).color, 'LineWidth', lineWidthEmu);
    end
    hPanel = plot(Vref, Iref, 'ko', ...
        'MarkerSize', markerSizeRef, 'MarkerFaceColor', 'k');
    % Mark characteristic points of the real panel.
    plot(0, Isc_ref, 'ks', 'MarkerSize', 9, ...
        'MarkerFaceColor', 'y', 'HandleVisibility', 'off');
    plot(Vmpp_ref, Impp_ref, 'kd', 'MarkerSize', 9, ...
        'MarkerFaceColor', 'g', 'HandleVisibility', 'off');
    plot(Voc_ref, 0, 'ko', 'MarkerSize', 9, ...
        'MarkerFaceColor', 'c', 'HandleVisibility', 'off');
    % Mark the same characteristic points for each emulator.
    for k = 1:nEmu
        plot(0, metrics(k).Isc_emu, 's', ...
            'MarkerSize', 8, 'MarkerEdgeColor', panel.emu(k).color, ...
            'MarkerFaceColor', panel.emu(k).color, 'HandleVisibility', 'off');
        plot(metrics(k).Vmpp_emu, metrics(k).Impp_emu, 'd', ...
            'MarkerSize', 8, 'MarkerEdgeColor', panel.emu(k).color, ...
            'MarkerFaceColor', panel.emu(k).color, 'HandleVisibility', 'off');
        plot(metrics(k).Voc_emu, 0, 'o', ...
            'MarkerSize', 8, 'MarkerEdgeColor', panel.emu(k).color, ...
            'MarkerFaceColor', panel.emu(k).color, 'HandleVisibility', 'off');
    end
    allV = [Vref(:); Voc_ref; Vmpp_ref];
    allI = [Iref(:); Isc_ref; Impp_ref];
    for k = 1:nEmu
        allV = [allV; panel.emu(k).V(:); metrics(k).Voc_emu; metrics(k).Vmpp_emu]; %#ok<AGROW>
        allI = [allI; panel.emu(k).I(:); metrics(k).Isc_emu; metrics(k).Impp_emu]; %#ok<AGROW>
    end
    xMax = max(allV);
    yMax = max(allI);
    xlim([0, 1.05*xMax]);
    ylim([0, 1.18*yMax]);
    xlabel('$V_{out}\,[\mathrm{V}]$', 'Interpreter','latex', 'FontSize', fontAxisLabel);
    ylabel('$I_{out}\,[\mathrm{A}]$', 'Interpreter','latex', 'FontSize', fontAxisLabel);
    legend([hEmu; hPanel], ...
        {'Open-source PV Emulator - Single Diode Model', ...
         'Open-source PV Emulator - Simplified Exp. Model', ...
         'Commercial PV emulator', ...
         'Real PV panel data'}, ...
        'Location','southwest', 'FontSize', fontLegend, 'Interpreter','none');
    % Build compact callout text for Isc, Voc and MPP errors.
    txtIsc = {'$E_{I_{sc}}$'};
    txtVoc = {'$E_{V_{oc}}$'};
    txtMPP = {'$E_{V_{mpp}}/E_{I_{mpp}}/E_{P_{mpp}}$ [\%]'};
    for k = 1:nEmu
        txtIsc{end+1} = sprintf('%s: %.2f\\%%', metrics(k).shortName, metrics(k).Isc_err); %#ok<AGROW>
        txtVoc{end+1} = sprintf('%s: %.2f\\%%', metrics(k).shortName, metrics(k).Voc_err); %#ok<AGROW>
        txtMPP{end+1} = sprintf('%s: %.2f / %.2f / %.2f\\%%', ...
            metrics(k).shortName, metrics(k).Vmpp_err, metrics(k).Impp_err, metrics(k).Pmpp_err); %#ok<AGROW>
    end
    drawFixedSpeechCallout(gca, ...
        0, Isc_ref, ...
        0.17*xMax, 0.78*yMax, ...
        txtIsc, 'nw', fontCallout);
    drawFixedSpeechCallout(gca, ...
        Vmpp_ref, Impp_ref, ...
        0.58*xMax, 0.67*yMax, ...
        txtMPP, 'ne', fontCallout);
    drawFixedSpeechCallout(gca, ...
        Voc_ref, 0, ...
        0.80*xMax, 0.31*yMax, ...
        txtVoc, 'se', fontCallout);
    set(gca, 'FontSize', fontTicks);
    set(findall(gcf, '-property','FontName'), 'FontName', 'Times New Roman');
    set(findall(gca, '-property','FontName'), 'FontName', 'Times New Roman');
end

function drawFixedSpeechCallout(ax, xTip, yTip, xCenter, yTop, str, tailCorner, fontSz)
    % Draw a fixed gray callout box with a triangular tail.
    faceCol = [0.94 0.94 0.94];
    edgeCol = [0.45 0.45 0.45];
    textCol = [0.10 0.10 0.10];
    lineW  = 0.9;
    axes(ax); %#ok<LAXES>
    xLims = xlim(ax);
    yLims = ylim(ax);
    xRange = diff(xLims);
    yRange = diff(yLims);
    % Measure text size so the box fits the content.
    htmp = text(ax, xCenter, yTop, str, ...
        'FontSize', fontSz, ...
        'Interpreter', 'latex', ...
        'HorizontalAlignment', 'center', ...
        'VerticalAlignment', 'middle', ...
        'BackgroundColor', 'none', ...
        'EdgeColor', 'none', ...
        'Color', textCol, ...
        'Visible', 'off', ...
        'Units', 'data', ...
        'HandleVisibility', 'off');
    drawnow;
    ext = get(htmp, 'Extent');
    delete(htmp);
    textW = ext(3);
    textH = ext(4);
    xPad = 0.09 * textW;
    yPad = 0.15 * textH;
    boxW = textW + 2*xPad;
    boxH = textH + 2*yPad;
    % Keep the callout inside the axes.
    xCenter = min(max(xCenter, xLims(1) + 0.025*xRange + boxW/2), ...
                  xLims(2) - 0.025*xRange - boxW/2);
    yTop = min(max(yTop, yLims(1) + 0.06*yRange + boxH), ...
               yLims(2) - 0.035*yRange);
    xl = xCenter - boxW/2;
    xr = xCenter + boxW/2;
    yt = yTop;
    yb = yTop - boxH;
    yCenter = 0.5*(yb + yt);
    triW = 0.10 * boxW;
    triH = 0.18 * boxH;
    switch lower(tailCorner)
        case 'nw'
            p1 = [xl,        yt - triH];
            p2 = [xl + triW, yt];
        case 'ne'
            p1 = [xr - triW, yt];
            p2 = [xr,        yt - triH];
        case 'sw'
            p1 = [xl,        yb + triH];
            p2 = [xl + triW, yb];
        case 'se'
            p1 = [xr - triW, yb];
            p2 = [xr,        yb + triH];
        otherwise
            error('tailCorner must be ''nw'', ''ne'', ''sw'' or ''se''.');
    end
    patch(ax, [p1(1), p2(1), xTip], [p1(2), p2(2), yTip], faceCol, ...
        'EdgeColor', edgeCol, ...
        'LineWidth', lineW, ...
        'HandleVisibility', 'off');
    patch(ax, [xl, xr, xr, xl], [yb, yb, yt, yt], faceCol, ...
        'EdgeColor', edgeCol, ...
        'LineWidth', lineW, ...
        'HandleVisibility', 'off');
    text(ax, xCenter, yCenter, str, ...
        'FontSize', fontSz, ...
        'Interpreter', 'latex', ...
        'HorizontalAlignment', 'center', ...
        'VerticalAlignment', 'middle', ...
        'BackgroundColor', 'none', ...
        'EdgeColor', 'none', ...
        'Color', textCol, ...
        'HandleVisibility', 'off');
end

