clear; clc; close all;

% Validates the photovoltaic emulator (PVE) against real PV-panel data.
% It computes relative errors at Isc, Voc and the MPP, then plots the I-V curves.
% The script does not calculate RMSE/MAPE metrics and does not generate a CSV file.

%% 1) Input data
data = struct([]);

%% CanadianSolar CS6P-250P
data(1).name = 'CanadianSolar CS6P-250P';
data(1).condition = 'G = 765 W/m^2, T = 44.5 degC';
data(1).V_emulator = [
    0;
    22.642405;
    23.414888;
    24.056635;
    24.833183;
    25.833637;
    26.700890;
    27.416254;
    28.057871;
    28.952494;
    29.632019;
    30.294250;
    30.687252;
    30.963587;
    31.429310;
    31.895897;
    32.502819;
    33.319305;
    34.069290;
    34.25
];
data(1).I_emulator = [
    7;
    6.816042;
    6.793228;
    6.770828;
    6.729959;
    6.648664;
    6.539245;
    6.408448;
    6.251021;
    5.944597;
    5.614259;
    5.214557;
    4.894783;
    4.659785;
    4.263731;
    3.582153;
    2.653505;
    1.404204;
    0.256652;
    0
];
data(1).V_panel = [
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
data(1).I_panel = [
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

%% Kyocera KB260-6BPA
data(2).name = 'Kyocera KB260-6BPA';
data(2).condition = 'G = 735 W/m^2, T = 33.6 degC';
data(2).V_emulator = [
    0;
    23.384464;
    24.569763;
    25.593412;
    26.668896;
    27.661718;
    28.555620;
    29.455330;
    30.391579;
    31.405592;
    31.972301;
    32.411503;
    32.913673;
    33.195103;
    33.651192;
    34.301086;
    35.331486;
    36.170616;
    36.3
];
data(2).I_emulator = [
    6.75;
    6.631636;
    6.605496;
    6.571879;
    6.513233;
    6.433278;
    6.326908;
    6.173790;
    5.944713;
    5.583774;
    5.297785;
    5.064976;
    4.688922;
    4.473986;
    4.120906;
    3.119885;
    1.532781;
    0.240290;
    0
];
data(2).V_panel = [
    0.0455969407525445;
    2.9271959778498773;
    4.902704170259113;
    6.88347582118527;
    9.038080417500383;
    11.176899184353912;
    13.510629760851112;
    15.860159814501406;
    18.0095873908246;
    20.332779674072395;
    22.6717646106255;
    24.857993191432502;
    26.98107843870138;
    29.141271766676834;
    30.932938702702867;
    32.3931145335085;
    33.45817143875639;
    34.191355106534374;
    34.68732138929362;
    35.25719902557649;
    35.78478474359848;
    36.30708425066834
];
data(2).I_panel = [
    6.727119597608984;
    6.729434019490265;
    6.725222853395946;
    6.725329760352628;
    6.736240521267449;
    6.7298792755918075;
    6.727846338739699;
    6.730132042305772;
    6.6546867392647755;
    6.654812127477331;
    6.6657328411244094;
    6.661533046823786;
    6.588245218439663;
    6.294752137400725;
    5.778873016808391;
    4.88948723260063;
    4.017351278763285;
    3.1344029398705855;
    2.4004055294681965;
    1.5152894827191918;
    0.7704993045386512;
    0.04297999950084552
];

%% Kyocera KC200GT
data(3).name = 'Kyocera KC200GT';
data(3).condition = 'G = 511 W/m^2, T = 54.3 degC';
data(3).V_emulator = [
    0;
    12.533917;
    14.705232;
    17.256174;
    18.053396;
    18.770275;
    19.421268;
    20.231714;
    21.174362;
    21.890760;
    22.502890;
    23.111118;
    23.679567;
    24.440063;
    24.939966;
    25.292833;
    25.808134;
    26.232281;
    26.693003;
    27.100122;
    27.427097;
    27.814808;
    28
];
data(3).I_emulator = [
    4.23;
    4.225679;
    4.217923;
    4.190397;
    4.175287;
    4.161449;
    4.134554;
    4.094382;
    4.022306;
    3.941084;
    3.847137;
    3.723294;
    3.574845;
    3.309085;
    3.073253;
    2.895774;
    2.549612;
    2.259961;
    1.693082;
    1.176412;
    0.761451;
    0.269414;
    0
];
data(3).V_panel = [
    0.06633956327079193;
    2.7458921458044943;
    4.750187578542191;
    6.738397879641979;
    8.55510640498975;
    10.398634093999334;
    12.247397477176047;
    14.085538251089517;
    15.92895053464986;
    17.767131103542148;
    19.471307382584936;
    21.31996331611866;
    22.868301124243413;
    24.28771513306349;
    25.433614258948253;
    26.2901123451705;
    26.830422064051454;
    27.41886779290256;
    27.846553737021324;
    28.24767918166411
];
data(3).I_panel = [
    4.14356789750856;
    4.1571047845181806;
    4.161687801527659;
    4.161814349519464;
    4.155246813078904;
    4.162047320620259;
    4.099788754649053;
    4.090994861739706;
    4.0331914003872384;
    4.046674736070509;
    4.044555482254624;
    3.922148395763573;
    3.6816528875172967;
    3.2696145187726584;
    2.7439448794896295;
    2.2115736529084318;
    1.6769545779212427;
    1.0888732202257323;
    0.5074647987806244;
    0.05749032916413199
];

%% Kyocera KC85TS - R load
data(4).name = 'Kyocera KC85TS - R load';
data(4).condition = 'STC';
data(4).V_emulator = [
    0;
    4.290837;
    10.134062;
    12.669167;
    13.826797;
    15.001449;
    15.612347;
    16.172581;
    16.606222;
    17.095009;
    17.511776;
    18.038937;
    18.432518;
    18.831503;
    19.143658;
    19.473244;
    19.819204;
    20.202765;
    20.301470;
    20.55
];
data(4).I_emulator = [
    2.14;
    2.135176;
    2.133931;
    2.131800;
    2.129210;
    2.120053;
    2.109718;
    2.092691;
    2.071157;
    2.033670;
    1.984805;
    1.889737;
    1.777257;
    1.632836;
    1.463260;
    1.284213;
    0.873231;
    0.406569;
    0.286478;
    0
];
data(4).V_panel = [
    20.981002130722764;
    0.06608783557910758;
    2.087423154611276;
    3.7422899997017183;
    5.503561760222222;
    7.335766446749726;
    9.179774931381964;
    10.834674829311304;
    12.678672292808942;
    14.39267378341276;
    15.598401507942457;
    16.49682293222942;
    17.17673915663834;
    17.732773312588492;
    18.17067586419277;
    18.620630027412773;
    19.04114353518865;
    19.485584047785522;
    19.894706778785416;
    20.286231201766142;
    20.54744023267662
];
data(4).I_panel = [
    0.020684459626193252;
    2.1848822023067904;
    2.1729592301689484;
    2.168795023167493;
    2.160178040675044;
    2.149337974661873;
    2.1418545624403293;
    2.130981453493248;
    2.1257343416974335;
    2.114872247505312;
    2.103915432139303;
    2.088428741983925;
    2.038238640744611;
    1.9410630969490947;
    1.829329572532087;
    1.6706359366326193;
    1.4895737887520175;
    1.2503722307386989;
    0.9820921567196246;
    0.6869731713502247;
    0.4510919206320487
];

%% Kyocera KC85TS - R+L load
data(5).name = 'Kyocera KC85TS - R+L load';
data(5).condition = 'STC';
data(5).V_emulator = [
    0;
    4.212106;
    10.206792;
    12.650108;
    13.811143;
    14.974236;
    15.629344;
    16.181227;
    16.588015;
    17.113245;
    17.581903;
    18.077917;
    18.393583;
    18.813585;
    19.156639;
    19.463346;
    19.813057;
    20.199619;
    20.309837;
    20.55
];
data(5).I_emulator = [
    2.14;
    2.135191;
    2.133915;
    2.131843;
    2.129245;
    2.120420;
    2.109370;
    2.092373;
    2.072266;
    2.031851;
    1.974965;
    1.881181;
    1.790574;
    1.642570;
    1.456207;
    1.289590;
    0.880713;
    0.410394;
    0.276297;
    0
];
data(5).V_panel = data(4).V_panel;
data(5).I_panel = data(4).I_panel;

%% ME Solar MESM-50W
data(6).name = 'ME Solar MESM-50W';
data(6).condition = 'G = 1000 W/m^2, T = 25 degC';
data(6).V_emulator = [
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
];
data(6).I_emulator = [
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
];
data(6).V_panel = [
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
data(6).I_panel = [
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

%% Renogy RNG-50DB-H
data(7).name = 'Renogy RNG-50DB-H';
data(7).condition = 'G = 1000 W/m^2, T = 25 degC';
data(7).V_emulator = [
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
];
data(7).I_emulator = [
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
];
data(7).V_panel = [
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
data(7).I_panel = [
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

%% Shell Solar SQ150-PC
data(8).name = 'Shell Solar SQ150-PC';
data(8).condition = 'G = 1000 W/m^2, T = 25 degC';
data(8).V_emulator = [
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
];
data(8).I_emulator = [
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
];
data(8).V_panel = [
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
data(8).I_panel = [
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

%% 2) Calculate characteristic-point errors and plot I-V curves
nCases = numel(data);
for c = 1:nCases
    Vemu = data(c).V_emulator(:);
    Iemu = data(c).I_emulator(:);
    Vref = data(c).V_panel(:);
    Iref = data(c).I_panel(:);
    % Prepare the real-panel data for interpolation.
    [Vref, idx] = sort(Vref);
    Iref = Iref(idx);
    [Vref, uniqueIdx] = unique(Vref, 'stable');
    Iref = Iref(uniqueIdx);
    % Interpolate the real-panel curve only to obtain a smooth reference line.
    Vinterp_plot = linspace(min(Vref), max(Vref), 500);
    Iinterp_plot = interp1(Vref, Iref, Vinterp_plot, 'pchip');
    % Estimate Isc, Voc and MPP for both the real panel and the PVE.
    Isc_ref_char = estimateIscFromCurve(Vref, Iref);
    Isc_pve_char = estimateIscFromCurve(Vemu, Iemu);
    Voc_ref_char = estimateVocFromCurve(Vref, Iref);
    Voc_pve_char = estimateVocFromCurve(Vemu, Iemu);
    [Vmpp_ref_char, Impp_ref_char, Pmpp_ref_char] = ...
        findMPPFromCurve(Vref, Iref, 3000);
    [Vmpp_pve_char, Impp_pve_char, Pmpp_pve_char] = ...
        findMPPFromCurve(Vemu, Iemu, 3000);
    % Relative error = |PVE - reference| / |reference| * 100.
    Isc_error = 100 * abs(Isc_pve_char - Isc_ref_char) / abs(Isc_ref_char);
    Voc_error = 100 * abs(Voc_pve_char - Voc_ref_char) / abs(Voc_ref_char);
    VMPP_error = 100 * abs(Vmpp_pve_char - Vmpp_ref_char) / abs(Vmpp_ref_char);
    IMPP_error = 100 * abs(Impp_pve_char - Impp_ref_char) / abs(Impp_ref_char);
    MPP_error = 100 * abs(Pmpp_pve_char - Pmpp_ref_char) / abs(Pmpp_ref_char);
    % Print the five characteristic-point errors in the Command Window.
    fprintf('\n%s\n', data(c).name);
    fprintf('  E_Isc  = %.3f %%\n', Isc_error);
    fprintf('  E_Voc  = %.3f %%\n', Voc_error);
    fprintf('  E_Vmpp = %.3f %%\n', VMPP_error);
    fprintf('  E_Impp = %.3f %%\n', IMPP_error);
    fprintf('  E_Pmpp = %.3f %%\n', MPP_error);
    % Plot the PVE curve, measured reference points and interpolated reference curve.
    fig1 = figure('Color','w','Position',[100 100 1000 560]);
    hold on; grid off; box on;
    plot(Vemu, Iemu, '-b', 'LineWidth', 2);
    plot(Vref, Iref, 'ko', 'MarkerSize', 6, 'MarkerFaceColor', 'k');
    plot(Vinterp_plot, Iinterp_plot, '--r', 'LineWidth', 1.5);
    % Mark the reference and PVE characteristic points.
    plot(0, Isc_ref_char, 'ks', 'MarkerSize', 9, ...
        'MarkerFaceColor', 'y', 'HandleVisibility','off');
    plot(Vmpp_ref_char, Impp_ref_char, 'kd', 'MarkerSize', 9, ...
        'MarkerFaceColor', 'g', 'HandleVisibility','off');
    plot(Voc_ref_char, 0, 'ko', 'MarkerSize', 9, ...
        'MarkerFaceColor', 'c', 'HandleVisibility','off');
    plot(0, Isc_pve_char, 'bs', 'MarkerSize', 8, ...
        'MarkerFaceColor', 'b', 'HandleVisibility','off');
    plot(Vmpp_pve_char, Impp_pve_char, 'bd', 'MarkerSize', 8, ...
        'MarkerFaceColor', 'b', 'HandleVisibility','off');
    plot(Voc_pve_char, 0, 'bo', 'MarkerSize', 8, ...
        'MarkerFaceColor', 'b', 'HandleVisibility','off');
    xlabel('$V_{out}$ [V]', 'Interpreter','latex', 'FontSize', 24);
    ylabel('$I_{out}$ [A]', 'Interpreter','latex', 'FontSize', 24);
    legend({'Open-source PVE', ...
        'Real PV panel data', ...
        'Interpolated real-panel reference'}, ...
        'Location','best');
    xMax = max([Vref; Vemu; Voc_ref_char; Voc_pve_char]);
    yMax = max([Iref; Iemu; Isc_ref_char; Isc_pve_char]);
    xlim([0, 1.05*xMax]);
    ylim([0, 1.12*yMax]);
    % Text shown in the three callout boxes.
    txtIsc = { ...
        sprintf('$E_{I_{sc}}=%.2f\\,\\%%$', Isc_error), ...
        sprintf('$I_{sc,\\mathrm{realPV}}=%.2f\\,\\mathrm{A}$', Isc_ref_char), ...
        sprintf('$I_{sc,\\mathrm{PVE}}=%.2f\\,\\mathrm{A}$', Isc_pve_char)};
    txtMPP = { ...
        sprintf('$E_{V_{mpp}}=%.2f\\,\\%%$', VMPP_error), ...
        sprintf('$E_{I_{mpp}}=%.2f\\,\\%%$', IMPP_error), ...
        sprintf('$E_{P_{mpp}}=%.2f\\,\\%%$', MPP_error), ...
        sprintf('$V_{mpp,\\mathrm{realPV}}=%.2f\\,\\mathrm{V}$', Vmpp_ref_char), ...
        sprintf('$I_{mpp,\\mathrm{realPV}}=%.2f\\,\\mathrm{A}$', Impp_ref_char)};
    txtVoc = { ...
        sprintf('$E_{V_{oc}}=%.2f\\,\\%%$', Voc_error), ...
        sprintf('$V_{oc,\\mathrm{realPV}}=%.2f\\,\\mathrm{V}$', Voc_ref_char), ...
        sprintf('$V_{oc,\\mathrm{PVE}}=%.2f\\,\\mathrm{V}$', Voc_pve_char)};
    Vcurve_all = [Vemu(:); Vinterp_plot(:); Vref(:)];
    Icurve_all = [Iemu(:); Iinterp_plot(:); Iref(:)];
    % Add callouts near Isc, MPP and Voc.
    drawSpeechCallout(gca, ...
        0, Isc_ref_char, ...
        0.13*xMax, 0.95*yMax, ...
        txtIsc, 'nw', Vcurve_all, Icurve_all);
    drawSpeechCallout(gca, ...
        Vmpp_ref_char, Impp_ref_char, ...
        max(0.16*xMax, Vmpp_ref_char - 0.10*xMax), ...
        min(0.78*yMax, Impp_ref_char + 0.05*yMax), ...
        txtMPP, 'ne', Vcurve_all, Icurve_all);
    drawSpeechCallout(gca, ...
        Voc_ref_char, 0, ...
        max(0.70*xMax, Voc_ref_char - 0.08*xMax), ...
        0.20*yMax, ...
        txtVoc, 'se', Vcurve_all, Icurve_all);
    set(gca, 'FontSize', 20);
    set(findall(gcf, '-property','FontName'), 'FontName', 'Times New Roman');
end

fprintf('\nCharacteristic-point validation completed.\n');

%% 3) Helper functions
function Isc = estimateIscFromCurve(V, I)
    % Return I at V = 0; extrapolate when the curve has no exact V = 0 point.
    V = V(:);
    I = I(:);
    [V, idx] = sort(V);
    I = I(idx);
    [V, uniqueIdx] = unique(V, 'stable');
    I = I(uniqueIdx);
    idxZero = find(abs(V) < 1e-9, 1, 'first');
    if ~isempty(idxZero)
        Isc = I(idxZero);
        return;
    end
    Isc = interp1(V, I, 0, 'linear', 'extrap');
end

function Voc = estimateVocFromCurve(V, I)
    % Return V at I = 0; extrapolate from the low-current region when needed.
    V = V(:);
    I = I(:);
    [V, idx] = sort(V);
    I = I(idx);
    [V, uniqueIdx] = unique(V, 'stable');
    I = I(uniqueIdx);
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
    [Ilow, idx2] = sort(Ilow);
    Vlow = Vlow(idx2);
    [Ilow, uniqueIdx2] = unique(Ilow, 'stable');
    Vlow = Vlow(uniqueIdx2);
    Voc = interp1(Ilow, Vlow, 0, 'linear', 'extrap');
end

function [Vmpp, Impp, Pmpp] = findMPPFromCurve(V, I, nPoints)
    % Interpolate the I-V curve on a dense grid and select the maximum-power point.
    if nargin < 3
        nPoints = 3000;
    end
    V = V(:);
    I = I(:);
    [V, idx] = sort(V);
    I = I(idx);
    [V, uniqueIdx] = unique(V, 'stable');
    I = I(uniqueIdx);
    Vgrid = linspace(min(V), max(V), nPoints);
    Igrid = interp1(V, I, Vgrid, 'pchip');
    Igrid(Igrid < 0) = 0;
    Pgrid = Vgrid .* Igrid;
    [Pmpp, idxMPP] = max(Pgrid);
    Vmpp = Vgrid(idxMPP);
    Impp = Igrid(idxMPP);
end

function drawSpeechCallout(ax, xTip, yTip, xCenterPref, yTopPref, str, tailCorner, curveX, curveY)
    % Draw a callout box and move it automatically to reduce overlap with the curve.
    if nargin < 7 || isempty(tailCorner)
        tailCorner = 'se';
    end
    if nargin < 8
        curveX = [];
        curveY = [];
    end
    % Callout appearance.
    faceCol = [0.94 0.94 0.94];
    edgeCol = [0.45 0.45 0.45];
    textCol = [0.10 0.10 0.10];
    lineW  = 0.9;
    fontSz = 13;
    axes(ax); %#ok<LAXES>
    xLims  = xlim(ax);
    yLims  = ylim(ax);
    xRange = diff(xLims);
    yRange = diff(yLims);
    % Measure the text to size the callout box.
    htmp = text(ax, xCenterPref, yTopPref, str, ...
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
    ext0 = get(htmp, 'Extent');   % [x_left, y_bottom, width, height]
    delete(htmp);
    textW = ext0(3);
    textH = ext0(4);
    xPad = 0.14 * textW;
    yPad = 0.22 * textH;
    boxW = textW + 2*xPad;
    boxH = textH + 2*yPad;
    % Test nearby box positions and keep the one with the least curve overlap.
    offsets = [ ...
         0.00,  0.00; ...
         0.00,  0.08; ...
         0.00, -0.08; ...
        -0.10,  0.00; ...
         0.10,  0.00; ...
        -0.10,  0.08; ...
         0.10,  0.08; ...
        -0.10, -0.08; ...
         0.10, -0.08; ...
        -0.16,  0.00; ...
         0.16,  0.00];
    bestScore = inf;
    bestGeom = struct();
    for k = 1:size(offsets,1)
        xCenter = xCenterPref + offsets(k,1)*xRange;
        yTop    = yTopPref    + offsets(k,2)*yRange;
        xCenter = min(max(xCenter, xLims(1) + 0.04*xRange + boxW/2), ...
                      xLims(2) - 0.04*xRange - boxW/2);
        yTop = min(max(yTop, yLims(1) + 0.08*yRange + boxH), ...
                   yLims(2) - 0.035*yRange);
        xl = xCenter - boxW/2;
        xr = xCenter + boxW/2;
        yt = yTop;
        yb = yTop - boxH;
        overlapScore = 0;
        if ~isempty(curveX)
            inBox = curveX >= xl & curveX <= xr & curveY >= yb & curveY <= yt;
            overlapScore = sum(inBox);
        end
        dPref = hypot((xCenter - xCenterPref)/xRange, (yTop - yTopPref)/yRange);
        dTip  = hypot((xCenter - xTip)/xRange, (yTop - yTip)/yRange);
        tooClosePenalty = 0;
        if dTip < 0.08
            tooClosePenalty = 100;
        end
        score = 1000*overlapScore + 10*dPref + tooClosePenalty;
        if score < bestScore
            bestScore = score;
            bestGeom.xCenter = xCenter;
            bestGeom.xl = xl;
            bestGeom.xr = xr;
            bestGeom.yb = yb;
            bestGeom.yt = yt;
        end
    end
    xCenter = bestGeom.xCenter;
    xl = bestGeom.xl;
    xr = bestGeom.xr;
    yb = bestGeom.yb;
    yt = bestGeom.yt;
    yCenter = 0.5*(yb + yt);
    triW = 0.12 * boxW;
    triH = 0.20 * boxH;
    % Select which corner of the box connects to the characteristic point.
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
    % Draw the tail, box and centered text.
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
