clc; clear; close all; format long

%% Module Parameters

ns = 72;       % number of series cells

Vmp_mod_ref = 34.0;     % voltage at maximum power point (V)
Imp_mod_ref = 4.4;     % current at maximum power point (A)
Voc_mod_ref = 43.4;     % open-circuit voltage (V)
Isc_mod_ref = 4.8;     % short-circuit current (A)

Tref = 25 + 273.15; % Reference temperature (K)
Gref = 1000; % Reference irradiance (W/m²)

alpha   = 0.0003;  % temperature coefficient of Isc (1/K or 1/ºC)
beta = -0.0037; % temperature coefficient of Voc (1/K or 1/ºC)

%% Operating Conditions

T = 25 + 273.15; % Current temperature (K)
G = 1000; % Current irradiance (W/m²)

%% Physical Constants

q = 1.60217662e-19; % Elementary charge (C)
k = 1.38064852e-23; % Boltzmann constant (J/K)
E_G0 = 1.166;            % Band gap energy at 0K (eV)
k1   = 4.73e-4;          % Coefficient k1 (eV/K)
k2   = 636;              % Coefficient k2 (K)

%% Parameter Estimation via fsolve

Rs_0 = (Voc_mod_ref - Vmp_mod_ref)/Imp_mod_ref;
Rp_0 = (Vmp_mod_ref)/(Isc_mod_ref - Imp_mod_ref);

x0 = [Isc_mod_ref; log10(1e-9); 0.5*ns; log10(Rs_0); Rp_0];
opts = optimoptions('fsolve', ...
    'Display', 'off', 'TolFun', 1e-10, 'TolX', 1e-10, 'MaxIter', 1000, 'MaxFunctionEvaluations', 2000);
fun = @(x) residuals_2_20(x, Voc_mod_ref, Isc_mod_ref, Vmp_mod_ref, Imp_mod_ref, q, k, Tref);
[xsol, ~, exitflag] = fsolve(fun, x0, opts);

if exitflag <= 0
    warning('fsolve did not converge (exitflag = %d)', exitflag);
end

Iph_ref = xsol(1);
Is0_ref = xsol(2);
A = xsol(3);
Rs = xsol(4);
Rp = xsol(5);

fprintf('Iph_ref = %.6f A\n', Iph_ref);
fprintf('Is0_ref = %.2e A\n', Is0_ref);
fprintf('A       = %.6f\n', A);
fprintf('Rs      = %.6f Ohms\n', Rs);
fprintf('Rp      = %.6f Ohms\n', Rp);

%% Characteristic Points (Adaptive Embedded Strategy - High Precision Knee)

Voc_estimation = Voc_mod_ref * (1 + beta * (T - Tref));

Vmp_estimation = Vmp_mod_ref * (Voc_estimation / Voc_mod_ref);

factors = [
    0.00;  
    0.20;   
    0.40;   
    0.60;   
    0.70;   
    0.80;   
    0.85; 
    0.88; 
    0.90; 
    0.92; 
    0.93; 
    0.94; 
    0.95; 
    0.96; 
    0.97; 
    0.98; 
    0.99; 
    1.00;   
    1.01; 
    1.02; 
    1.03; 
    1.04;
    1.05;   
    1.08;
    1.12;   
];
V_mod = factors * Vmp_estimation;

V_mod = V_mod(V_mod < Voc_estimation); 
V_mod = V_mod(V_mod >= 0);

V_mod = [V_mod; Voc_estimation];

V_mod = unique(V_mod);

I_mod = arrayfun(@(V) solve_I_V_2_11(V, Iph_ref, Is0_ref, A, Rs, Rp, ...
    q, k, G, Gref, alpha, T, Tref, E_G0, k1, k2, ns), V_mod);

points_V = V_mod;
points_I = I_mod;

%% Piecewise PV Model

% Coefficient initialization
coeffs = zeros(length(points_V) - 1, 2);
for i = 1:length(points_V) - 1
    a = (points_I(i + 1) - points_I(i)) / (points_V(i + 1) - points_V(i));
    b = points_I(i) - a * points_V(i);
    coeffs(i, :) = [a, b]; % Store [a, b] in each row
end

% MATLAB function for the PV model (accepts a voltage vector)
modele_pv = @(V) piecewise_pv_model(V, points_V, coeffs);

function I_out = piecewise_pv_model(V_in, points_V, coeffs)
    I_out = zeros(size(V_in)); % Initialize the output current array
    for k = 1:numel(V_in) % Iterate over each voltage input
        V = V_in(k);
        found = false;
        for i = 1:length(points_V) - 1
            % Use a small tolerance for floating-point comparisons
            tol = 1e-9;
            if (points_V(i) - tol <= V && V <= points_V(i + 1) + tol)
                a = coeffs(i, 1);
                b = coeffs(i, 2);
                I_out(k) = a * V + b;
                found = true;
                break;
            end
        end
        if ~found
            I_out(k) = 0; % Default to 0 if V is outside defined segments
        end
    end
end

%% Provided Data: (R, V_test, I_test, V*, I*)
mesures = [
    
    % ---------------------CS6P-250P---------------------
    
    % --- Sref = 1000, Tref = 25, S = 765, T = 44.5 ---
    % 36.470619	10.97875	22.642405	6.816042	3.321928621;
    % 36.469147	10.580585	23.414888	6.793228	3.446798488;
    % 36.011593	10.135595	24.056635	6.770828	3.552982737;
    % 36.008194	9.758461	24.833183	6.729959	3.689945659;
    % 36.011219	9.268014	25.833637	6.648664	3.885538057;
    % 36.012985	8.819846	26.70089	6.539245	4.083176269;
    % 36.00943	8.417071	27.416254	6.408448	4.278142539;
    % 36.009499	8.022567	28.057871	6.251021	4.488526114;
    % 36.011665	7.394004	28.952494	5.944597	4.870388018;
    % 36.014778	6.823575	29.632019	5.614259	5.277992875;
    % 36.012032	6.198761	30.29425	5.214557	5.809553908;
    % 36.01263	5.744209	30.687252	4.894783	6.26937946;
    % 36.012844	5.41966	30.963587	4.659785	6.644853142;
    % 36.013351	4.885607	31.42931	4.263731	7.371316342;
    % 36.014511	4.044705	31.895897	3.582153	8.904113532;
    % 36.019165	2.940577	32.502819	2.653505	12.24901366;
    % 36.026188	1.518282	33.319305	1.404204	23.72825102;
    % 36.04068	0.271503	34.06929	0.256652	132.7450789;

    % ---------------------Kyocera KB260-6BPA---------------------
    
    % --- Sref = 1000, Tref = 25, S = 735, T = 33.6 ---
    % 38.970139	10.476988	24.569763	6.605496	3.719593956;
    % 38.970062	11.051579	23.384464	6.631636	3.526198362;
    % 38.9701	10.006746	25.593412	6.571879	3.894382718;
    % 38.478722	9.397497	26.668896	6.513233	4.09457116;
    % 38.481834	8.949709	27.661718	6.433278	4.299785895;
    % 38.481712	8.526176	28.55562	6.326908	4.513361029;
    % 38.97081	8.168219	29.45533	6.17379	4.771028817;
    % 38.484024	7.527627	30.391579	5.944713	5.11237111;
    % 38.485474	6.842545	31.405592	5.583774	5.624438238;
    % 38.480385	6.376169	31.972301	5.297785	6.03503181;
    % 38.480625	6.013403	32.411503	5.064976	6.399142464;
    % 38.485367	5.482672	32.913673	4.688922	7.019454152;
    % 38.481773	5.186515	33.195103	4.473986	7.419581331;
    % 38.482933	4.712598	33.651192	4.120906	8.165969328;
    % 38.483788	3.500327	34.301086	3.119885	10.99434306;
    % 38.970074	1.690633	35.331486	1.532781	23.0505767;
    % 38.513306	0.255853	36.170616	0.24029	150.5290108;


    % ---------------- KC85TS - charge R ----------------- 
    
    % --- Sref = 1000, Tref = 25, S = 400, T = 25 ---
    % 22.491886	11.192251	4.290837	2.135176	2.009594057;
    % 22.785145	4.797872	10.134062	2.133931	4.749011097;
    % 22.49477	3.785123	12.669167	2.1318	5.94294352;
    % 22.494724	3.463998	13.826797	2.12921	6.493862512;
    % 22.785276	3.220089	15.001449	2.120053	7.075978289;
    % 22.492458	3.039437	15.612347	2.109718	7.40020562;
    % 22.494221	2.910696	16.172581	2.092691	7.728126608;
    % 22.784504	2.841723	16.606222	2.071157	8.017847995;
    % 22.493919	2.67594	17.095009	2.03367	8.405989664;
    % 22.494114	2.549509	17.511776	1.984805	8.822920136;
    % 22.494652	2.356512	18.038937	1.889737	9.545739434;
    % 22.495331	2.168992	18.432518	1.777257	10.37132953;
    % 22.497145	1.950675	18.831503	1.632836	11.53300331;
    % 22.500511	1.719844	19.143658	1.46326	13.08288206;
    % 22.504717	1.484131	19.473244	1.284213	15.16356243;
    % 22.512615	0.991903	19.819204	0.873231	22.6964045;
    % 22.518124	0.453164	20.202765	0.406569	49.69086428;
    % 22.51823	0.317759	20.30147	0.286478	70.86572093;

  
    % --------------- KC85TS - charge R+L --------------- 
    
    % % ---- Sref = 1000, Tref = 25, S = 400, T = 25 ----
    % 22.491726	11.401453	4.212106	2.135191	1.972706891;
    % 22.784645	4.763543	10.206792	2.133915	4.7831296;
    % 22.494654	3.790882	12.650108	2.131843	5.933883499;
    % 22.494001	3.46787	13.811143	2.129245	6.486403866;
    % 22.494892	3.18538	14.974236	2.12042	7.061919808;
    % 22.493298	3.035744	15.629344	2.10937	7.409484348;
    % 22.49297	2.908536	16.181227	2.092373	7.733433284;
    % 22.495255	2.81023	16.588015	2.072266	8.004771106;
    % 22.493977	2.670704	17.113245	2.031851	8.422490133;
    % 22.493221	2.52665	17.581903	1.974965	8.902387131;
    % 22.494623	2.340781	18.077917	1.881181	9.609876455;
    % 22.493725	2.189714	18.393583	1.790574	10.27245062;
    % 22.785147	1.989317	18.813585	1.64257	11.45374931;
    % 22.501472	1.710468	19.156639	1.456207	13.15516201;
    % 22.504978	1.49112	19.463346	1.28959	15.09266201;
    % 22.512096	1.000689	19.813057	0.880713	22.49661013;
    % 22.518711	0.457511	20.199619	0.410394	49.22006413;
    % 22.517979	0.306336	20.309837	0.276297	73.50726573;


    % ---------------------KC200GT---------------------
    
    % --- Sref = 1000, Tref = 25, S = 511, T = 54,3 ---
    % 30.799538	10.383743	12.533917	4.225679	2.966130887;
    % 30.407265	8.72176	14.705232	4.217923	3.486368054;
    % 30.411446	7.384953	17.256174	4.190397	4.118028435;
    % 30.412125	7.033544	18.053396	4.175287	4.323869473;
    % 30.408638	6.741723	18.770275	4.161449	4.510514246;
    % 30.799961	6.556941	19.421268	4.134554	4.69730665;
    % 30.411461	6.154503	20.231714	4.094382	4.941335225;
    % 30.410652	5.776842	21.174362	4.022306	5.264234496;
    % 30.410862	5.474993	21.89076	3.941084	5.554502264;
    % 30.411789	5.199258	22.50289	3.847137	5.849256213;
    % 30.411108	4.899352	23.111118	3.723294	6.20716978;
    % 30.409916	4.590909	23.679567	3.574845	6.623942297;
    % 30.410358	4.117438	24.440063	3.309085	7.385746513;
    % 30.410532	3.747369	24.939966	3.073253	8.11516852;
    % 30.412104	3.481879	25.292833	2.895774	8.734394673
    % 30.413876	3.004618	25.808134	2.549612	10.12237705;
    % 30.411911	2.620044	26.232281	2.259961	11.60740429;
    % 30.799515	1.95355	26.693003	1.693082	15.76592451;
    % 30.428387	1.320891	27.100122	1.176412	23.03625091;
    % 30.434181	0.844936	27.427097	0.761451	36.01951669;
    % 30.799423	0.298323	27.814808	0.269414	103.2418805;

        
    % --------------------- ME Solar MESM-50W ---------------------
    
    % --- Sref = 1000, Tref = 25, S = 1000, T = 25 ---
    % 23.415052	11.967575	5.922083	3.026813	1.956540758;
    % 23.112965	5.560591	12.551384	3.019652	4.156566386;
    % 23.116062	5.030003	13.841137	3.0118	4.595636164;
    % 23.115747	4.57807	15.125275	2.995558	5.049234567;
    % 23.115692	4.439414	15.540461	2.984576	5.2069242;
    % 23.114132	4.291178	15.992638	2.96906	5.386431396;
    % 23.116524	4.171403	16.356771	2.951598	5.541666243;
    % 23.114021	4.05022	16.715565	2.929032	5.706856395;
    % 23.115362	3.948944	17.007088	2.905429	5.853554845;
    % 23.114918	3.798939	17.418226	2.862687	6.084572292;
    % 23.114498	3.65326	17.789234	2.811599	6.327087896;
    % 23.114141	3.503831	18.138689	2.749611	6.596820059;
    % 23.114452	3.3509	18.464882	2.676852	6.8979839;
    % 23.11492	3.15309	18.834526	2.569205	7.330877061;
    % 23.116348	2.903193	19.264559	2.419445	7.962387655;
    % 23.115911	2.630473	19.65144	2.236233	8.787742601;
    % 23.116068	2.396255	19.971132	2.070245	9.64674809;
    % 23.123291	1.748199	20.552086	1.553807	13.22692329;
    % 23.136288	0.80102	21.464041	0.743124	28.8835255;
    % 23.140648	0.309762	21.969181	0.29408	74.70477761;

    % ---------------------Renogy RNG-50DB-H – 50 W---------------------

    % --- Sref = 1000, Tref = 25, S = 1000, T = 25 ---
    % 23.424229	11.58681	5.894531	2.915734	2.021628516;
    % 23.72953	5.693411	12.125871	2.909352	4.167894088;
    % 23.426764	4.922009	13.816704	2.902917	4.759593195;
    % 23.424603	4.170353	16.145973	2.874516	5.616936208;
    % 23.426395	4.002863	16.715784	2.856222	5.852410632;
    % 23.423584	3.848856	17.231913	2.831469	6.08585614;
    % 23.424973	3.715216	17.664825	2.801653	6.305143785;
    % 23.426899	3.545583	18.180428	2.751547	6.607347794;
    % 23.427967	3.344581	18.727219	2.673501	7.004754814;
    % 23.427265	3.14929	19.188766	2.579515	7.4389046;
    % 23.427248	2.939382	19.602093	2.459445	7.970128627;
    % 23.426811	2.739423	19.970377	2.335244	8.551730355;
    % 23.426422	2.418296	20.420042	2.107949	9.687161312;
    % 23.429234	2.159837	20.759869	1.91376	10.84768675;
    % 23.437695	1.558239	21.242069	1.412265	15.04113534;
    % 23.447426	0.686716	21.980999	0.643768	34.14428645;
    % 23.450655	0.307881	22.318256	0.293014	76.16788276;

   
    % ---------------------Shell Solar SQ150-PC---------------------

    % --- Sref = 1000, Tref = 25, S = 1000, T = 25 ---
    45.001869	9.165687	23.428795	4.771824	4.909819599;
    45.002346	8.243192	25.957676	4.754732	5.459335247;
    45.00293	7.68873	27.723593	4.736563	5.8531034;
    45.002434	7.429503	28.600252	4.721648	6.057260516;
    45.006149	7.207458	29.371901	4.70373	6.244384988;
    45.009201	6.976861	30.186039	4.679128	6.451210354;
    45.00576	6.623515	31.42337	4.624589	6.794845985;
    45.007038	6.370541	32.285755	4.569901	7.064869677;
    45.006256	6.092615	33.184868	4.492322	7.387019007;
    45.003723	5.867533	33.870052	4.415938	7.669956417;
    45.006004	5.55717	34.74757	4.290498	8.098726535;
    45.006401	5.199488	35.65852	4.119548	8.655930214;
    45.007046	4.942261	36.230217	3.97847	9.106570365;
    45.009201	4.642657	36.874702	3.803591	9.694707449;
    45.006428	4.315445	37.490471	3.594777	10.42915068;
    45.005337	3.995132	38.099937	3.382138	11.26504507;
    45.006454	3.076956	39.200233	2.680002	14.6269417;
    45.010605	1.953921	40.635674	1.764004	23.03604414;
    45.032486	0.160229	43.159355	0.153565	281.0494253;
    
    ];

if ~isempty(mesures)
    V_violet = mesures(:, 1); % Purple Points (Test Point)
    I_violet = mesures(:, 2); % Purple Points (Test Point)
    V_noir = mesures(:, 3);   % Black Points (Intersection)
    I_noir = mesures(:, 4);   % Black Points (Intersection)
    R_data = mesures(:, 5);   % Black dotted-line
else
    V_violet = []; I_violet = []; V_noir = []; I_noir = []; R_data = [];
end

%% Plotting (Figure Generation)
figure; 
hold on; 

% 1. Plot "Load Lines" FIRST
if isempty(R_data)
    Rs_to_plot = [0.1, 0.5, 1:1:15, 20:5:50, 60:10:200, 300, 400, 500, 1000];
    Rs_to_plot = sort([Rs_to_plot, inf]);
else
    Rs_to_plot = R_data; 
    Rs_to_plot = sort(unique(Rs_to_plot));
end

V_line = linspace(0, max(V_violet) * 1.25, 100);
h_load = []; % Inicializa o handle para a legenda

for i = 1:length(Rs_to_plot)
    r_val = Rs_to_plot(i);
    
    % Definimos HandleVisibility como 'on' apenas para a primeira linha para não poluir a legenda
    if i == 1
        vis = 'on';
    else
        vis = 'off';
    end

    if isinf(r_val) 
        h_temp = plot(V_line, zeros(size(V_line)), 'k--', 'LineWidth', 0.8, 'Color', [0.5 0.5 0.5], 'HandleVisibility', vis);
    elseif r_val == 0 
        h_temp = plot(zeros(size(V_line)), linspace(0, max(I_violet) * 1.1, 100), 'k--', 'LineWidth', 0.8, 'Color', [0.5 0.5 0.5], 'HandleVisibility', vis);
    else
        I_line = V_line ./ r_val;
        h_temp = plot(V_line, I_line, 'k--', 'LineWidth', 0.8, 'Color', [0.5 0.5 0.5], 'HandleVisibility', vis);
    end
    
    if i == 1, h_load = h_temp; end % Guarda o handle da primeira linha
end

% 2. Plot the I-V Model (P-V Curve)
V_plot = linspace(0, max(points_V), 500); 
I_plot = modele_pv(V_plot);
h_model = plot(V_plot, I_plot, 'r-', 'LineWidth', 2, 'DisplayName', 'I-V model');

% 3. Plot "Test Points" (Purple Points)
h_test_point = []; 
if ~isempty(V_violet)
    h_test_point = scatter(V_violet, I_violet, 40, 'm', 'filled', 'DisplayName', 'Test point');
end

% 4. Plot "Intersections" (Black Points)
h_intersection = []; 
if ~isempty(V_noir) && ~isempty(I_noir)
    h_intersection = scatter(V_noir, I_noir, 70, 'ks', 'filled', 'DisplayName', 'Intersection');
end

%% Plot Configurations

xlabel('$V_{out}^{OT}\,[\mathrm{V}]$', 'Interpreter','latex','FontSize',30);
ylabel('$I_{out}^{OT}\,[\mathrm{A}]$', 'Interpreter','latex','FontSize',30);
xlim([0 1.1 *Voc_estimation]);
ylim([0 1.1 * max(I_violet)]);
set(gca, 'FontSize', 30); 
grid off;

legend_handles = [h_model];
legend_labels = {'Single-diode model I–V curve'};

if ~isempty(h_load)
    legend_handles = [legend_handles, h_load];
    legend_labels = [legend_labels, 'Load line'];
end

if ~isempty(h_test_point) && ishandle(h_test_point)
    legend_handles = [legend_handles, h_test_point];
    legend_labels = [legend_labels, 'Test point'];
end

if ~isempty(h_intersection) && ishandle(h_intersection)
    legend_handles = [legend_handles, h_intersection];
    legend_labels = [legend_labels, 'Emulation Points'];
end

leg = legend(legend_handles, legend_labels, 'Location', 'NorthWest', 'FontSize', 20);
set(leg, 'Box', 'on'); 
set(findall(gcf, '-property', 'FontName'), 'FontName', 'Times New Roman');
set(findall(gca, '-property', 'FontName'), 'FontName', 'Times New Roman');
hold off;

%% --- Auxiliary Functions ---

function F = residuals_2_20(x, Voc, Isc, Vmp, Imp, q, k, Tref)
    Iph = x(1); Is0 = x(2); A = x(3); Rs = x(4); Rp = x(5);
    C = q / (A * k * Tref);
    
    F = zeros(5,1);
    F(1) = Iph - Is0*(exp(C*Isc*Rs)-1) - (Isc*Rs)/Rp - Isc;
    F(2) = Iph - Is0*(exp(C*Voc)-1) - Voc/Rp;
    F(3) = Iph - Is0*(exp(C*(Vmp + Imp*Rs))-1) - (Vmp + Imp*Rs)/Rp - Imp;
    F(4) = Iph - 2*Vmp/Rp - Is0*((1 + C*(Vmp - Imp*Rs))*exp(C*(Vmp + Imp*Rs)) - 1);
    F(5) = Rs + Is0 * C * Rp * (Rs - Rp) * exp(C * Isc * Rs);
end

function I = solve_I_V_2_11(V, Iph_ref, Is0_ref, A, Rs, Rp, ...
    q, k, S, Sref, alpha, T, Tref, E_G0, k1, k2, ns) 
    
    % 1. Photocurrent (Iph) 
    Iph = Iph_ref * (S / Sref) * (1 + alpha * (T - Tref));
    
    % 2. Calculation of Saturation Current (Is)
    
    % 2a. Calculation of band gap energy in eV at temperature T: Eg(T)
    Eg_T_eV = E_G0 - (k1 * T^2) / (T + k2);
    
    % 2b. Calculation of the temperature difference term
    temp_diff = (1 / Tref) - (1 / T);
    
    % 2c. Calculation of the exponent 
    exponent_term = (ns*q / (A * k)) * Eg_T_eV * temp_diff; 
    
    % 2d. Final calculation of Is(T)
    Is = Is0_ref * (T / Tref)^3 * exp(exponent_term);
    
    % 3. Iterative Solver (Newton-Raphson)
    I = Iph; % Initial guess
    for it = 1:30
        V_diode = V + I * Rs;
        arg_exp = min(q * V_diode / (A * k * T), 700); 
        expo = exp(arg_exp);
        
        f = Iph - Is * (expo - 1) - V_diode / Rp - I;
        df = -Is * expo * (q * Rs / (A * k * T)) - Rs / Rp - 1;
        
        dI = -f / df;
        I = I + dI;
        
        if abs(dI) < 1e-6, break;end
    end
end
