clc; clear; close all; format long

%% Module Parameters

ns = 72;       % number of series cells

Vmp_mod_ref = 17.4;     % voltage at maximum power point (V)
Imp_mod_ref = 5.02;     % current at maximum power point (A)
Voc_mod_ref = 21.7;     % open-circuit voltage (V)
Isc_mod_ref = 5.34;     % short-circuit current (A)

Tref = 25 + 273.15; % Reference temperature (K)
Gref = 1000; % Reference irradiance (W/m²)

alpha   = 0.0004;  % temperature coefficient of Isc (1/K or 1/ºC)
beta = -0.38; % temperature coefficient of Voc (1/K or 1/ºC)

%% Operating Conditions

T = 25 + 273.15; % Current temperature (K)
G = 400; % Current irradiance (W/m²)

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
    % 35.784859	12.165794	20.089266	6.829757	2.94143203;
    % 35.786652	11.396839	21.406616	6.817284	3.140050495;
    % 35.784744	10.764368	22.612379	6.802004	3.324370142;
    % 35.78397	10.236521	23.688374	6.776401	3.495716089;
    % 35.786598	9.671537	24.878904	6.723669	3.700197615;
    % 35.785908	9.176556	25.916357	6.645714	3.899709948;
    % 35.783024	8.804281	26.656096	6.558634	4.064275579;
    % 35.787361	8.393621	27.416744	6.430365	4.26363729;
    % 35.783806	7.966219	28.124405	6.261077	4.491943638;
    % 35.786224	7.636677	28.616514	6.106682	4.686098605;
    % 35.786362	7.083742	29.345022	5.808709	5.051900861;
    % 35.785271	6.442036	30.058086	5.411033	5.554962611;
    % 36.241776	5.862995	30.656363	4.959418	6.181443669;
    % 35.78756	5.102286	31.281769	4.459889	7.014024116;
    % 35.788845	4.317395	31.720699	3.826633	8.289454202;
    % 35.790173	2.6143	32.716629	2.389792	13.69015755;
    % 35.804443	1.283223	33.539886	1.202061	27.90198334;
    % 35.815018	0.287514	34.182873	0.274412	124.5677048;

    % ---------------------Kyocera KB260-6BPA---------------------
    
    % --- Sref = 1000, Tref = 25, S = 735, T = 33.6 ---
    % 38.478035	10.264172	24.649357	6.575316	3.748771466;
    % 38.478245	9.852852	25.582468	6.550722	3.905289829;
    % 38.480354	9.395959	26.636724	6.504035	4.09541523;
    % 38.482475	8.908169	27.760311	6.426134	4.319908517;
    % 38.480614	8.570052	28.505875	6.348569	4.490126043;
    % 38.481716	8.149579	29.370399	6.220003	4.721926822;
    % 38.479599	7.588441	30.389585	5.993034	5.070818053;
    % 38.48497	7.885333	29.870871	6.120357	4.880576574;
    % 38.486076	7.263933	30.910389	5.834084	5.298242021;
    % 38.486244	6.853367	31.497862	5.608923	5.615670245;
    % 38.48027	6.395314	32.049358	5.326514	6.016948045;
    % 38.480522	6.130786	32.373322	5.157776	6.276604878;
    % 38.482475	5.848843	32.652031	4.96269	6.579502447;
    % 38.482391	5.435788	33.04081	4.667143	7.079450962;
    % 38.482742	4.888101	33.570869	4.264192	7.87273861;
    % 38.969883	3.59669	34.342472	3.169607	10.83493064;
    % 38.486595	2.187879	35.130131	1.997071	17.59082727;
    % 38.511208	0.261213	36.306252	0.246257	147.4323654;


    % ---------------- KC85TS - charge R ----------------- 
    
    % --- Sref = 1000, Tref = 25, S = 400, T = 25 ---
    % 21.917355	12.392923	3.771689	2.13266	1.768537413;
    % 21.916889	3.961433	11.757546	2.125153	5.532564479;
    % 21.917103	3.453911	13.46219	2.121503	6.345590838;
    % 21.635138	3.04005	15.029916	2.11192	7.116707072;
    % 21.635174	2.889083	15.72958	2.100471	7.488596605;
    % 21.635191	2.762043	16.312876	2.082573	7.833039226;
    % 21.635189	2.648232	16.803375	2.056799	8.169672875;
    % 21.63496	2.591394	17.031826	2.04004	8.348770612;
    % 21.636127	2.528483	17.26795	2.018001	8.556958099;
    % 21.636707	2.436547	17.583126	1.980066	8.880070664;
    % 21.637428	2.318445	17.934851	1.921715	9.332731961;
    % 21.635872	2.168547	18.304523	1.834649	9.977125325;
    % 21.640064	1.933421	18.780472	1.677932	11.19262998;
    % 21.643158	1.67558	19.152412	1.48275	12.91681807;
    % 21.647425	1.443203	19.492191	1.299517	14.99956599;
    % 21.917015	0.969186	19.843899	0.877512	22.61382067;
    % 21.658646	0.473057	20.207394	0.441359	45.78448383;
    % 21.662161	0.321711	20.323677	0.301833	67.33417817;


    % --------------- KC85TS - charge R+L --------------- 
    
    % ---- Sref = 1000, Tref = 25, S = 400, T = 25 ----
    22.491859	8.839152	5.422981	2.131196	2.544571687;
    22.785051	4.69319	10.325287	2.126768	4.854919295;
    22.491726	3.834816	12.457455	2.123983	5.865138751;
    22.494526	3.482104	13.70118	2.120913	6.460038672;
    22.493938	3.150763	15.073689	2.111396	7.139205057;
    22.496943	3.016575	15.675177	2.101856	7.457778744;
    22.492147	2.923689	16.085083	2.090854	7.693068478;
    22.493511	2.841831	16.439856	2.077012	7.915147337;
    22.493788	2.71131	16.966587	2.045084	8.296278784;
    22.492807	2.604371	17.350945	2.009011	8.636560477;
    22.494133	2.4528	17.817669	1.94287	9.170798355;
    22.494324	2.267365	18.28035	1.842608	9.920911013;
    22.494551	2.06133	18.67819	1.711611	10.91263728;
    22.500011	1.816658	19.048414	1.537975	12.38538598;
    22.502754	1.546906	19.429449	1.335638	14.54694236;
    22.785137	1.093687	19.7838	0.949623	20.83332017;
    22.514206	0.725954	20.036785	0.646071	31.01328647;
    22.517738	0.315571	20.33769	0.285019	71.35555875;


    % ---------------------KC200GT---------------------
    
    % --- Sref = 1000, Tref = 25, S = 511, T = 54,3 ---
    % 30.409756	12.985798	9.868175	4.213981	2.341770169
    % 30.408054	9.382547	13.614193	4.200723	3.240916623;
    % 30.79991	7.846811	16.422968	4.184036	3.925149784;
    % 30.412991	7.081144	17.876022	4.162126	4.294925718;
    % 30.409344	6.461602	19.423361	4.127219	4.706161946;
    % 30.411236	6.173776	20.185144	4.09778	4.925873034;
    % 30.411804	5.800905	21.167236	4.037549	5.242595446;
    % 30.411945	5.398189	22.170599	3.935331	5.633731699;
    % 30.411777	5.028286	23.001207	3.803022	6.048139348;
    % 30.413933	4.639174	23.76317	3.624704	6.555892564;
    % 30.411734	4.145412	24.581303	3.350668	7.33623952;
    % 30.412529	3.61688	25.310093	3.010062	8.408495573;
    % 30.412239	3.119902	25.86165	2.653071	9.747816775;
    % 30.414455	2.73973	26.282421	2.367517	11.10125967;
    % 30.413723	2.225109	26.626436	1.948026	13.66841921;
    % 30.427082	1.440454	27.169167	1.28622	21.12326585;
    % 30.434351	0.807039	27.623262	0.732496	37.71114382;
    % 30.437847	0.320416	27.982399	0.294568	94.99470071;
    
    % --------------------- ME Solar MESM-50W ---------------------
    
    % --- Sref = 1000, Tref = 25, S = 1000, T = 25 ---
    % 22.2323	12.606372	5.335613	3.025451	1.763576075;
    % 22.522175	6.00836	11.316692	3.019014	3.748472846;
    % 22.522625	5.191639	13.071528	3.013088	4.338249663;
    % 22.23243	4.632921	14.41359	3.003587	4.798792244;
    % 22.233393	4.410621	15.08972	2.993472	5.040875612;
    % 22.233898	4.184878	15.801759	2.974217	5.312913953;
    % 22.234447	3.957017	16.524055	2.940751	5.618991543;
    % 22.234411	3.798619	17.004604	2.905137	5.853288158;
    % 22.232758	3.570288	17.643137	2.833255	6.227161692;
    % 22.234035	3.413467	18.036154	2.76899	6.513621934;
    % 22.234131	3.221753	18.466652	2.675841	6.901251607;
    % 22.234028	2.954685	18.971693	2.521153	7.525006614;
    % 22.233418	2.640574	19.479301	2.313479	8.41991693;
    % 22.23427	2.323803	19.947374	2.08479	9.56804954;
    % 22.239458	1.824127	20.410786	1.674135	12.19183997;
    % 22.253723	0.893098	21.333824	0.856181	24.91742283;
    % 22.262161	0.32815	21.935129	0.32333	67.84130455;


    % ---------------------Renogy RNG-50DB-H – 50 W---------------------

    % --- Sref = 1000, Tref = 25, S = 1000, T = 25 ---
    % 22.825958	12.842241	5.179413	2.914019	1.777412227;
    % 22.825714	5.812697	11.411573	2.906022	3.926870822;
    % 22.53413	4.284508	15.187758	2.887712	5.259443463;
    % 22.533112	4.006109	16.149439	2.871171	5.624687279;
    % 22.533306	3.873092	16.624014	2.857385	5.817911832;
    % 22.532803	3.694992	17.250757	2.828827	6.098201481;
    % 22.534595	3.617272	17.517166	2.81187	6.229721146;
    % 22.531389	3.49121	17.929047	2.778083	6.453747782;
    % 22.532784	3.313553	18.462984	2.71507	6.800187104;
    % 22.532482	3.218438	18.722519	2.67424	7.00106161;
    % 22.533628	3.074749	19.081297	2.603673	7.328607317;
    % 22.533056	2.874862	19.507757	2.488882	7.837959775;
    % 22.535362	2.673573	19.885263	2.359168	8.428930453;
    % 22.533543	2.343695	20.378786	2.11958	9.614539673;
    % 22.535172	2.148954	20.665606	1.970672	10.48657818;
    % 22.548851	1.268318	21.433184	1.205564	17.77855344;
    % 22.826054	0.470768	22.157703	0.456983	48.48693059;
   
    % ---------------------Shell Solar SQ150-PC---------------------

    % --- Sref = 1000, Tref = 25, S = 1000, T = 25 ---
    % 45.006153	9.481608	22.64146	4.769958	4.746679111;
    % 45.006439	8.779803	24.412127	4.762289	5.126133042;
    % 45.006214	8.201163	26.068161	4.750216	5.487784345;
    % 45.006348	7.660948	27.799574	4.732023	5.874775757;
    % 45.005699	7.415041	28.636826	4.718141	6.069514667;
    % 45.004837	7.178912	29.45788	4.698951	6.269033237;
    % 45.005585	6.936991	30.314299	4.672532	6.487767018;
    % 45.007683	6.64842	31.332647	4.628379	6.769680486;
    % 45.003468	6.374206	32.268707	4.570478	7.060247747;
    % 45.00351	6.078027	33.229164	4.487823	7.404294688;
    % 45.006283	5.854822	33.908855	4.411169	7.687045089;
    % 45.008648	5.596772	34.641445	4.307623	8.041893406;
    % 45.00668	5.223126	35.597893	4.131216	8.616807497;
    % 45.004463	4.908363	36.30035	3.959058	9.168936146;
    % 45.004627	4.754821	36.647137	3.871837	9.465051602;
    % 45.01001	4.395983	37.328506	3.645755	10.23889592;
    % 45.005409	4.020102	38.051414	3.398937	11.19509247;
    % 45.008011	3.115784	39.146133	2.709982	14.44516347;
    % 45.569969	1.69742	41.002605	1.527292	26.84660497;
    % 45.033173	0.202371	43.096004	0.193666	222.5274648;

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
xlabel('Voltage (V)', 'FontSize', 30); 
ylabel('Current(A)', 'FontSize', 30);
xlim([0 1.05*Voc_estimation]);
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
    C   = q/(A*k*Tref);
    cap = 700;
    a_sc = min(C*Rs*Isc, cap);
    a_oc = min(C*Voc, cap);
    a_mp = min(C*(Rs*Imp+Vmp), cap);
    a_eq5 = min(C*Isc, cap); 
    F = zeros(5,1);
    F(1) = Isc - ( Iph - Is0*(exp(a_sc)-1) - (Rs*Isc)/Rp );
    F(2) = 0   - ( Iph - Is0*(exp(a_oc)-1) - (Voc)/Rp );
    F(3) = Imp - ( Iph - Is0*(exp(a_mp)-1) - (Rs*Imp+Vmp)/Rp );
    F(4) = Rs + (q*Is0*Rp*(Rs-Rp)/(A*k*Tref))*exp(a_eq5);
    term_coeff = 1 + q*(Vmp - Rs*Imp)/(A*k*Tref);
    F(5) = Iph - 2*Vmp/Rp + Is0 - Is0 * term_coeff * exp(a_mp);
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
