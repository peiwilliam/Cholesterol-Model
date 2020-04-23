tic
close all;
%% Constants
MWc = 386.65; %Molecular weight of cholesterol, [g/mol]
MWbs = 408.57; %Molecular weight of bile salt (cholic acid) [g/mol]

%% ODE Solver Conditions
options = odeset('NonNegative',1); %Ensure non-negative concentrations do not exist

%% Loop Tracker %% UNUSED
global count;

%% Initialize Dependent Variables In Regulatory Network:
%(1)C_FCG, (2)C_CMG, (3)C_BSG, (4)C_CMB, (5)C_VLDLB, (6)C_LDLB, (7)C_HDLeB, 
%(8)HDLfB, (9)C_CML, (10)C_VLDLL, (11)C_LDLL, (12)C_HDLeL, (13)C_HDLfL,  
%(14)C_FCL, (15)C_BSL, (16)C_LDLP, (17)C_HDLeP, (18)C_HDLfP, (19)C_FCP

run Define_Constants_Diet %Initialize constants

global V_liver;
global V_peripheral;
global V_gut;

C = zeros(20,1);
C(1) = 3150/MWc/V_gut; %C_FCG, [mM] (Source: A whole-body mathematical model of cholesterolmetabolism and its age-associated dysregulation)
C(2) = 7.12889485549129; %C_CMG, [mM] Obtained from fitting
C(3) = 467/MWbs/V_gut; %C_BSG [mM] (Source: A whole-body mathematical model of cholesterolmetabolism and its age-associated dysregulation)
C(4) = 0.76; %C_CMB, (Source: High throughput prediction of chylomicron triglycerides in human plasma by nuclear magnetic resonance and chemometrics)[mM]
C(5) = 26/MWc*10; %C_VLDLB, [mM] (Source: Hypercholesterolemic Effect of Dietary Cholesterol in Diets Enriched in Polyunsaturated and Saturated Fat, normal VLDL level is 0.1 to 1.7 mM)
C(6) = 4.48; %C_LDLB, [mM] (Source: Assessment of the longer-term effects of a dietary portfolio of cholesterol-lowering foods in hypercholesterolemia)
C(7) = 54/MWc; %C_HDLeB, [mM] (Source: An apolipoprotein A-I mimetic dose-dependently increases the formation of preb1 HDL in human plasma)
C(8) = 1.24; %HDLfB, [mM] (Source: Assessment of the longer-term effects of a dietary portfolio of cholesterol-lowering foods in hypercholesterolemia, normal HDL level is 2.2 to 3.3 mM)
C(9) = 1.49084621815203; %C_CML, [mM] Obtained form fitting
C(10) = 4.17906517693227; %C_VLDLL, [mM] Obtained from fitting
C(11) = 0.683016370321173; %C_LDLL [mM] Obtained from fitting
C(12) = 3.20620527934023; %C_HDLeL [mM] Obtained from fitting
C(13) = 0.266292711361245; %C_HDLfL [mM] Obtained from fitting
% C(14) = 60000/MWc/V_liver; %C_FCL, [mM]
C(14) = 8; %C_FCL, [mM] (Source: A physiologically based in silico kinetic model predicting plasma cholesterol concentrations in humans)
C(15) = 400/MWbs/V_liver; %C_BSL [mM] (Source: A whole-body mathematical model of cholesterolmetabolism and its age-associated dysregulation)
C(16) = 0.00413095730723442; %C_LDLP [mM] obtained from fitting
C(17) = 0.150886276369151; %C_HDLeP {mM] obtained from fitting
C(18) = 0.00907520469090437; %C_HDLfP [mM] obtained from fitting
C(19) = 57516/MWc/V_peripheral; %C_FCP, [mM] (Source: A whole-body mathematical model of cholesterolmetabolism and its age-associated dysregulation)
C(20) = 0.0015; %C_HMGCoA, [mM] % unused

%% Model initialization UNUSED

global Model;
%Set the FBA execution requirement
global FBA_Size;
%Daily Dosage Intake Model

%% Model inputs
global Q_diet; %Diet cholesterol, [mmol/day]
global Q_syn_liver; % Liver cholesterol synthesis [mmol/day]

%% Declare input and output

% Clear all local x, y variables
clear x;
clear y;
% Declare global variables x, y
global x;
global y;
% Clear the cell arrays
x = [];
y = [];

% Declare result set's column index tracker
col = 1; %Initialize to 1

%% David Jenkins Data

tLDL = [0; 2; 4; 8; 12; 24; 32; 42; 52];
tLDL = tLDL*7;
tHDL = [0; 2; 4; 8; 12; 18; 24; 42; 52];
tHDL = tHDL*7;

concsLDL = [4.48; 3.6736; 3.6288; 3.7184; 3.7632; 3.808; 3.8528; 3.8976; 3.87];
LDLoutlier = 4.1664;
tLDLoutlier = 18*7;
concsHDL = [1.24; 1.2338; 1.2524; 1.2462; 1.271; 1.302; 1.2896; 1.3268; 1.28];
HDLoutlier = 1.395;
tHDLoutlier = 32*7;

%% Run baseline model

Q_diet = 94.99/MWc; % [mmol/day] Baseline (Source: Assessment of the longer-term effects of a dietary portfolio of cholesterol =lowering foods in hypercholesterolemia) taken as an average of all the diets
Q_syn_liver =  29.8*24/MWc;

[x{1,col}, y{1,col}] = ode15s(@Cholesterol_Model, [0 365], C, options);

titles = ["FCG" "CMG" "BSG" "CMB" "VLDLB" "LDLB" "HDLeB" "HDLfB" "CML" "VLDLL" "LDLL" "HDLeL" "HDLfL" "FCL" "BSL" "LDLP" "HDLeP" "HDLfP" "FCP" "HMGCoA"];

figure(1)
for i = 1:(size(y{1,col}, 2) - 1)
    subplot(5, 4, i);
    hold on
    plot(x{1,col}, y{1,col}(:, i));
    if i == 6
        plot(tLDL, concsLDL, 'x');
        plot(tLDLoutlier, LDLoutlier, 'rx');
    end
    if i == 8
        plot(tHDL, concsHDL, 'x');
        plot(tHDLoutlier, HDLoutlier, 'rx');
    end
    if i == 3
        ylim([0.92 0.925]);
    end
    ylabel('Concentration (mM)');
    xlabel('Time (days)');
    title(titles(i));
end

figure(2)
hold on
plot(x{1,col}, y{1,col}(:, 6), 'LineWidth', 2);
plot(tLDL, concsLDL, 'kx', 'MarkerSize', 10, 'LineWidth', 2);
plot(tLDLoutlier, LDLoutlier, 'rx', 'MarkerSize', 10, 'LineWidth', 2);
ylabel('LDL-C Concentration (mM)');
xlabel('Time (Days)');
title('Model versus Clinical Data');
legend('Model Fit', "Jenkins et al.");
set(gca,'FontSize',15);

figure(3)
hold on
plot(x{1,col}, y{1,col}(:, 8), 'LineWidth', 2);
plot(tHDL, concsHDL, 'kx', 'MarkerSize', 10, 'LineWidth', 2);
plot(tHDLoutlier, HDLoutlier, 'rx', 'MarkerSize', 10, 'LineWidth', 2);
ylabel('HDL-C Concentration (mM)');
xlabel('Time (Days)');
title('Model versus Clinical Data');
legend('Model Fit', "Jenkins et al.");
set(gca,'FontSize',15);

disp("Section 1 Done")

col = col + 1;

%% Sensitivity analysis wrt. changes in different variables

% Liver cholesterol synthesis rate

liverSyn = [29.8*24/MWc*1.5 29.8*24/MWc*1.25 29.8*24/MWc 29.8*24/MWc*0.75 29.8*24/MWc*0.5];

figure(4)
for i = 1:length(liverSyn)
    Q_syn_liver = liverSyn(i);
    [x{1,col}, y{1,col}] = ode15s(@Cholesterol_Model, [0 1095], C, options); %Column 1
    subplot(2, 1, 1)
    hold on
    plot(x{1,col}, y{1,col}(:, 6), 'LineWidth', 2);
    subplot(2, 1, 2)
    hold on
    plot(x{1,col}, y{1,col}(:, 8), 'LineWidth', 2);
end

subplot(2, 1, 1)
hold on
ylabel('LDL-C Concentration (mM)');
xlabel('Time (Days)');
title('Sensitivity of LDL-C to Changes in Liver Cholesterol Synthesis');
legend('1.5x synthesis rate', '1.25x synthesis rate', 'Baseline synthesis rate', '0.75x synthesis rate', '0.5x synthesis rate', 'Location', 'bestoutside');
set(gca,'FontSize',15);

subplot(2, 1, 2)
hold on
ylabel('HDL-C Concentration (mM)');
xlabel('Time (Days)');
title('Sensitivity of HDL-C to Changes in Liver Cholesterol Synthesis');
legend('1.5x synthesis rate', '1.25x synthesis rate', 'Baseline synthesis rate', '0.75x synthesis rate', '0.5x synthesis rate', 'Location', 'bestoutside');
set(gca,'FontSize',15);

col = col + 1;

% LDL-C transport rate into liver

Q_syn_liver = 29.8*24/MWc;

LDLBL = [0.036*1.5 0.036*1.25 0.036 0.036*0.75 0.036*0.5];
global K_LDLBL;

figure(5)
for i = 1:length(LDLBL)
    K_LDLBL = LDLBL(i);
    [x{1,col}, y{1,col}] = ode15s(@Cholesterol_Model, [0 1095], C, options); %Column 1
    subplot(2, 1, 1)
    hold on
    plot(x{1,col}, y{1,col}(:, 6), 'LineWidth', 2);
    subplot(2, 1, 2)
    hold on
    plot(x{1,col}, y{1,col}(:, 8), 'LineWidth', 2);
end

subplot(2, 1, 1)
hold on
ylabel('LDL-C Concentration (mM)');
xlabel('Time (Days)');
title('Sensitivity of LDL-C to Changes in LDL-C Transport from the Blood into the Liver');
legend('1.5x baseline transport rate', '1.25x baseline transport rate', 'Baseline transport rate', '0.75x baseline transport rate', '0.5x baseline transport rate', 'Location', 'bestoutside');
set(gca,'FontSize',15);

subplot(2, 1, 2)
hold on
ylabel('HDL-C Concentration (mM)');
xlabel('Time (Days)');
title('Sensitivity of HDL-C to Changes in LDL-C Transport from the Blood into the Liver');
legend('1.5x baseline transport rate', '1.25x baseline transport rate', 'Baseline transport rate', '0.75x baseline transport rate', '0.5x baseline transport rate', 'Location', 'bestoutside');
set(gca,'FontSize',15);

col = col +1;

% Synthesis rate of VLDL in the liver from free cholesterol

run Define_Constants_Diet;

FCVLDL = [0.115*1.5 0.115*1.25 0.115 0.115*0.75 0.115*0.5];
global K_FCVLDL;

figure(6)
for i = 1:length(FCVLDL)
    K_FCVLDL = FCVLDL(i);
    [x{1,col}, y{1,col}] = ode15s(@Cholesterol_Model, [0 1095], C, options); %Column 1
    subplot(2, 1, 1)
    hold on
    plot(x{1,col}, y{1,col}(:, 6), 'LineWidth', 2);
    subplot(2, 1, 2)
    hold on
    plot(x{1,col}, y{1,col}(:, 8), 'LineWidth', 2);
end

subplot(2, 1, 1)
hold on
ylabel('LDL-C Concentration (mM)');
xlabel('Time (Days)');
title('Sensitivity of LDL-C to Changes in Liver VLDL Synthesis Rate');
legend('1.5x VLDL synthesis rate', '1.25x VLDL synthesis rate', 'Baseline VLDL synthesis rate', '0.75x VLDL synthesis rate', '0.5x VLDL synthesis rate', 'Location', 'bestoutside');
set(gca,'FontSize',15);

subplot(2, 1, 2)
hold on
ylabel('HDL-C Concentration (mM)');
xlabel('Time (Days)');
title('Sensitivity of HDL-C to Changes in Liver VLDL Synthesis Rate');
legend('1.5x VLDL synthesis rate', '1.25x VLDL synthesis rate', 'Baseline VLDL synthesis rate', '0.75x VLDL synthesis rate', '0.5x VLDL synthesis rate', 'Location', 'bestoutside');
set(gca,'FontSize',15);

disp("Section 2 Done")

col = col +1;

%% Run multiple types of diets

run Define_Constants_Diet

Diets = [1.0685*1000/MWc 0.4764*1000/MWc 94.99/MWc 0.0403*1000/MWc 0.0014*1000/MWc]; %348 taken from a different source not VMH
ratio = [];
ssLDL = [];
ssHDL = [];

Q_syn_liver =  29.8*24/MWc;
% Q_syn_peripheral = 441/MWc; %Baseline
% Q_syn_intestine = 49/MWc; %Baseline

figure(7)
for i = 1:length(Diets)
    Q_diet = Diets(i);
    [x{1,col}, y{1,col}] = ode15s(@Cholesterol_Model, [0 365], C, options); %Column 1
    subplot(2, 1, 1)
    hold on
    plot(x{1,col}, y{1,col}(:, 6), 'LineWidth', 2);
    subplot(2, 1, 2)
    hold on
    plot(x{1,col}, y{1,col}(:, 8), 'LineWidth', 2);
    ratio(i) = y{1, 2}(end, 6)/y{1, 2}(end, 8);
    ssLDL(i) = y{1, 2}(end, 6);
    ssHDL(i) = y{1, 2}(end, 8);
end

subplot(2, 1, 1)
ylabel('LDL-C Concentration (mM)');
xlabel('Time (Days)');
title('Comparison of LDL-C Concentration of Different Diets');
legend('Unhealthy', 'Average', "Jenkins et al.", 'Vegetarian', 'Vegan', 'Location', 'bestoutside');
set(gca,'FontSize',15);

subplot(2, 1, 2)
ylabel('HDL-C Concentration (mM)');
xlabel('Time (Days)');
title('Comparison of HDL-C Concentration of Different Diets');
legend('Unhealthy', 'Average', "Jenkins et al.", 'Vegetarian', 'Vegan', 'Location', 'bestoutside');
set(gca,'FontSize',15);

disp("Section 3 Done")

col = col +1;

disp(toc)
% %sound(sin(1:3000))