%% Define Global Variables
% Bio-data:
global M_liver;
global V_blood;
global V_liver;
global V_peripheral;
global V_gut;
global BW;

% Diet Model Constants:
global K_FCEX; %unused
global K_FCCM;
global K_CMGB; %unused
global K_BSLG;
global K_FCLG;
global K_BSGL; %unused
global K_BSEX; %unused
global K_CMBL;
global K_VLDLLB;
global K_LDLBP;
global K_LDLBL;
global K_HDLeLB;
global K_HDLfBL;
global K_HDLfFC;
global K_HDLeBP;
global K_HDLfPB;
global K_CMFC;
global K_FCVLDL;
global K_FCHDLe;
global K_LDLFCL;
global K_FCBS;
global K_LDLFCP;
global K_HDLefP;
global K_FCPloss;

global V_FCGmax;
global K_FCGt;
global V_FCBmax;
global K_FCBt;
global V_FCPmax;
global K_FCPt;

global K_CETP; %Cholesterylester transfer protein rate, [mM/day]
global V_max_LPL; %Lipoprotein lipase V_max (M/s)
global Km_LPL; %Lipoprotein lipase Km (M)

% HMG-CoA Regulation Model Constants ALL UNUSED
global V_m; %HMG-CoA synthesis rate constant, [mM/day]
global delta; %HMG-CoA removal rate, [1/day]
global V_m_HMGCoA; %Michaelis for HMG-CoA enzyme (reductase) activity, [mM/day]
global K_m_HMGCoA; %Michaelis constant for HMG-CoA enzyme (reductase), [mM]
global C_s_liver; %Steady state liver cholesterol pool, [mg/(g/mol)] = [mmol]

%% Bio-data
%Body Parameters FROM JENKINS DATA
age = 59.2; %Age, [years old]
sex = 1; %Sex, 1 = Male; 2 = Female
BW = Human_Body_Weight(age,sex); %Body Weight, [kg] 
height = 180; %Height, [cm]

%BSA = body surface area 
%Du Bois D, Du Bois EF (Jun 1916). "A formula to estimate the approximate surface area if height and weight be known". Archives of Internal Medicine 17 (6): 863-71. PMID 2520314. Retrieved 2012-09-09.
BSA = 0.007184*BW^.425*height^.725; %Body surface area, [m^2]  

%Volume calculations
m = organVolM(age, sex-1, height, BW, BSA);

% densities (kg/L)
p = [1.06;  % 1.blood, google search lol
     1.05;  % 2.lung
     1.08;  % 3.liver, [Heinemann et al]
     1.05;  % 4.stomach
     1.054; % 5.spleen
     1.045; % 6.pancreas
     1.044; % 7.gut
     1.025; % 8.adrenal
     1.99;  % 9.bone
     1.116; %10.skin, dermis (epidermis & hypodermis average out to ~same)
     1.05;  %11.kidney
     1.03;  %12.heart
     .916;  %13.adipose
     1.041; %14.muscle
     1.035];%15.brain
 
v = m./p;

%Liver Mass
if sex == 1
    liver_mass_coefficient = [4.25E-2, -1.01E-6, 1.99E-11, -1.66E-16, 4.83E-22, 0];
else
    liver_mass_coefficient = [3.34E-2, -1.89E-7, 5.34E-13, 0, 0, 0];
end
weight_matrix = [(BW*1000)^0; (BW*1000)^1; (BW*1000)^2; (BW*1000)^3; (BW*1000)^4; (BW*1000)^5];
liver_mass_fraction = liver_mass_coefficient*weight_matrix; %Liver mass fraction, [unitless]
M_liver = BW*liver_mass_fraction; %Liver mass, [kg]

%Blood Volume
% A Compendium of Drugs Used for Laboratory Animal Anesthesia, Analgesia, Tranquilization and Restraint Archived June 6, 2011, at the Wayback Machine. at Drexel University College of Medicine.
V_blood = 77*BW/1000; %Blood volume (77mL/kg BW), [L]          

%Liver volume
%Johnson, T.N., et al.(2005).Changes in Liver Volume from Birth to Adulthood:A Meta-Analysis
%Authors recommend equation 3 for liver volume prediction using BSA;
V_liver = (sqrt(BSA)*0.72+0.171)^3; %Liver volume, [L]

V_gut = v(7);

v(1) = V_blood;
v(3) = V_liver;
 
V_peripheral = sum(v) - v(3) - v(7) - v(1);

%% LDL-C Regulation Model Constants UNUSED

% V_max_CETP = 76E-9*60*24; %Cholesterylester transfer protein V_max, [mmol/day] (Source: Physical and kinetic characterization of recombinant human cholesteryl ester transfer protein, original value is 76 [pmol/min])
% K_m_HDL = 2000E-6; %Cholesterylester transfer protein K_m_HDL, [mM] (Source: Physical and kinetic characterization of recombinant human cholesteryl ester transfer protein, original value is 2000 [nM])
% K_m_LDL = 700E-6; %Cholesterylester transfer protein K_m_LDL, [mM] (Source: Physical and kinetic characterization of recombinant human cholesteryl ester transfer protein, original value is 700 [nM])

% V_max_CETP = 37/1000*24; %Cholesterylester transfer protein V_max, [mM/day] (Source: High density lipoprotein is an inappropriate substrate for hepatic lipase in postmenopausal women)
% K_m_HDL = 45/1000; %Cholesterylester transfer protein K_m_HDL, [mM] (Source: High density lipoprotein is an inappropriate substrate for hepatic lipase in postmenopausal women)

K_CETP = 0.0125;

V_max_LPL = 0.7; 
Km_LPL = 6.09E-3*1000; %Lipoprotein lipase Km [mM] (M) (Source: Biochemical Analysis of the Lipoprotein Lipase Truncation Variant, LPLS447X, Reveals Increased Lipoprotein Uptake.)

% K_FCGL = 1.68; %k6/C(2) 1.68E-7 initial in mg (Source: A whole-body mathematical model of cholesterol metabolism and its age-associated dysregulation) [1/day];
% K_FCEX = 0.02; %Intestinal cholesterol elimination constant 0.3, [1/day] (According to Masood) not really a source
K_FCCM = 0.23; % 0.0444 Assume same value as K_CMGB NO SOURCE [1/day]
% K_CMGB = 0.6; %B1 = 0.0444 intestine-related VLDL production rate constant NO SOURCE [1/day]
K_BSLG = 0.7; % [1/day]
K_FCLG = 0.5; % [1/day]
% K_BSGL = 0.1; %k3/C(3) initial in mg 0.0091863 (Source: A whole-body mathematical model of cholesterol metabolism and its age-associated dysregulation) [1/day]
% K_BSEX = 0.1; %k4/C(3) initial in mg 0.00184 (Source: A whole-body mathematical model of cholesterol metabolism and its age-associated dysregulation) [1/day]
K_CMBL = 0.44; % [1/day]
K_VLDLLB = 0.1; % [1/day]
K_LDLBP = 0.0065; % 0.11650969 Peripheral LDL uptake from blood (Source: A physiologically based in silico kinetic model predicting plasma cholesterol concentrations in humans, using formula 13 assuming 1st order rate) [1/day]
K_LDLBL = 0.036; % Uptake of LDL into liver initial: 2.7E-3 (Source: An Integrated Mathematical Model of Cellular Cholesterol Biosynthesis and Lipoprotein Metabolism) [1/day]
K_HDLeLB = 0.05;
K_HDLfBL = 0.0175;
K_HDLfFC = 0.4;
K_HDLeBP = 0.285;
K_HDLfPB = 0.4; %Transport of HDL out of peripheral tissues and into the blood stream
K_CMFC = 0.8;
K_FCVLDL = 0.115; %k10/C(13)1E-6 VLDL synthesis rate (Source: A whole-body mathematical model of cholesterol metabolism and its age-associated dysregulation) [1/day]
K_LDLFCL = 0.65; %Conversion of LDL to free cholesterol (Source: An Integrated Mathematical Model of Cellular Cholesterol Biosynthesis and Lipoprotein Metabolism) [1/day]
K_FCHDLe = 0.032; %k11 5E-5 (Source: A whole-body mathematical model of cholesterol metabolism and its age-associated dysregulation) [1/day]
K_FCBS = 0.00092; %k5/C(14) = 0.069 initial in mg (Source: A whole-body mathematical model of cholesterol metabolism and its age-associated dysregulation) [1/day]
K_LDLFCP = 0.9;
K_HDLefP = 0.0099;
K_FCPloss = 0.009;

% Source for all of below: A whole-body mathematical model of cholesterolmetabolism and its age-associated dysregulation
V_FCGmax = 1E2/MWc;
K_FCGt = 3.12E2/MWc;
V_FCBmax = 2E3/MWc;
K_FCBt = 5.55E4/MWc;
V_FCPmax = 5E2/MWc;
K_FCPt = 8.0342E4/MWc;

%% HMG-CoA Regulation Model Constants
V_m = 0.0175; %HMG-CoA synthesis rate constant, [mM/day]
delta = 2.60*60*24; %HMG-CoA removal rate, [1/day] (Source: A mathematical model of the sterol regulatory element binding protein 2 cholesterol biosynthesis pathway, original value in [1/min])
V_m_HMGCoA = 17.28; %Michaelis for HMG-CoA enzyme (reductase) activity, [mM/day]
K_m_HMGCoA = 4/1000; %Michaelis constant for HMG-CoA enzyme (reductase), [mM] (Source: Corsini, A., Maggi, F. M., & Catapano, A. L. (1995). Pharmacology of competitive inhibitors of HMg-CoA reductase. Pharmacological Research, 31(1), 9-27., original value is 4 µM)
C_s_liver = 60000/386.65/V_liver; %Steady state liver cholesterol pool, [mg/(g/mol)] = [mmol]