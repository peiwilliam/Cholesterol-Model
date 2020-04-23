function dC = Cholesterol_Model(t, C)

MWc = 386.65; %Molecular weight of cholesterol, [g/mol]

%% Define global variables
% Model input
global Q_diet;
global Q_syn_intestine;
global Q_syn_peripheral;
global Q_FCLBi;
global LiverModel;
global Model;
global count;
global FBA_Size;
% Bio-data:
global M_liver;
global V_blood;
global V_gut;
global V_liver;
global V_peripheral;
global BW;
% LDL-C Regulation Model Constants:
global K_FCEX;
global K_FCCM;
global K_CMGB;
global K_BSLG;
global K_FCLG;
global K_BSGL;
global K_BSEX;
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

global Q_syn_liver; %Cholesterol synthesis rate in liver, [µmol/kg body weight/day] --> Update from FBA
global K_CETP; %Cholesterylester transfer protein rate, [mM/day]
global V_max_LPL; %Lipoprotein lipase V_max (M/s)
global Km_LPL; %Lipoprotein lipase Km (M)
% HMG-CoA Regulation Model Constants
global V_m; %HMG-CoA synthesis rate constant, [mM/day]
global delta; %HMG-CoA removal rate, [1/day]
global V_m_HMGCoA; %Michaelis for HMG-CoA enzyme (reductase) activity, [mM/day]
global K_m_HMGCoA; %Michaelis constant for HMG-CoA enzyme (reductase), [mM]
global C_s_liver; %Steady state liver cholesterol pool, [mg/(g/mol)] = [mmol]

%% Define vector for system of ODE

dC = zeros(20, 1);

%% Genome-Scale Model: UNUSED

%HMG-CoA regulation in liver
%Dependent variable: (20)C_HMGCoA
%HMGCoA is used as an upper-bound constraint to feed into FBA

%Q_syn_HMGCoA = V_m/(1+(C_liver/C_s_liver))  [mM/day]
Q_syn_HMGCoA = V_m/(1+(C(14)/C_s_liver));

%ODE1: HMG-CoA Concentration in liver
%d(C_HMGCoA)/dt = [mM/day]
%Q_syn_HMGCoA - [mM/day]
%delta*C_HMGCoA [mM/day]
dC(20) = Q_syn_HMGCoA - delta*C(20);

%HMG-CoA Reductase Michaelis-Menten Kinetics
%V_HMGCoA_Reductase = V_m_HMGCoA*C_HMGCoA/(K_m_HMGCoA+C_HMGCoA) [mM/day]
V_HMGCoA_Reductase = V_m_HMGCoA*C(20)/(K_m_HMGCoA + C(20));
%Unit of flux in HMR database is [µmol/100g cell/min] (Source: Integration of clinical data with a genome-scale metabolic model of the human adipocyte)
V_constraint_ub = (V_HMGCoA_Reductase*1000/(24*60))*V_liver/(M_liver*10); %Convert upper bound unit to [µmol/min/100g liver cell]
%Ensure V_constraint_ub cannot be less than 0 (Biological limit)
if V_constraint_ub < 0
    V_constraint_ub = 0;
end

%Reaction HMR_1440 corresponds to HMG-CoA Reductase's reaction. It's also row 2997 in the vector.
%Reaction HMR_1526 corresponds to cholesterol synthesis reaction. It's also row 3023 in the vector.
%Vector c is a parameter vector that linearly combines one or more reaction
%fluxes to form what is termed the objective function (c = # of reactions)
%Vector b is a vector of known metabolic exchanges (b = # of metabolites)
%Output of FBA is a particular flux distribution (v) which maximizes or minimizes the objective function and stands between upper (ub) and lower (lb) bounds [µmol/min/gDW].
%Objective function: Maximize cholesterol synthesis rate

%To improve runtime
% if mod(count,FBA_Size) == 0 || count == 1
%     %Modify Constraints 
%     LiverModel = changeRxnBounds(LiverModel, 'HMR_1440', V_constraint_ub, 'u'); %Set upper bound constraint for HMG-CoA Reductase reaction rate    
%     %Optimize Model
%     FBALiverModel = optimizeCbModel(LiverModel, 'max'); %Run FBA to get the flux
%     %Use the fluxes to update the model: liver cholesteral synthesis rate
%     Q_syn_liver = FBALiverModel.x(3023)*(M_liver*10)*(24*60)/BW; %[µmol/kg body weight/day]
%     %hold on
%     %plot(count, Q_syn_liver, '-o'); for debugging purposes
% end

% Q_syn_liver =  29.8*24/MWc; % original 29.8 mg/h (Source: Biliary lipids, cholesterol and bile synthesis: different adaptivemechanisms to dietary cholesterol in lean and obese subjects)

%% Modified Michaelis-Menten rate equation for CETP (Source: Physical and kinetic characterization of recombinant human cholesteryl ester transfer protein) UNUSED
%CE-HDL + LDL --> HDL + LDL-CE, donor substrate is HDL containing CE. The acceptor substrate is LDL.
%In this model, assume CE-HDL is HDL and LDL-CE is LDL

%V_CETP/V_max_CETP = (HDL*LDL)/(K_m_LDL*HDL+K_m_HDL*LDL+HDL*LDL) %[mmol/day]
%R_CETP = V_max_CETP*C(8)*C(6)/(K_m_LDL*C(6)+K_m_HDL*C(8)+C(8)*C(6))/V_blood; %[mM/day]

%% Michelis Menten rate for CETP (Source: High density lipoprotein is an inappropriate substrate for hepatic lipase in postmenopausal women) UNUSED

% R_CETP = V_max_CETP*C(8)/(K_m_HDL + C(8)); %[mM/day]

%% Lipoprotein Lipase Michaelis-Menten rate equation (Source: Biochemical Analysis of the Lipoprotein Lipase Truncation Variant, LPLS447X, Reveals Increased Lipoprotein Uptake.)

R_LPL = V_max_LPL*C(5)/(Km_LPL + C(5));

%% Cholesterol synthesis in the intestines

Q_syn_intestine = V_FCGmax/(1 + C(1)*V_gut/K_FCGt)/V_gut;

%% Cholesterol synthesis in the peripheral organs

Q_syn_peripheral = V_FCPmax/(1 + C(19)*V_peripheral/K_FCPt)/V_peripheral;

%% Cholesterol release into bile UNUSED

% Q_FCLBi = V_FCBmax/(1 + K_FCBt/C(14))/V_liver;

%% K_FCEX

K_FCEX = 0.5/C(1)*(Q_diet/V_gut + K_FCLG*C(14)*V_liver/V_gut);
K_CMGB = K_FCEX;

%% K_BSEX

K_BSEX = 0.015/C(3)*(K_BSLG*C(15)*V_liver/V_gut);
K_BSGL = 0.985/C(3)*(K_BSLG*C(15)*V_liver/V_gut);
%% Model
%(1)C_FCG, (2)C_CMG, (3)C_BSG, (4)C_CMB, (5)C_VLDLB, (6)C_LDLB, (7)C_HDLeB, 
%(8)HDLfB, (9)C_CML, (10)C_VLDLL, (11)C_LDLL, (12)C_HDLeL, (13)C_HDLfL,  
%(14)C_FCL, (15)C_BSL, (16)C_LDLP, (17)C_HDLeP, (18)C_HDLfP, (19)C_FCP

dC(1:19) = [Q_diet/V_gut + Q_syn_intestine  + K_FCLG*C(14)*V_liver/V_gut - K_FCEX*C(1) - K_FCCM*C(1);  %1. dFCG
            K_FCCM*C(1) - K_CMGB*C(2);                                                      %2. dCMG
            K_BSLG*C(15)*V_liver/V_gut - K_BSGL*C(3) - K_BSEX*C(3)                          %3. dBSG
            K_CMGB*C(2)*V_gut/V_blood - K_CMBL*C(4);                                        %4. dCMB
            K_VLDLLB*C(10)*V_liver/V_blood - R_LPL;                                         %5. dVLDLB
            R_LPL + K_CETP*C(8) - K_LDLBP*C(6) - K_LDLBL*C(6);                              %6. dLDLB
            K_HDLeLB*C(12)*V_liver/V_blood - K_HDLeBP*C(7);                                 %7. dHDLeB
            K_HDLfPB*C(18)*V_peripheral/V_blood - K_HDLfBL*C(8) - K_CETP*C(8);              %8. dHDLfB
            K_CMBL*C(4)*V_blood/V_liver - K_CMFC*C(9);                                      %9. dCML
            K_FCVLDL*C(14) - K_VLDLLB*C(10);                                                %10. dVLDLL
            K_LDLBL*C(6)*V_blood/V_liver - K_LDLFCL*C(11);                                  %11. dLDLL
            K_FCHDLe*C(14) - K_HDLeLB*C(12);                                                %12. dHDLeL
            K_HDLfBL*C(8)*V_blood/V_liver - K_HDLfFC*C(13);                                 %13. dHDLfL
            Q_syn_liver/V_liver + K_HDLfFC*C(13) + K_CMFC*C(9) + K_LDLFCL*C(11) - K_FCHDLe*C(14) - K_FCVLDL*C(14) - K_FCBS*C(14)/C(15) - K_FCLG*C(14); %14. dFCL
            K_FCBS*C(14)/C(15) + K_BSGL*C(3)*V_gut/V_liver - K_BSLG*C(15);                  %15. dBSL 
            K_LDLBP*C(6)*V_blood/V_peripheral - K_LDLFCP*C(16);                             %16. dLDLP
            K_HDLeBP*C(7)*V_blood/V_peripheral - K_HDLefP*C(17)*C(19);                      %17. dHDLeP
            K_HDLefP*C(17)*C(19) - K_HDLfPB*C(18);                                          %18. dHDLfP
            Q_syn_peripheral + K_LDLFCP*C(16) - K_HDLefP*C(17)*C(19) - K_FCPloss*C(19)];    %19. dFCP