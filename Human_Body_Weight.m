function HBW = Human_Body_Weight(age,sex)
%Leucke, R, et al.(2007).Postnatal Growth Considerations for PBPK Modeling
%Body Weight Coefficients for Growth Curves Using Modified Gompertz Equations

%C = [t0 t1 t2 t3 W01 R01 alpha1 beta1 W03 R03 alpha3 beta3] --> 2 x 12 matrix
C = [0 2.949004 11 21.66022 3.42826 2.188189 1.506303 2.462428 34.35453 0.179781 0.221175 -0.0911; % Male, sex = 1
     0 2.949004 11 21.86533 3.4 2.07402 1.504154 2.672665 34.8806 0.212837 0.443786 0.31127]; %Female, sex = 2

if (age>C(sex,1)) && (age<=C(sex,2))
    HBW = C(sex,5) * exp( (C(sex,6)/C(sex,7)) * (1-exp(-C(sex,7)*age)) );
elseif (age>C(sex,2)) && (age<=C(sex,3))
    HBW = C(sex,5) * exp( (C(sex,6)/C(sex,7)) * (1-exp(-C(sex,7)*age)) ) + C(sex,8)*(age-C(sex,2));
elseif (age>C(sex,3)) && (age<=C(sex,4))
    HBW = C(sex,9) * exp( (C(sex,10)/C(sex,11)) * (1-exp(-C(sex,11)*age)) );
else
    HBW = C(sex,9) * exp( (C(sex,10)/C(sex,11)) * (1-exp(-C(sex,11)*age)) ) + C(sex,12)*(age-C(sex,4));
end

end