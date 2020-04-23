function m = organVolM(age, sex, height, bodyMass, BSA)
%UNTITLED6 Summary of this function goes here
%   Detailed explanation goes here
 
 % mass (kg)
 m =[exp(.067*BSA-.0025*age-.38*sex+1.7);  % 1.blood
     exp(0.028*height+0.0077*age-5.6);     % 2.lung
     exp(0.87*BSA - 0.0014*age - 1);       % 3.liver
     1.05;                                 % 4.stomach, not from Stader
     exp(1.13*BSA-3.93);                   % 5.spleen
     0.103;                                % 6.pancreas
     3E-6*height^2.49;                     % 7.gut
     1.025;                                % 8.adrenal, not from Stader
     exp(0.024*height-1.9);                % 9.bone
     exp(-0.0058*age-0.37*sex+1.13);       %10.skin, dermis
     -0.00038*age-0.056*sex+0.33;          %11.kidney
     0.34*BSA+0.0018*age-0.36;             %12.heart
     0.68*bodyMass-0.56*height+6.1*sex+65; %13.adipose
     17.9*BSA-0.0667*age-5.68*sex-1.22;    %14.muscle
     exp(-0.0075*age+0.0078*height-0.97)]; %15.brain
 
 m=(bodyMass/sum(m)).*m; %adjusting for mass balance wrt ratios calculated by m
end