function [F_design] = FoS( Max_load, Uncertainty )
prompt = 'Enter desired failure probability: ';
P_dsr = input(prompt);
F_dsr = icdf('normal', P_dsr, Max_load, Uncertainty);
n = (Max_load - F_dsr)/Uncertainty;
FOS = 1/(1-n*(Uncertainty/Max_load));
F_design = Max_load/FOS;

end

