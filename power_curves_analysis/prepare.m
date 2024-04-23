% Prepare data

clear
load ../VASCA/blood_hemo_vasca.mat
F = [F_ASCA id'];
var_l = var_final;
save blood_hemo_vasca.mat F var_l

clear
load ../VASCA/biomark_vasca.mat
F = F_ASCA_2;
var_l = var_final;
save biomark_vasca.mat F var_l

clear
load ../VASCA/bacteria_vasca.mat
F = F_ASCA;
var_l = var_final;
save bacteria_vasca.mat F var_l