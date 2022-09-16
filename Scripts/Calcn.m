%Code to calculate charge stored on device active layers

num_stop_HTL = CV_solutions_ion{2}.par.layer_points(1);
num_start_pero = sum(CV_solutions_ion{2}.par.layer_points(1:2))+1;
num_stop_pero = num_start_pero + CV_solutions_ion{2}.par.layer_points(3)-1;
num_start_ETL = sum(CV_solutions_ion{2}.par.layer_points(1:4))+1;

x = CV_solutions_ion{2}.x;
area = 0.045; 
e = 1.6e-19;

Qn = zeros(3,2);
Qp = zeros(3,2);

Qn(1,1) = e*area*trapz(x(1:num_stop_HTL), eqm_solutions_dark{2}.ion.u(end,1:num_stop_HTL,2),2);
Qp(1,1) = e*area*trapz(x(1:num_stop_HTL), eqm_solutions_dark{2}.ion.u(end,1:num_stop_HTL,3),2);
Qn(2,1) = e*area*trapz(x(num_start_pero:num_stop_pero), eqm_solutions_dark{2}.ion.u(end,num_start_pero:num_stop_pero,2),2);
Qp(2,1) = e*area*trapz(x(num_start_pero:num_stop_pero), eqm_solutions_dark{2}.ion.u(end,num_start_pero:num_stop_pero,3),2);
Qn(3,1) = e*area*trapz(x(num_start_ETL:end), eqm_solutions_dark{2}.ion.u(end,num_start_ETL:end,2),2);
Qp(3,1) = e*area*trapz(x(num_start_ETL:end), eqm_solutions_dark{2}.ion.u(end,num_start_ETL:end,3),2);

Qn(1,2) = e*area*trapz(x(1:num_stop_HTL), CV_solutions_ion{2}.u(139,1:num_stop_HTL,2),2);
Qp(1,2) = e*area*trapz(x(1:num_stop_HTL), CV_solutions_ion{2}.u(139,1:num_stop_HTL,3),2);
Qn(2,2) = e*area*trapz(x(num_start_pero:num_stop_pero), CV_solutions_ion{2}.u(139,num_start_pero:num_stop_pero,2),2);
Qp(2,2) = e*area*trapz(x(num_start_pero:num_stop_pero), CV_solutions_ion{2}.u(139,num_start_pero:num_stop_pero,3),2);
Qn(3,2) = e*area*trapz(x(num_start_ETL:end), CV_solutions_ion{2}.u(139,num_start_ETL:end,2),2);
Qp(3,2) = e*area*trapz(x(num_start_ETL:end), CV_solutions_ion{2}.u(139,num_start_ETL:end,3),2);
