%Script to calculate recombination losses through SRH versus surface
%recombination
%Load in data first!
%Using v5_doped_vsr_offset_0p2eV
%Need to get the stats array before running this script

Ion_Conc = [5e15 1e16 5e16 1e17 5e17 1e18 0];
n_ion_concs = length(Ion_Conc);

% v_sr = [0.1 1 10 100 1000 1e4 1e5];
% n_recom  = length(v_sr);

Delta_TL = linspace(0, 0.5, n_values);
Delta_TL(1) = 1e-3;
n_offset  = length(Delta_TL);

J_surf = zeros(n_ion_concs, n_offset);
J_SRH = zeros(n_ion_concs, n_offset);
J_SAM1 = zeros(n_ion_concs, n_offset);
J_SAM2 = zeros(n_ion_concs, n_offset);

x = solCV{1,1}.par.x_sub;

for i = 1:n_ion_concs
    for j = 1:n_offset
        loss_currents = dfana.calcr(solCV{i,j},'sub');
        Vapp = dfana.calcVapp(solCV{i,j});
        end_value = cast((length(Vapp)-1)/2, 'int32');
        J_srh = e*trapz(x, loss_currents.srh+loss_currents.btb, 2)';
        J_vsr = e*trapz(x, loss_currents.vsr, 2)';
        J_l1 = dfana.calcj_surf_rec_lucy(solCV{i,j}).l;
        J_l2 = dfana.calcj_surf_rec(solCV{i,j}).l;
        J_SRH(i,j) = interp1(Vapp(1:end_value), J_srh(1:end_value), Stats_array(i,j,2));
        J_surf(i,j) = interp1(Vapp(1:end_value), J_vsr(1:end_value), Stats_array(i,j,2));
        J_SAM1(i,j) = -e*interp1(Vapp(1:end_value), J_l1(1:end_value), Stats_array(i,j,2));
        J_SAM2(i,j) = e*interp1(Vapp(1:end_value), J_l2(1:end_value), Stats_array(i,j,2));
    end
end