%Script to calculate recombination losses through SRH versus surface
%recombination
%Load in data first!
%Using v5_doped_vsr_offset_0p2eV
%Need to get the stats array before running this script

Ion_Conc = [1e15 5e15 1e16 5e16 1e17 5e17 1e18 0];
n_ion_concs = length(Ion_Conc);

v_sr = [0.1 1 10 100 1000 1e4 1e5];
n_recom  = length(v_sr);

J_surf = zeros(n_ion_concs, n_recom);
J_SRH = zeros(n_ion_concs, n_recom);

x = solCV{1,1}.par.x_sub;

for i = 1:n_ion_concs
    for j = 1:n_recom
        loss_currents = dfana.calcr(solCV{i,j},'sub');
        Vapp = dfana.calcVapp(solCV{i,j});
        end_value = cast((length(Vapp)-1)/2, 'int32');
        J_srh = e*trapz(x, loss_currents.srh, 2)';
        J_vsr = e*trapz(x, loss_currents.vsr, 2)';
        J_SRH(i,j) = interp1(Vapp(1:end_value), J_srh(1:end_value), Stats_array(i,j,2));
        J_surf(i,j) = interp1(Vapp(1:end_value), J_vsr(1:end_value), Stats_array(i,j,2));
    end
end