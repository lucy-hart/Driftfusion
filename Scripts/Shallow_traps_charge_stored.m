par = pc('Input_files/SAM_MAFACsPbIBr_C60.csv');
Et = [-4.715 -4.5 -4.25 -4.05];
suns = logspace(-1,1,5);
sols = cell(length(Et),length(suns));
n_SC = zeros(length(Et),length(suns));
p_SC = zeros(length(Et),length(suns));

for i = 1:length(Et)
    par.Et(1) = Et(i);
    par = refresh_device(par);
    eqm = equilibrate(par);
    n_SC_eqm = trapz(x, eqm.ion.u(end,:,2));
    p_SC_eqm = trapz(x, eqm.ion.u(end,:,3));
    for j = 1:length(suns)
        sols{i,j} = doCV(eqm.ion, suns(j), -0.2, 1.2, -0.2, 1e-4, 1, 281);
        x = sols{i,j}.x;
        n_SC(i,j) = trapz(x, sols{i,j}.u(21,:,2))-n_SC_eqm;
        p_SC(i,j) = trapz(x, sols{i,j}.u(21,:,3))-p_SC_eqm;
    end
end

%%
figure('Name','n vs suns')
hold on

for i = 1:length(Et)
    loglog(suns, (n_SC(i,:).*p_SC(i,:)).^0.5, 'Marker', 'o')
end
hold off