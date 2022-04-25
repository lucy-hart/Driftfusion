Ncat = 4;
y_var = 9;

J_ion = dfana.calcJ(solCV_ion{y_var, Ncat});
J_el = dfana.calcJ(solCV_el{y_var, Ncat});
Vapp = dfana.calcVapp(solCV_ion{y_var, Ncat});

figure(009)
plot(Vapp, J_el.tot(:,1), 'blue', Vapp, J_ion.tot(:,1), 'red')
ylim([-0.03,0.01])
xlim([-0.15,1.2])

%%
path = 'C:\Users\ljh3218\OneDrive - Imperial College London\PhD\IonEfficiency\';
filename = 'tld_symmetric_Vbi_vs_Ncat_DopedTLs_srh.mat';
load(append(path,filename));

%%

Ncat = 1;
y_var = 4;

solCV_ion_single = doCV(soleq{y_var, Ncat}.ion, 1, -0.1, 0.4, -0.1, 1e-3, 1, 241);
%solCV_el_single = doCV(soleq{y_var, Ncat}.el, 1, -0.2, 0.4, -0.2, 1e-4, 1, 241);

  