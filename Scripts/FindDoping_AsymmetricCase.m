%File to plot the dark doping (?) for the asmmetic offset case
%Assumes you have already run IonEfficiency_v6
%Need to have actually run it thoug as I take the doping density from the
%dark eqm file and I don't save these 

%%
N = zeros(n_ion_concs, n_values);
P = zeros(n_ion_concs, n_values);

for i = 1:n_ion_concs
    if i == n_ion_concs
        par = soleq{n_ion_concs,j}.el.par;
    else
        par = soleq{i,j}.ion.par;
    end
    %Shift x coordiantes so centre is at the centre of the active layer
    %Do this as I want to find the dipole moment of the ionic charge
    xstart = sum(par.layer_points(1:2)) + 1;
    xstop = sum(par.layer_points(1:3));
    x = par.xx(xstart:xstop+1);
    for j = 1:n_values 
        if i == n_ion_concs
            n = soleq{i,j}.el.u(:,xstart:xstop+1,2);
            p = soleq{i,j}.el.u(:,xstart:xstop+1,3);
        else
            n = soleq{i,j}.ion.u(:,xstart:xstop+1,2);
            p = soleq{i,j}.ion.u(:,xstart:xstop+1,3);
        end
        N(i,j) = trapz(x,n(end,:));
        P(i,j) = trapz(x,p(end,:));
    end
end


%%
figure('Name', 'Doping')
Colours = parula(n_ion_concs-1);

for i = 1:n_ion_concs
    if i == n_ion_concs
        semilogy(Delta_TL, N(i,:), 'Color', 'Black', 'Marker', 'o')
        hold on
        semilogy(Delta_TL, P(i,:), 'Color', 'Black', 'Marker', 'x')
    else
        semilogy(Delta_TL, N(i,:), 'Color', Colours(i,:), 'Marker', 'o')
        hold on
        semilogy(Delta_TL, P(i,:), 'Color', Colours(i,:), 'Marker', 'x')
    end
end

xlabel('HTL Offset (eV)')
ylabel('Integrated Bulk Carriers')

        


