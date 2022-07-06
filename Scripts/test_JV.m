%Function which returns the stats from a JV sweep for use with ga

%%
function stats = test_JV(par, ions, mob, LUMO, v_surf)

    %reset device parameters
    par.Phi_EA(5) = -LUMO;
    par.Phi_IP(5) = -(LUMO+2);
    par.EF0(5) = -(LUMO+1);
    par.Et(5) = -(LUMO+1);
    par.mu_n(5) = mob;
    par.mu_p(5) = mob;
    par.sp(4) = v_surf;
    if abs(par.Phi_right) < LUMO
        par.Phi_right = -LUMO;
    end

    par = refresh_device(par);

    %find eqm solution
    try
        eqm = equilibrate(par);
    catch
        warning('Parameter combination unsucessful. Could not find eqm solution');
        stats = 0;
        return
    end

    %select eqm_sol with or without ion motion
    if ions == 1
        eqm = eqm.ion;
    elseif ions == 0 
        eqm = eqm.el;
    end

    %do JV
    try
        CV_sol = doCV(eqm, 1.15, -0.3, 1.3, -0.3, 1e-3, 1, 321);
    catch
        warning('Parameter combination unsucessful. doCV failed.');
        stats = 0;
        return
    end

    %Find CV stats
    try
        stats = CVstats(CV_sol);
    catch
        warning('Could not extract JV parameters')
        stats = 0;
        return
    end

    %Calculate QFLS at SC
    num_start = sum(CV_sol.par.layer_points(1:2))+1;
    num_stop = num_start + CV_sol.par.layer_points(3)-1;
    x = CV_sol.par.x_sub;
    d = CV_sol.par.d(3);
    [~, ~, Efn, Efp] = dfana.calcEnergies(CV_sol);
    QFLS_SC = trapz(x(num_start:num_stop), Efn(31, num_start:num_stop)-Efp(31,num_start:num_stop))/d;
    stats.QFLS_SC = QFLS_SC;
  
    % Calculate QFLS at OC 
    %Need to find time point where voltage is closest to Voc first
    Vapp = dfana.calcVapp(CV_sol);
    Voc = stats.Voc_f;
    OC_time = find(abs(Voc-Vapp) == min(abs(Voc-Vapp)),1);
    
    QFLS_OC = trapz(x(num_start:num_stop), Efn(OC_time, num_start:num_stop)-Efp(OC_time,num_start:num_stop))/d;
    stats.QFLS_OC = QFLS_OC;

end


