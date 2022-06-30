function Framefile = n_trap_Et_movie(sol, movie_name, yrange, zrange, Vcounter)
% find the time
par = sol.par;
E_CBtrap = par.dev.E_CBtrap(2:end, :);
n_trap_Et = dfana.calcn_trap_Et(sol, 'whole');
Vapp = dfana.calcVapp(sol);

figure(400)
for i = 1:length(sol.t)
    clf
    surf(sol.x, E_CBtrap, log10(squeeze(n_trap_Et(i,:,:)))')
    s1 = gca;
    ylabel('Energy (eV)')
    xlabel('Position (cm)')
    zlabel('log10(n-trap)')
    %set(s1,'ZScale','log');
    ylim([yrange(1), yrange(2)])
    zlim([log10(zrange(1)), log10(zrange(2))])
    set(s1,'color','w');
    
    % Voltage counter
    if Vcounter
        dim = [0.7 0.7 .3 .2];
        Vnow = round(Vapp(i), 2, 'decimal');
        anno = ['V = ', num2str(Vnow), ' V'];
        T = annotation('textbox', dim, 'String', anno, 'FitBoxToText','on');
        T.EdgeColor = 'none';
        T.FontSize = 16;
        drawnow
    end
    
    Framefile(i) = getframe(gcf);
end

moviewrite(Framefile, movie_name);

function moviewrite(Framefile, movie_name)
% Write a frame file to an avi movie
% name is a string with the desired filename- do NOT include .avi extension

% name is a string for the final video
movie_name = [movie_name, '.mp4'];

% Write to file
myVideo = VideoWriter(movie_name, 'MPEG-4');
myVideo.FrameRate = 10;
open(myVideo);
writeVideo(myVideo, Framefile);
close(myVideo);

end

end