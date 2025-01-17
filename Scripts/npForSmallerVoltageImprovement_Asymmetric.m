%Code to calculate charge stored on device active layers

num = 2;
nbar = cell(num,2);
pbar = cell(num,2);
pinterface = cell(num,2);
ninterface = cell(num,2);

for j = 1:2
    for i = 1:num
        sol = solCV{i,j};
        w = sol.par.d_active;
        num_start = sum(sol.par.layer_points(1:2))+1;
        num_stop = num_start + sol.par.layer_points(3)-1;
        x = sol.x;
        nbar{i,j} = trapz(x(num_start:num_stop), sol.u(:,num_start:num_stop,2),2)/w;
        pbar{i,j} = trapz(x(num_start:num_stop), sol.u(:,num_start:num_stop,3),2)/w;
        pinterface{i,j} = sol.u(:,num_stop,3);
        ninterface{i,j} = sol.u(:,num_start,2);
    end
end

%%
figure('Name', 'carirrier bar')

%which = 21+round(100*Stats_array(:,1,2));
which = [21+110, 21+100];%, 21+100];
values = zeros(3,4,num);
r = zeros(num,1);
hold on 
for i = 1:2
    for j = 1:2
        plot(j+(0.25*(i-1)), log10(nbar{i,j}(which(j))), 'Marker', 'o')
        plot(j+(0.25*(i-1)), log10(pbar{i,j}(which(j))), 'Marker', 'x')
        plot(j+(0.25*(i-1)), log10(ninterface{i,j}(which(j))), 'Marker', 's')
        plot(j+(0.25*(i-1)), log10(pinterface{i,j}(which(j))), 'Marker', 'd')
        values(j,1,i) = nbar{i,j}(which(j));
        values(j,2,i) = pbar{i,j}(which(j));
        values(j,3,i) = pinterface{i,j}(which(j));
        values(j,4,i)= ninterface{i,j}(which(j));
    end
end

hold off
 





