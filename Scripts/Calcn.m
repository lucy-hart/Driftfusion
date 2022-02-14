%Code to calculate charge stored on device active layers

num_start = sum(CV_solutions_ion{1}.par.layer_points(1:2))+1;
num_stop = num_start + CV_solutions_ion{1}.par.layer_points(3)-1;
x = CV_solutions_ion{1}.par.x_sub;
d = CV_solutions_ion{1}.par.d(3);
n = zeros(num_values,3);
p = zeros(num_values,3);

for z=1:3
    n(:,z) = trapz(x(num_start:num_stop), CV_solutions_ion{z}.u(:,num_start:num_stop,2),2)/d;
    p(:,z) = trapz(x(num_start:num_stop), CV_solutions_ion{z}.u(:,num_start:num_stop,3),2)/d;
end

n_bar = (n.*p).^0.5;