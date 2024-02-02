ps = parallel.Settings;
ps.Pool.AutoCreate = false; %do not autocreate parpool when encountering a |parfor|
ps.Pool.IdleTimeout = Inf;  %do not shutdown parpool after Inf idle time
%Best to set this to the number of PHYSICAL cores, which is 4 on Intel i7
%(But have 8 logical cores due to hyperthreading)
parpool('Processes', 4);