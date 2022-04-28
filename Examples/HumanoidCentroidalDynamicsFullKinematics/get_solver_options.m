function s_opts = get_solver_options(solver_name)
switch solver_name       
    case 'ipopt'
        s_opts = struct('max_iter',1000,...
            'max_cpu_time',1200,...
            'tol', 1e-4,... % (1e-6), 1e-4 works well
            'acceptable_tol', 1e-4,... % (1e-4)
            'constr_viol_tol', 1e-4,... % (1e-6), 1e3 works well
            'acceptable_iter', 5,... % (15), % 5 works well
            'linear_solver','ma27',...
            'mu_strategy','adaptive',...
            'print_level',3,...% 5 is default
            'file_print_level',5,...% 5 is default %2 prints nothing
            'bound_frac',1e-2,... % 1e-2
            'bound_push',1e-2,...
            'bound_relax_factor', 1e-10,...
            'warm_start_init_point','no');%  
    otherwise
        
end
end