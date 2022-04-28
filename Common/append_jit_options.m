function nlp_opts = append_jit_options(nlp_opts,optim_flag)

nlp_opts.jit=true;
jit_options = struct('flags', optim_flag, 'verbose', true, 'compiler', 'ccache gcc');
nlp_opts.compiler='shell';
nlp_opts.jit_options=jit_options;
nlp_opts.jit_temp_suffix=false;

end