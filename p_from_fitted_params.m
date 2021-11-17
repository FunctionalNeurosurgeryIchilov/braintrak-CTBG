%
function p = p_from_fitted_params(bt_model,p,fitted_params)

params = bt_model.fullgab2full_params(fitted_params);

if any(strcmp(bt_model.param_names,'Gsn'))
    p.gab = fitted_params(1:8);
else
    p.gab = [fitted_params(1:5) bt_model.params_typical.gab(6) fitted_params(6:7)]; %set Gsn (instead on Gsn=1). Affects the total power level, provides better nus roots
end
p.gabcd = params(1:5);
p.alpha(:) = params(6);
p.beta(:) = params(7);
p.t0 = params(8);
p.taues = p.t0/2;
p.tause = p.t0/2;
p.emg_a = params(9);
