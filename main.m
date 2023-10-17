
load('mat/pdeparams_fg.mat')
ind_obj = Indicator(gd, a, b, c, e, f, p, t);


err = ind_obj.calcError(du_dx, du_dy);
figure
err_m = Indicator.marker(err, type, kappa);
ind_obj.plotFld(err_m);
title('Error')

indr0 = ind_obj.getIndicator('AG');
figure
indr0_m = Indicator.marker(indr0);
ind_obj.plotFld(indr0_m);
title('Average Gradient indicator')

indr = ind_obj.getIndicator('MP');
figure
indr_m = Indicator.marker(indr);
ind_obj.plotFld(indr_m);
title('M_+ indicator')


N_t = length(t);
complare_flds = @(name, fld) sprintf('%s = %d/%d = %.2f%%', name,...
                            length(find(err_m == indr0_m)), N_t,...
                            length(find(err_m == indr0_m))/N_t*100);
disp('fraction of correct marking')
disp(complare_flds('AD', indr0_m))
disp(complare_flds('MP', indr_m))
