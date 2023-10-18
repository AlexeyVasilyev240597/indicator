
%% loading task params exported from pdetool
load('mat/pdeparams_lg.mat')
ind_obj = Indicator(gd, a, b, c, e, f, p, t);

%% estimating the error
err   = ind_obj.calcError(du_dx, du_dy);
indr0 = ind_obj.getIndicator('AG');
indr  = ind_obj.getIndicator('MP');

%% plotting marked fields of the error and indicators
err_m = Indicator.marker(err);
ind_obj.plotFld(err_m, 'Error');

indr0_m = Indicator.marker(indr0);
ind_obj.plotFld(indr0_m, 'Average Gradient indicator');

indr_m = Indicator.marker(indr);
ind_obj.plotFld(indr_m, 'M_+ indicator');

%% comparing indicators
N_t = length(t);
complare_flds = @(name, fld) sprintf('%s = %d/%d = %.2f%%', name,...
                             length(find(err_m == fld)), N_t,...
                             length(find(err_m == fld))/N_t*100);
disp('fraction of correct marking')
disp(complare_flds('AD', indr0_m))
disp(complare_flds('MP', indr_m))
