% get linear approximation on triangle
% PARAM_IN:
%   fld      - [M, 1], field values (error or indicator) on simplexes,
%   type     - STR, marker type:
%              'mean'     - marking by mean value with kappa weight,
%              'elem_num' - marking by proportion (kappa) of elements with the greatest error,
%              'greedy'   - marking by greedy algorithm (by default)
%   kappa    - DBL, marker weight parameter, by default = 0.3
% PARAM_OUT:
%   fld_bool - [M, 1], marked field
function fld_bool = marker(fld, type, kappa)
    arguments
        fld (1, :) {mustBeNumeric}
        type (1, :) char {mustBeMember(type,{'mean','elem_num','greedy'})} = 'greedy'
        kappa (1, :) {mustBeNumeric} = 0.3
    end

    fld_bool = zeros(size(fld));
    
    if strcmp(type, 'mean')
        i1 = find(fld > kappa*mean(fld(:)));            
    elseif strcmp(type, 'elem_num')
        if ~(0 < kappa && kappa < 1)
            error('Error: for marker type "elem_num" kappa should be in [0,1)!') 
        end
        
        [~, sorted_inds] = sort(fld(:), 'descend');
        i1 = sorted_inds(1:round(kappa*length(fld(:))));
    elseif strcmp(type, 'greedy')
        if ~(0 < kappa && kappa < 1)
            error('Error: for marker type "greedy" kappa should be in [0,1)!') 
        end
        
        S = kappa*sum(fld(:));
        [sorted_fld, sorted_inds] = sort(fld(:), 'descend');
        k = 1;
        while sum(sorted_fld(1:k)) < S
            k = k + 1;
        end
        i1 = sorted_inds(1:(k-1));
    else
        error('Unknown marker type')
    end
    
    fld_bool(i1) = ones(size(i1));

end