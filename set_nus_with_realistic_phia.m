%
%assign nus according to the most realistic set of phia
%set empty if no realistic phia found
%
%inputs:
%   params
%   typical_phia - phia vals of a typical fit
%   prints_flg
%
function [params, is_realistic] = set_nus_with_realistic_phia(params,typical_phia,prints_flg)

    max_phia_literature = [20 40 20]; %values are based on the literature
    if nargin<3
        prints_flg = 0;
    end
    if nargin<2 
        typical_phia = max_phia_literature/2;
    end
    lim_phia = [eps eps eps; max([max_phia_literature*1.5;typical_phia*2])];
    
    is_realistic = true;
    
    if isempty(params.phia)
        params.nus = [];
        is_realistic = false;
        return;
    end
        
    for iPhia = 1:size(params.phia,1)
        err(iPhia) = sum((params.phia(iPhia,:)-typical_phia).^2);
        if any(params.phia(iPhia,:)<lim_phia(1,:)) || any(params.phia(iPhia,:)>lim_phia(2,:))
            err(iPhia) = inf;
        end
    end
    [M,inx] = min(err);
    params.phia = params.phia(inx,:);
    params.nus = params.nus(inx,:);
    if prints_flg
        disp(['Phi_e, Phi_r, Phi_s: ' num2str(params.phia)]);
        disp(['XYZ: ' num2str(params.xyz)]);
    end
    
    if isinf(M)
        if prints_flg
            warning('Can''t find realistic phia!')
        end
        params.phia = [];
        params.nus = [];
        is_realistic = false;
    end

end
