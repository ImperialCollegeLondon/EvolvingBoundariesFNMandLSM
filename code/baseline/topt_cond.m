%% TOPOPT
%
%   CALCULATION OF CONDITIONING NUMBER

function [cond] = topt_cond( s_mesh, s_geo, s_param, s_sol )

    nel = s_mesh.nel;
    l_u = s_sol.l_u;
    l_dofel = s_mesh.l_dofel;
    con = s_mesh.con;
    x = s_mesh.x;
    ngpt = s_mesh.ngpt;
    xgp = s_mesh.xgp;
    ndim = s_mesh.ndim;
    nnodel = s_mesh.nnodel;


    % conditioning number
    gradv = zeros(ngpt * nel, 1);
    numgp = 1;
    for iel = 1:nel
        xel = x(con(iel,:), :) ;
        phiel = l_u(l_dofel(iel, :), 1);
        for igp = 1:ngpt
            xi = xgp(igp,1);
            eta = xgp(igp,2);
            % shape functions
            [n, dn] = shapefuncquad( ndim, nnodel, xi, eta ) ;
            % jacobian
            J = dn * xel ; dnxy = J\dn;
            
            % condition number of the level set
            gradv(numgp,1) = norm( dnxy*phiel );
            
            numgp = numgp + 1;
        end
    end
    
    cond = 0;
    for igp = 1:size(gradv,1)
        cond = cond + (gradv(igp,1) - 1)^2 ;
    end
    cond = sqrt(cond/size(gradv,1)) ;

end

%%
function [n, dn] = shapefuncquad( ndim, nnodel, xi, eta )
% calculates the shape functions and derivatives
% for a isoparametric quadrilateral
    
    n = zeros(1, nnodel) ;
    n(:,1) = (1/4) * (1 - xi) * (1 - eta) ;
    n(:,2) = (1/4) * (1 + xi) * (1 - eta) ;
    n(:,3) = (1/4) * (1 + xi) * (1 + eta) ;
    n(:,4) = (1/4) * (1 - xi) * (1 + eta) ;
    
    dn = zeros(ndim, nnodel) ;
    dn(1,1) = -(1/4) * (1 - eta);
    dn(1,2) = +(1/4) * (1 - eta) ;
    dn(1,3) = +(1/4) * (1 + eta) ;
    dn(1,4) = -(1/4) * (1 + eta) ;
    dn(2,1) = -(1/4) * (1 - xi) ;
    dn(2,2) = -(1/4) * (1 + xi) ;
    dn(2,3) = +(1/4) * (1 + xi) ;
    dn(2,4) = +(1/4) * (1 - xi) ;

end