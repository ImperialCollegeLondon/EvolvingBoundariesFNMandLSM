%% 2D FEM
%
%   POST PROCESSING

function [sig, eps] = topt_post( s_mesh, s_geo, s_param, s_sol)
    
    D = constitutivematrix( s_geo.E, s_geo.nu ) ;
    nel = s_mesh.nel ;
    nstr = s_mesh.nstr ;
    dofel = s_mesh.dofel ;
    nnode = s_mesh.nnode ;
    x = s_mesh.x;
    con = s_mesh.con;
    ndim = s_mesh.ndim ;
    ngpt = s_mesh.ngpt ;
    xgp = s_mesh.xgp ;
    u = s_sol.u ;
    nnodel = s_mesh.nnodel;
    ndofel = s_mesh.ndofel;
    
    
    sig = zeros(nnode, nstr);
    eps = zeros(nnode, nstr);
    ndcount = zeros(nnode, 1) ;
    for iel = 1:nel
        
        xel = x( con(iel, :), :) ;
        uel = u( dofel(iel, :), :) ;

            
        siggp = zeros(ngpt, nstr);
        epsgp = zeros(ngpt, nstr);
        for igp = 1:ngpt
            xi = xgp(igp,1);
            eta = xgp(igp,2);
            [n, dn] = sfquad4( ndim, nnodel, xi, eta ) ;
            J = dn * xel ; dnxy = J\dn ;
            b = bmatrixquad( nnodel, dnxy, nstr, ndofel) ;
            
            epsgp(igp,:) = (b * uel)' ;
            siggp(igp,:) = (D * epsgp(igp,:)' )' ;
        end
        
        % nodal extrapolation
        for ind = 1:nnodel
            xi = sqrt(3) * sign( xgp(ind,1) );
            eta = sqrt(3) * sign( xgp(ind,2) );
            
            [n, dn] = sfquad4( ndim, nnodel, xi, eta ) ;
            
            for ist = 1:nstr
                eps(con(iel,ind), ist) = eps(con(iel,ind), ist) + n * epsgp(:, ist) ;
                sig(con(iel,ind), ist) = sig(con(iel,ind), ist) + n * siggp(:, ist) ;
            end
        end
        
        ndcount(con(iel,:), 1) = ndcount(con(iel,:), 1) + 1 ;
        
    end
    
    % averaging ;
    nds = find(ndcount);
    for ind = (nds)'
        sig(ind,:) = sig(ind,:) / ndcount(ind,1) ;
        eps(ind,:) = eps(ind,:) / ndcount(ind,1) ;
    end

end

%%
function [b] = bmatrixquad( nnodel, dnxy, nstr, ndofel )
% calculates the b matrix for a quadrilateral element
    
    b = zeros(nstr, ndofel) ;
    for ind = 1:nnodel
        b(:, (2*ind - 1):(2*ind) ) = [ dnxy(1,ind),           0 ;
                                                 0, dnxy(2,ind) ;
                                       dnxy(2,ind), dnxy(1,ind)]; 
    end
    b = sparse(b);

end

%%
function [n, dn] = sfquad4( ndim, nnodel, xi, eta )
% calculates the shape functions and derivatives
% for a isoparametric quadrilateral
    
    n = zeros(1, nnodel) ;
    n(1,1) = 0.25 * (1 - xi) * (1 - eta) ;
    n(1,2) = 0.25 * (1 + xi) * (1 - eta) ;
    n(1,3) = 0.25 * (1 + xi) * (1 + eta) ;
    n(1,4) = 0.25 * (1 - xi) * (1 + eta) ;
    
    dn = zeros(ndim, nnodel) ;
    dn(1,1) = -0.25 * (1 - eta) ;
    dn(1,2) = +0.25 * (1 - eta) ;
    dn(1,3) = +0.25 * (1 + eta) ;
    dn(1,4) = -0.25 * (1 + eta) ;
    dn(2,1) = -0.25 * (1 - xi) ;
    dn(2,2) = -0.25 * (1 + xi) ;
    dn(2,3) = +0.25 * (1 + xi) ;
    dn(2,4) = +0.25 * (1 - xi) ;

end

%%
function [D] = constitutivematrix( E, nu )
% constructs the plane stress constitutive matrix
 
    c = E/(1 - nu^2);
    D = c*[ 1, nu,          0 ;
           nu,  1,          0 ;
            0,  0, 0.5*(1-nu)];

end
