%% TOPOPT
%
%   CALCULATE THE INTERSECTIONS AND UPDATE THE DENSITY FIELD

function [rho, G] = topt_intersect( s_mesh, s_geo, s_param, s_sol )

    nel = s_mesh.nel ;
    nnode = s_mesh.nnode;
    l_u = s_sol.l_u ;
    dofel = s_mesh.dofel;
    l_dofel = s_mesh.l_dofel;
    ndim = s_mesh.ndim ;
    nnodel = s_mesh.nnodel;
    con = s_mesh.con ;
    f = s_sol.fe;

    [G] = intersectmat( nel, nnode, l_u, l_dofel, ndim, nnodel, con ) ;
    [rho] = elementdensity( nel, l_u, f, dofel, l_dofel ) ;

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

%%
function [G] = intersectmat( nel, nnode, l_u, l_dofel, ndim, nnodel, con )
% creates the intersection matrix for the level set

    G = sparse(5*nel, nnode);
    intcount = 1;
    corr = [2, 4; % 1 is neighbour of 2 and 4
            1, 3; % 2 is neighbour of 1 and 3 and so on..
            4, 2;
            3, 1];
    xieta = [-1, -1;
              1, -1;
              1,  1;
             -1,  1];
    for iel = 1:nel
        phiel = l_u(l_dofel(iel,:), 1) ;
        if any(phiel <= 0) && any(phiel > 0)
            indx = find(phiel <= 0) ;
            
            if isequal(indx,[1]) || isequal(indx,[2]) || isequal(indx,[3]) || isequal(indx,[4])
                xi = xieta(indx,1) - phiel(indx,1) * ( (xieta(indx,1) - xieta(corr(indx,1),1) ) / (phiel(indx,1) - phiel(corr(indx,1),1)) ) ;
                eta = xieta(indx,2) - phiel(indx,1) * ( (xieta(indx,2) - xieta(corr(indx,2),2) ) / (phiel(indx,1) - phiel(corr(indx,2),1)) ) ;
                [n, dn] = shapefuncquad( ndim, nnodel, xi, xieta(indx, 2) ) ;
                G(intcount, con(iel,:)) = n ;
                intcount = intcount + 1;
                [n, dn] = shapefuncquad( ndim, nnodel, xieta(indx,1), eta ) ;
                G(intcount, con(iel,:)) = n ;
                intcount = intcount + 1;
                
            elseif isequal(indx,[1;2;3]) || isequal(indx,[2;3;4]) || isequal(indx,[1;3;4])
                indx = setdiff([1;2;3;4], indx) ;
                xi = xieta(indx,1) - phiel(indx,1) * ( (xieta(indx,1) - xieta(corr(indx,1),1) ) / (phiel(indx,1) - phiel(corr(indx,1),1)) ) ;
                eta = xieta(indx,2) - phiel(indx,1) * ( (xieta(indx,2) - xieta(corr(indx,2),2) ) / (phiel(indx,1) - phiel(corr(indx,2),1)) ) ;
                [n, dn] = shapefuncquad( ndim, nnodel, xi, xieta(indx, 2) ) ;
                G(intcount, con(iel,:)) = n ;
                intcount = intcount + 1;
                [n, dn] = shapefuncquad( ndim, nnodel, xieta(indx,1), eta ) ;
                G(intcount, con(iel,:)) = n ;
                intcount = intcount + 1;
                
            elseif isequal(indx,[1;2]) || isequal(indx,[3;4])
                for iint = 1:size(indx,1)
                    eta = xieta(indx(iint),2) - phiel(indx(iint),1) * ( (xieta(indx(iint),2) - xieta(corr(indx(iint),2),2) ) / (phiel(indx(iint),1) - phiel(corr(indx(iint),2),1)) ) ;
                    [n, dn] = shapefuncquad( ndim, nnodel, xieta(indx(iint),1), eta ) ;
                    G(intcount, con(iel,:)) = n ;
                    intcount = intcount + 1;
                end
            elseif isequal(indx,[2;3]) || isequal(indx,[1;4])
                for iint = 1:size(indx,1)
                    xi = xieta(indx(iint),1) - phiel(indx(iint),1) * ( (xieta(indx(iint),1) - xieta(corr(indx(iint),1),1) ) / (phiel(indx(iint),1) - phiel(corr(indx(iint),1),1)) ) ;
                    [n, dn] = shapefuncquad( ndim, nnodel, xi, xieta(indx(iint), 2) ) ;
                    G(intcount, con(iel,:)) = n ;
                    intcount = intcount + 1;
                end
            elseif isequal(indx,[1;3]) || isequal(indx,[2;4])
                for iint = 1:size(indx,1)
                    xi = xieta(indx(iint),1) - phiel(indx(iint),1) * ( (xieta(indx(iint),1) - xieta(corr(indx(iint),1),1) ) / (phiel(indx(iint),1) - phiel(corr(indx(iint),1),1)) ) ;
                    eta = xieta(indx(iint),2) - phiel(indx(iint),1) * ( (xieta(indx(iint),2) - xieta(corr(indx(iint),2),2) ) / (phiel(indx(iint),1) - phiel(corr(indx(iint),2),1)) ) ;
                    [n, dn] = shapefuncquad( ndim, nnodel, xi, xieta(indx(iint), 2) ) ;
                    G(intcount, con(iel,:)) = n ;
                    intcount = intcount + 1;
                    [n, dn] = shapefuncquad( ndim, nnodel, xieta(indx(iint),1), eta ) ;
                    G(intcount, con(iel,:)) = n ;
                    intcount = intcount + 1;
                end  
            end
            
            
        end
    end
    G(intcount:end, :) = [];

end

%%
function [rho] = elementdensity( nel, l_u, f, dofel, l_dofel )
% updates the element density that affects the stifness value
    
    rho = ones(nel,1);
    for iel = 1:nel
        phiel = l_u(l_dofel(iel, :), 1) ;
        fel = f(dofel(iel, :), 1) ;
        if any(phiel <= 0) && ~any(fel ~= 0)
            avg = sum(phiel,1)/size(phiel,1) ;
            if avg <= 0
                rho(iel,1) = 0.0001;
            end
        end
    end
end
















