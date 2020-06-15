%% TOPOPT
%
%   CALCULATION OF THE EXTERNAL FORCE VECTOR

function [fe] = topt_fext(mode, s_mesh, s_geo, s_param, s_sol)

    ndof = s_mesh.ndof;
    x = s_mesh.x;
    Lx = s_geo.Lx;
    Ly = s_geo.Ly;
    divy = s_mesh.divy;
    dofnode = s_mesh.dofnode;
    F = s_param.F;
    nel = s_mesh.nel;
    l_u = s_sol.l_u ;
    l_dofel = s_mesh.l_dofel ;
    l_ndofel = s_mesh.l_ndofel ;
    l_ndof = s_mesh.l_ndof;
    ngpt = s_mesh.ngpt;
    xgp = s_mesh.xgp;
    wgp = s_mesh.wgp;
    ndim = s_mesh.ndim ;
    nnodel = s_mesh.nnodel ;
    con = s_mesh.con;
    h = s_param.h;
    vn = s_sol.vn;
    dt = s_param.dt;
    dtau = s_param.dtau;
    rinit = s_param.rinit;
    ubc = s_param.ubc;
    
    if mode == 1
        
        fe = pointload( ndof, x, Lx, Ly, dofnode, F, ubc );
%         fe = lineload( ndof, x, Lx, Ly, divy, dofnode );
    
    elseif mode == 2
        
        fe = zeros(l_ndof, 1);
        for iel = 1:nel
            phiel = l_u(l_dofel(iel,:), 1) ;
            xel = x(con(iel,:), :) ;
            vnel = vn(con(iel,:), :);
            
            [fel] = levelsetforce( l_ndofel, ngpt, iel, xgp, wgp, ndim, nnodel, xel, phiel, h, vnel, dt, dtau, rinit ) ;
            fe( l_dofel(iel,:) , 1) = fe( l_dofel(iel,:) , 1) + fel;
        end
        
    end

end

%%
function [fp] = pointload( ndof, x, Lx, Ly, dofnode, F, ubc )
% apply a point load in the y direction of the
% middle node of the right boundary

    fp = zeros(ndof,1);

    if ubc == 1
        
        nodeindx = find( x(:,1) == Lx );
        nodeindx2 = find( abs(x(nodeindx,2) - Ly/2) == min(abs(x(nodeindx,2) - Ly/2)) );
        nodeindx = nodeindx(nodeindx2);
        nodeindx = nodeindx(1);
        dofindx = dofnode(nodeindx,2);
        fp(dofindx,1) = F;
        
    elseif ubc == 2
        
        nodeindx = find( x(:,1) == Lx );
        nodeindx2 = find( abs(x(nodeindx,2) - 0) == min(abs(x(nodeindx,2) - 0)) );
        nodeindx = nodeindx(nodeindx2);
        nodeindx = nodeindx(1);
        dofindx = dofnode(nodeindx,2);
        fp(dofindx,1) = F;
        
    elseif ubc == 3
        
        nodeindx = find( x(:,2) == 0 );
        nodeindx2 = find( abs(x(nodeindx,1) - Lx/2) == min(abs(x(nodeindx,1) - Lx/2)) );
        nodeindx = nodeindx(nodeindx2);
        nodeindx = nodeindx(1);
        dofindx = dofnode(nodeindx,2);
        fp(dofindx,1) = F;
        
    elseif ubc == 4
        
        nodeindx = find( x(:,2) == 0 );
        nodeindx2 = find( abs(x(nodeindx,1) - Lx) == min(abs(x(nodeindx,1) - Lx)) );
        nodeindx = nodeindx(nodeindx2);
        nodeindx = nodeindx(1);
        dofindx = dofnode(nodeindx,2);
        fp(dofindx,1) = F;
        
    end

    
end

%%
function [fl] = lineload( ndof, x, Lx, Ly, divy, dofnode )
    
    h = Ly/divy;

    fl = zeros(ndof,1);
    nodeindx = find( x(:,1) == Lx );
    dofindx = dofnode(nodeindx,1);
    fl(dofindx,1) = h;
    
    nodeindx2 = find( x(:,2) == Ly ) ;
    nodeindx3 = find( x(:,2) == 0 ) ;
    nodeindx = [intersect(nodeindx, nodeindx2); intersect(nodeindx, nodeindx3)];
    dofindx = dofnode(nodeindx,1);
    fl(dofindx,1) = h/2;

end

%%
function [fel] = levelsetforce( l_ndofel, ngpt, iel, xgp, wgp, ndim, nnodel, xel, phiel, h, vnel, dt, dtau, rinit )
% creates the force vector for the level set part

    fel = zeros(l_ndofel,1);
    for igp = 1:ngpt
        gpnum = ngpt*(iel - 1) + igp ;
        xi = xgp(igp,1);
        eta = xgp(igp,1);
        % shape functions
        [n, dn] = shapefuncquad( ndim, nnodel, xi, eta ) ;
        % jacobian
        J = dn * xel ; detJ = det(J) ; dnxy = J\dn;
        % values from previous iteration
        phi = n*phiel ; dphi = dnxy*phiel ; normdphi = norm(dphi) ;
        % normal vector
        normal = -dphi / normdphi ;
        % s function
        s = phi / sqrt(phi^2 + h^2 * normdphi^2) ;
        % w function
        w = s * (dphi / normdphi) ;
        % velocity vector
        vn = n * vnel;
        v = vn * normal ;
        % beta values
        beta1 = 0;
        beta2 = 0;
        if norm(v) > 0
            beta1 = 1/( 2 * sqrt( dt^(-2) + norm(J\v)^2 ) ) ;
        end
        if norm(w) > 0
            beta2 = 1/( 2 * sqrt( dtau^(-2) + norm(J\w)^2 ) ) ;
        end
        % W shape functions
        W = n + beta1 * transpose(v) * dnxy ;
        Wt = n + beta2 * transpose(w) * dnxy ;
        
        r3 = transpose(n) * phi * detJ * wgp(igp,1) * wgp(igp,2) ;
        if rinit == 0
            r1 = transpose(W) * vn * normdphi * detJ * wgp(igp,1) * wgp(igp,2) ;
            
            fel = fel + dt*r1 + r3;
        else
            r1 = transpose(Wt) * s * detJ * wgp(igp,1) * wgp(igp,2) ;
            r2 = transpose(Wt) * (transpose(w) * dphi) * detJ * wgp(igp,1) * wgp(igp,2) ;
            
            fel = fel + dtau*(r1 - r2) + r3;
        end 
    end

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
