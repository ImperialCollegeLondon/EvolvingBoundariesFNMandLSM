%% TOPOPT
%
%   VELOCITY CALCULATION

function [vn, vol, obj, lmult, rpen, elU, ndU] = topt_velocity( s_mesh, s_geo, s_param, s_sol )

    nel = s_mesh.nel;
    u = s_sol.u;
    dofel = s_mesh.dofel;
    rho = s_param.rho;
    nnode = s_mesh.nnode;
    con = s_mesh.con;
    Ael = s_param.Ael;
    lmult = s_param.lmult;
    rpen = s_param.rpen;
    Lx = s_geo.Lx;
    Ly = s_geo.Ly;
    f = s_sol.fe;
    h = s_param.h;
    l_u = s_sol.l_u;
    divx = s_mesh.divx;
    divy = s_mesh.divy;
    ndofel = s_mesh.ndofel;
    thck = s_geo.thck;
    ngpt = s_mesh.ngpt;
    xgp = s_mesh.xgp;
    wgp = s_mesh.wgp;
    ndim = s_mesh.ndim;
    nnodel = s_mesh.nnodel;
    nstr = s_mesh.nstr;
    x = s_mesh.x;
    E = s_geo.E;
    nu = s_geo.nu;
    l_ndofel = s_mesh.l_ndofel;
    l_dofel = s_mesh.l_dofel;
    dt = s_param.dt;
    ndof = s_mesh.ndof;
    
    
    % elemental strain energy
    elU = zeros(nel,1);
    [D] = constitutivematrix( E, nu );
    xel = x(con(1,:), :) ;
    kel = kquad( ndofel, thck, ngpt, xgp, wgp, ndim, nnodel, xel, D, nstr ) ;
    for iel = 1:nel
        uel = u(dofel(iel, :), 1);
        elU(iel,1) = rho(iel,1) * transpose(uel) * kel * uel;
    end

    
    % nodal strain energy density
    ndU = zeros(nnode,1);
    for ind = 1:nnode
        els = [find(con(:,1) == ind); find(con(:,2) == ind); find(con(:,3) == ind); find(con(:,4) == ind)] ;
        ndU(ind,1) = sum(elU(els,1))/size(els,1);
    end

    
    
    % current volume
    vol = sum(rho(:,1)) * Ael ;
    
    
    
    % velocity field
    vfrac = 0.5 ;
    vmax = vfrac*Lx*Ly;
    lmult = lmult + (1/rpen) * (vol - vmax);
    if vol - vmax > 0.1
        rpen = 0.9*rpen;
    end
    lamb = max([0, lmult + (1/rpen) * (vol - vmax)]);
    vn = (ndU/Ael - lamb);
    
    
    
    % regularize velocity field
    vr = regvel( h, l_u, divx, divy, Lx, Ly, nnode, vn ) ;
    for i = 1:5
        f = zeros(nnode, 1);
        rows = zeros(ceil(0.001*ndof*ndof), 1) ;
        cols = rows ;
        vals = rows ;
        count = 1;
        for iel = 1:nel
            xel = x(con(iel,:), :) ;
            vel = vr(con(iel,:), :) ;
            phiel = l_u(l_dofel(iel,:),:) ;
            glbdfs = l_dofel(iel,:)';
            %stifness matrix
            kel = l_kquad( l_ndofel, ngpt, xgp, wgp, ndim, nnodel, xel) ;
            [i,j,v] = find(kel);
            rows(count:(count - 1 + size(i,1)), 1) = glbdfs(i,1);
            cols(count:(count - 1 + size(j,1)), 1) = glbdfs(j,1);
            vals(count:(count - 1 + size(v,1)), 1) = v ;
            count = count + size(i,1);
            %force vector
            fel = levelsetforce( l_ndofel, ngpt, xgp, wgp, ndim, nnodel, xel, phiel, vel, h, dt ) ;
            f(l_dofel(iel,:),1) = f(l_dofel(iel,:),1) + fel ;
        end
        rows(count:end,:) = [];
        cols(count:end,:) = [];
        vals(count:end,:) = [];
        k = sparse(rows,cols,vals);
        
        %new vel
        vr = k\f;
    end   
    vn = vr;
    
    

    % objective function
    obj = sum(elU(:,1)) ;
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
function [b] = bmatrixquad( nnodel, dnxy, nstr, ndofel )
% calculates the b matrix for a quadrilateral element
    
    b = zeros(nstr, ndofel) ;
    for ind = 1:nnodel
        b(:, (2*ind - 1):(2*ind) ) = [ dnxy(1,ind),           0 ;
                                                 0, dnxy(2,ind) ;
                                       dnxy(2,ind), dnxy(1,ind)]; 
    end

end

%%
function [D] = constitutivematrix( E, nu )
% constructs the plane stress constitutive matrix
 
    c = E/(1 - nu^2);
    D = c*[ 1, nu,          0 ;
           nu,  1,          0 ;
            0,  0, 0.5*(1-nu)];

end

%%
function [kel] = kquad( ndofel, thck, ngpt, xgp, wgp, ndim, nnodel, xel, D, nstr )
% constructs the elemental stiffness matrix
% for a quadrilateral isoparametric element

    kel = sparse(ndofel,ndofel);
    for igp = 1:ngpt
        xi = xgp(igp,1);
        eta = xgp(igp,1);
        % shape functions
        [n, dn] = shapefuncquad( ndim, nnodel, xi, eta ) ;
        % jacobian
        J = dn * xel ; detJ = det(J) ;
        % b matrix
        dnxy = J\dn ;
        [b] = bmatrixquad( nnodel, dnxy, nstr, ndofel ) ;
        
        % mechanical part of the stiffness matrix
        kmech = thck * transpose(b) * D * b * detJ * wgp(igp,1) * wgp(igp,2) ;
        
        % atribution
        kel = kel + kmech;
    end

end

%%
function [kel] = l_kquad( l_ndofel, ngpt, xgp, wgp, ndim, nnodel, xel)
% constructs the elemental stiffness matrix
% for a quadrilateral isoparametric element

    kel = sparse(l_ndofel, l_ndofel);
    for igp = 1:ngpt
        xi = xgp(igp,1);
        eta = xgp(igp,1);
        % shape functions
        [n, dn] = shapefuncquad( ndim, nnodel, xi, eta ) ;
        % jacobian
        J = dn * xel ; detJ = det(J) ;
        
        % level set part of the stiffness matrix 
        kphi = transpose(n) * n * detJ * wgp(igp,1) * wgp(igp,2) ;
        kphi = eye(nnodel) .* sum(kphi,2) ;
        
        kel = kel + kphi  ;
    end

end

%%
function [fel] = levelsetforce( l_ndofel, ngpt, xgp, wgp, ndim, nnodel, xel, phiel, vel, h, dt )
% creates the force vector for the level set part

    fel = sparse(l_ndofel,1);
    for igp = 1:ngpt
        xi = xgp(igp,1);
        eta = xgp(igp,1);
        % shape functions
        [n, dn] = shapefuncquad( ndim, nnodel, xi, eta ) ;
        % jacobian
        J = dn * xel ; detJ = det(J) ; dnxy = J\dn;
        % values from previous iteration
        phi = n*phiel ; dphi = dnxy*phiel ; normdphi = norm(dphi) ;
        vn = n*vel; dvn = dnxy*vel;
        % normal vector
        normal = -dphi / normdphi ;
        % s function
        s = phi / sqrt(phi^2 + h^2 * normdphi^2) ;
        % velocity vector
        v = -s * normal ;
        % beta values
        beta1 = 0;
        if norm(v) > 0
            beta1 = 1/( 2 * sqrt( dt^(-2) + norm(J\v)^2 ) ) ;
        end
        % W shape functions
        W = n + beta1 * transpose(v) * dnxy ;
        
        r3 = transpose(n) * vn * detJ * wgp(igp,1) * wgp(igp,2) ;
        r1 = transpose(W) * (transpose(v) * dvn) * detJ * wgp(igp,1) * wgp(igp,2) ;
        
        fel = fel - dt*r1 + r3;
    end

end

%%
function [v] = regvel( h, l_u, divx, divy, Lx, Ly, nnode, vn )

    eps = h/20;
    
    phi = l_u;
    
    s = phi ./ sqrt( phi.^2 + eps^2 );
    
    
    % finite difference shifts
    nodx = divx + 1;
    nody = divy + 1;
    smat = reshape(s,[nodx, nody]);
    smat = transpose(smat);
    sr = 0*smat; sl = 0*smat; su = 0*smat; sd = 0*smat;
    % right w
    sr(:, 1:(end-1)) = smat(:, 2:end) ;
    sr(:, end) = 2*smat(:,end) - smat(:,end-1) ;
    sr = transpose(sr);
    sr = reshape(sr, [nnode,1]);
    % left e
    sl(:, 2:end) = smat(:, 1:(end-1)) ;
    sl(:, 1) = 2*smat(:,1) - smat(:,2) ;
    sl = transpose(sl);
    sl = reshape(sl, [nnode,1]);
    % up s
    su(2:end, :) = smat(1:(end-1), :) ;
    su(1, :) = 2*smat(1,:) - smat(2,:) ;
    su = transpose(su);
    su = reshape(su, [nnode,1]);
    % down n
    sd(1:(end-1), :) = smat(2:end, :) ;
    sd(end, :) = 2*smat(end,:) - smat(end-1,:) ;
    sd = transpose(sd);
    sd = reshape(sd, [nnode,1]);
    
    dx = Lx/divx;
    dy = Ly/divy;
    dsx = (sl - sr)/(2*dx);
    dsy = (sd - su)/(2*dy);
    
    delta = 0.5 * sqrt( dsx.^2 + dsy.^2 );
    
    b = -vn .* delta ;
    
    v = -b;

end





