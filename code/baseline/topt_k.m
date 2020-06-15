%% TOPOPT
%
%   CALCULATION OF THE STIFFNESS MATRIX

function [k] = topt_k(mode, s_mesh, s_geo, s_param, s_sol)

    
    nel = s_mesh.nel;
    con = s_mesh.con;
    dofel = s_mesh.dofel ;
    ndof = s_mesh.ndof ;
    ndofel = s_mesh.ndofel;
    thck = s_geo.thck;
    ngpt = s_mesh.ngpt;
    xgp = s_mesh.xgp;
    wgp = s_mesh.wgp;
    ndim = s_mesh.ndim;
    nnodel = s_mesh.nnodel;
    x = s_mesh.x;
    E = s_geo.E;
    nu = s_geo.nu;
    nstr = s_mesh.nstr ;
    c = s_param.c;
    h = s_param.h;
    rinit = s_param.rinit ;
    G = s_param.G;
    rho = s_param.rho;
    l_ndof = s_mesh.l_ndof;
    l_ndofel = s_mesh.l_ndofel;
    l_dofel = s_mesh.l_dofel;
    
    if mode == 1
        
        
        rows = zeros(ceil(0.001*ndof*ndof), 1) ;
        cols = rows ;
        vals = rows ;
        count = 1;
        [D] = constitutivematrix( E, nu ) ;
        for iel = 1:nel
            xel = x(con(iel,:), :) ;
            glbdfs = (dofel(iel,:))' ;
            
            % elemental stiffness matrix
            [kel] = kquad( ndofel, thck, ngpt, xgp, wgp, ndim, nnodel, xel, rho(iel,1)*D, nstr ) ;
            
            
            % assembling
            [i,j,v] = find(kel);
            rows(count:(count - 1 + size(i,1)), 1) = glbdfs(i,1);
            cols(count:(count - 1 + size(j,1)), 1) = glbdfs(j,1);
            vals(count:(count - 1 + size(v,1)), 1) = v ;
            count = count + size(i,1);
            
            
        end
        rows(count:end,:) = [];
        cols(count:end,:) = [];
        vals(count:end,:) = [];
        k = sparse(rows,cols,vals);
        
        
    elseif mode == 2
       
        
        rows = zeros(ceil(0.001*l_ndof*l_ndof), 1) ;
        cols = rows ;
        vals = rows ;
        count = 1;
        for iel = 1:nel
            xel = x(con(iel,:), :) ;
            glbdfs = (l_dofel(iel,:))' ;
            
            % elemental stiffness matrix
            [kel] = l_kquad( l_ndofel, ngpt, xgp, wgp, ndim, nnodel, xel, c, h, rinit ) ;
            
            % assembling
            [i,j,v] = find(kel);
            rows(count:(count - 1 + size(i,1)), 1) = glbdfs(i,1);
            cols(count:(count - 1 + size(j,1)), 1) = glbdfs(j,1);
            vals(count:(count - 1 + size(v,1)), 1) = v ;
            count = count + size(i,1);
            
        end
        rows(count:end,:) = [];
        cols(count:end,:) = [];
        vals(count:end,:) = [];
        k = sparse(rows,cols,vals);
        if rinit == 1
            k = k + 10^5 * transpose(G) * G ;
        end
        
        
        
    end
    
    

end

%%
function [kel] = kquad( ndofel, thck, ngpt, xgp, wgp, ndim, nnodel, xel, D, nstr )
% constructs the elemental stiffness matrix
% for a quadrilateral isoparametric element

    kel = zeros(ndofel,ndofel);
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
function [kel] = l_kquad( l_ndofel, ngpt, xgp, wgp, ndim, nnodel, xel, c, h, rinit )
% constructs the elemental stiffness matrix
% for a quadrilateral isoparametric element

    kel = zeros(l_ndofel, l_ndofel);
    for igp = 1:ngpt
        xi = xgp(igp,1);
        eta = xgp(igp,1);
        % shape functions
        [n, dn] = shapefuncquad( ndim, nnodel, xi, eta ) ;
        % jacobian
        J = dn * xel ; detJ = det(J) ;
        % b matrix
        dnxy = J\dn ;
        
        % level set part of the stiffness matrix 
        if rinit == 1
            kphia = transpose(n) * n * detJ * wgp(igp,1) * wgp(igp,2) ;
            kphib = transpose(dnxy) * (dnxy) * detJ * wgp(igp,1) * wgp(igp,2) ;
            kphi = kphia + c * h^2 * kphib ;
        else
            kphi = transpose(n) * n * detJ * wgp(igp,1) * wgp(igp,2) ;
            kphi = eye(nnodel) .* sum(kphi,2) ;
        end
        
        kel = kel + kphi  ;
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










