% **************************************************************************
%  MODELLING EVOLVING BOUNDARIES USING FNM AND LSM (April 2020)
% **************************************************************************
%
%
%   ls_velocity.m
%
%
%   FUNCTIONS:
%       - ls_velocity - computes the velocity field, objective function and
%                       volume of the domain
%       - lsforcequad4 - modified LSFEM force vector for the velocity
%                        extension. Computes for quadrilateral elements
%
%
% **************************************************************************
%                   developed by Rui O.S.S. da Costa and Silvestre T. Pinho
%                                      Imperial College London, Aeronautics
% **************************************************************************

function [vn, obj, vol, lmult, pen] = ls_velocity( s_mesh, s_geo, s_param, s_sol )

    % **************************************************************************
    %
    %   INPUTS:
    %       s_mesh - structure arry with mesh data
    %       s_geo - structure arry with geometric data
    %       s_param - structure arry with generic parameters
    %       s_sol - structure array with solution data
    %
    %   OUTPUTS:
    %       vn - nodal sensitivity/velocity field for the mean compliance
    %            objective function
    %       obj - mean compliance value
    %       vol - volume of the domain
    %       lmult - lagrange multiplier for the volume constraint
    %       pen - penalty parameter for the volume constraint
    %
    % **************************************************************************

    % VARIABLE INITIALISATION
        nel = s_mesh.nel ;
        ndofnode = s_mesh.ndofnode ;
        x = s_mesh.x;
        con = s_mesh.con ;
        thck = s_geo.thck ;
        ngptquad = s_mesh.ngptquad;
        xgpquad = s_mesh.xgpquad;
        wgpquad = s_mesh.wgpquad;
        ndim = s_mesh.ndim;
        D = s_geo.D ;
        nstr = s_mesh.nstr ;
        u = s_sol.u ;
        dofel = s_mesh.dofel ;
        ndofel = s_mesh.ndofel ;
        nfloatdofel = s_mesh.nfloatdofel ;
        nnodel = s_mesh.nnodel ;
        nfloatnodel = s_mesh.nfloatnodel ;
        nnode = s_mesh.nnode ;
        nfloatnode = s_mesh.nfloatnode ;
        Ael = s_geo.Ael ;
        phi = s_param.p ;
        elstat = s_mesh.elstat ;
        lmult = s_param.lmult ;
        pen = s_param.pen ;
        vfrac = s_param.vfrac ;
        Lx = s_geo.Lx ;
        Ly = s_geo.Ly ;
        h = s_param.h;
        c = s_param.c;
        eps = s_sol.eps;
        sig = s_sol.sig;
        elsin = s_mesh.elsin ;
        elsbnd = s_mesh.elsbnd ;
        fbnd = s_param.fbnd;
        elsoutdom = s_mesh.elsoutdom;
        ndsoutdom = s_mesh.ndsoutdom;
        elsfixbnd = s_mesh.elsfixbnd;
        ubc = s_param.ubc;
        ndindx = [];

    
    % MODIFY BOUNDARY NODES VECTOR FOR DIFFERENT TEST CASES
        if ubc == 3 % L bracket
            xs = ((0+h/2):h:Lx)';
            ys = ((0+h/2):h:Ly)';
            dx = abs(xs - (2*Lx/5 + h/100));
            xs = xs(find(dx == min(dx)),1);
            dy = abs(ys - (2*Ly/5 + h/100));
            ys = ys(find(dy == min(dy)),1);

            ndindx = find( ( x(:,1) > xs - h & x(:,1) < xs + h & x(:,2) < Ly + h & x(:,2) > ys - h ) | ...
                           ( x(:,2) > ys - h & x(:,2) < ys + h & x(:,1) < Lx + h & x(:,1) > xs - h ) );
        elseif ubc == 4 % bracket with 2 holes
            ndindaux = con(elsfixbnd,:);
            ndindx = reshape(ndindaux, [size(ndindaux,1)*size(ndindaux,2), 1]);
            ndindx = unique(ndindx);
        end
        fbnd = setdiff(fbnd,ndindx);
    
    
    % NODAL MEAN COMPLIANCE
        ndU = zeros(nnode,1);
        vol = 0;
        obj = 0;
        for ind = fbnd'
            ndU(ind,1) = sig(ind,:) * eps(ind,:)' ;
        end
        
    % ELEMENTAL OBJECTIVE AND VOLUME
        % elemental stiffness matrix
            [kel] = fem_kquad4( ndofel - nfloatdofel, thck, ngptquad, xgpquad, wgpquad, ndim, nnodel - nfloatnodel, x(con(elsin(1),1:4),:), D, nstr ) ;
    
        % quadrilaterals fully inside the domain
            for iel = elsin'
                nsub = fnm_partitioninfo(elstat(iel,1)) ;
                for isub = 1:nsub
                    [scon] = fnm_subcon( elstat(iel,1), isub ) ;
                    nnodesub = size(scon,2) ;
                    sdofel = zeros(1, ndofnode*nnodesub) ;
                    for ind = 1:nnodesub
                        for idof = 1:ndofnode
                            sdofel(1, ndofnode * ind + idof - ndofnode) = ndofnode * scon(1, ind) + idof - ndofnode ;
                        end
                    end
                    dofs = dofel(iel, sdofel) ;

                    uel = u(dofs, 1) ;

                    elU = transpose(uel) * kel * uel ;
                    obj = obj + elU ;
                    vol = vol + Ael ;
                end
            end
        

    % LAGRANGE MULTIPLIER CALCULATION
        vmax = vfrac * Lx * Ly;
        lmult = lmult + (1/pen) * (vol - vmax);
        if vol - vmax > 0.1
            pen = 0.9*pen;
        end
        lamb = max([0, lmult + (1/pen)*(vol - vmax)]);
        vn = zeros(nnode, 1) - lamb;
        vn = vn + ndU;
    
    
    % LIMITING THE VELOCITY FIELD TO THE BOUNDARY NODES
        nds = (1:nnode)' ;
        nds = setdiff(nds, fbnd);
        vn(nds,1) = 0;
    
    
    % EXTRAPOLATION OF BOUNDARY VELOCITY TO NEIGHBOURING NODES
    % projecting the velocity of the floating nodes normally to a line
    % containing the real node and extrapolating on that line to the
    % position of the real node
        ndcount = zeros(nnode - nfloatnode,1);
        for iel = elsbnd'

            if elstat(iel,1) == 11

                xel = x(con(iel,:),:) ;
                vel = vn(con(iel,:),1) ;

                vecta = xel(8, :)' - xel(5, :)' ;
                vectb = xel(1, :)' - xel(5, :)' ;
                vectc = xel(2, :)' - xel(5, :)' ;
                vectd = xel(3, :)' - xel(5, :)' ;
                vecte = xel(4, :)' - xel(5, :)' ;

                vecti = vecta / norm(vecta) ;

                x5 = 0 ;
                x1 = x5 + dot(vecti, vectb) ;
                x2 = x5 + dot(vecti, vectc) ;
                x3 = x5 + dot(vecti, vectd) ;
                x4 = x5 + dot(vecti, vecte) ;
                x8 = x5 + norm(vecta) ;

                v1 = ( (vel(5,1) - vel(8,1)) / (x5 - x8) ) * ( x1 - x5 ) + vel(5,1) ;
                v2 = ( (vel(5,1) - vel(8,1)) / (x5 - x8) ) * ( x2 - x5 ) + vel(5,1) ;
                v3 = ( (vel(5,1) - vel(8,1)) / (x5 - x8) ) * ( x3 - x5 ) + vel(5,1) ;
                v4 = ( (vel(5,1) - vel(8,1)) / (x5 - x8) ) * ( x4 - x5 ) + vel(5,1) ;

                vn(con(iel,1:4), 1) = vn(con(iel,1:4), 1) + [v1; v2; v3; v4];
                ndcount(con(iel,1:4), 1) = ndcount(con(iel,1:4), 1) + 1 ;

            elseif elstat(iel,1) == 12

                xel = x(con(iel,:),:) ;
                vel = vn(con(iel,:),1) ;

                vecta = xel(8, :)' - xel(7, :)' ;
                vectb = xel(1, :)' - xel(7, :)' ;
                vectc = xel(2, :)' - xel(7, :)' ;
                vectd = xel(3, :)' - xel(7, :)' ;
                vecte = xel(4, :)' - xel(7, :)' ;

                vecti = vecta / norm(vecta) ;

                x7 = 0 ;
                x1 = x7 + dot(vecti, vectb) ;
                x2 = x7 + dot(vecti, vectc) ;
                x3 = x7 + dot(vecti, vectd) ;
                x4 = x7 + dot(vecti, vecte) ;
                x8 = x7 + norm(vecta) ;

                v1 = ( (vel(7,1) - vel(8,1)) / (x7 - x8) ) * ( x1 - x7 ) + vel(7,1) ;
                v2 = ( (vel(7,1) - vel(8,1)) / (x7 - x8) ) * ( x2 - x7 ) + vel(7,1) ;
                v3 = ( (vel(7,1) - vel(8,1)) / (x7 - x8) ) * ( x3 - x7 ) + vel(7,1) ;
                v4 = ( (vel(7,1) - vel(8,1)) / (x7 - x8) ) * ( x4 - x7 ) + vel(7,1) ;

                vn(con(iel,1:4), 1) = vn(con(iel,1:4), 1) + [v1; v2; v3; v4];
                ndcount(con(iel,1:4), 1) = ndcount(con(iel,1:4), 1) + 1 ;

            elseif elstat(iel,1) == 13

                xel = x(con(iel,:),:) ;
                vel = vn(con(iel,:),1) ;

                vecta = xel(7, :)' - xel(6, :)' ;
                vectb = xel(1, :)' - xel(6, :)' ;
                vectc = xel(2, :)' - xel(6, :)' ;
                vectd = xel(3, :)' - xel(6, :)' ;
                vecte = xel(4, :)' - xel(6, :)' ;

                vecti = vecta / norm(vecta) ;

                x6 = 0 ;
                x1 = x6 + dot(vecti, vectb) ;
                x2 = x6 + dot(vecti, vectc) ;
                x3 = x6 + dot(vecti, vectd) ;
                x4 = x6 + dot(vecti, vecte) ;
                x7 = x6 + norm(vecta) ;

                v1 = ( (vel(6,1) - vel(7,1)) / (x6 - x7) ) * ( x1 - x6 ) + vel(6,1) ;
                v2 = ( (vel(6,1) - vel(7,1)) / (x6 - x7) ) * ( x2 - x6 ) + vel(6,1) ;
                v3 = ( (vel(6,1) - vel(7,1)) / (x6 - x7) ) * ( x3 - x6 ) + vel(6,1) ;
                v4 = ( (vel(6,1) - vel(7,1)) / (x6 - x7) ) * ( x4 - x6 ) + vel(6,1) ;

                vn(con(iel,1:4), 1) = vn(con(iel,1:4), 1) + [v1; v2; v3; v4];
                ndcount(con(iel,1:4), 1) = ndcount(con(iel,1:4), 1) + 1 ;

            elseif elstat(iel,1) == 14

                xel = x(con(iel,:),:) ;
                vel = vn(con(iel,:),1) ;

                vecta = xel(6, :)' - xel(5, :)' ;
                vectb = xel(1, :)' - xel(5, :)' ;
                vectc = xel(2, :)' - xel(5, :)' ;
                vectd = xel(3, :)' - xel(5, :)' ;
                vecte = xel(4, :)' - xel(5, :)' ;

                vecti = vecta / norm(vecta) ;

                x5 = 0 ;
                x1 = x5 + dot(vecti, vectb) ;
                x2 = x5 + dot(vecti, vectc) ;
                x3 = x5 + dot(vecti, vectd) ;
                x4 = x5 + dot(vecti, vecte) ;
                x6 = x5 + norm(vecta) ;

                v1 = ( (vel(5,1) - vel(6,1)) / (x5 - x6) ) * ( x1 - x5 ) + vel(5,1) ;
                v2 = ( (vel(5,1) - vel(6,1)) / (x5 - x6) ) * ( x2 - x5 ) + vel(5,1) ;
                v3 = ( (vel(5,1) - vel(6,1)) / (x5 - x6) ) * ( x3 - x5 ) + vel(5,1) ;
                v4 = ( (vel(5,1) - vel(6,1)) / (x5 - x6) ) * ( x4 - x5 ) + vel(5,1) ;

                vn(con(iel,1:4), 1) = vn(con(iel,1:4), 1) + [v1; v2; v3; v4];
                ndcount(con(iel,1:4), 1) = ndcount(con(iel,1:4), 1) + 1 ;

            elseif elstat(iel,1) == 21

                xel = x(con(iel,:),:) ;
                vel = vn(con(iel,:),1) ;

                vecta = xel(7, :)' - xel(5, :)' ;
                vectb = xel(1, :)' - xel(5, :)' ;
                vectc = xel(2, :)' - xel(5, :)' ;
                vectd = xel(3, :)' - xel(5, :)' ;
                vecte = xel(4, :)' - xel(5, :)' ;

                vecti = vecta / norm(vecta) ;

                x5 = 0 ;
                x1 = x5 + dot(vecti, vectb) ;
                x2 = x5 + dot(vecti, vectc) ;
                x3 = x5 + dot(vecti, vectd) ;
                x4 = x5 + dot(vecti, vecte) ;
                x7 = x5 + norm(vecta) ;

                v1 = ( (vel(5,1) - vel(7,1)) / (x5 - x7) ) * ( x1 - x5 ) + vel(5,1) ;
                v2 = ( (vel(5,1) - vel(7,1)) / (x5 - x7) ) * ( x2 - x5 ) + vel(5,1) ;
                v3 = ( (vel(5,1) - vel(7,1)) / (x5 - x7) ) * ( x3 - x5 ) + vel(5,1) ;
                v4 = ( (vel(5,1) - vel(7,1)) / (x5 - x7) ) * ( x4 - x5 ) + vel(5,1) ;

                vn(con(iel,1:4), 1) = vn(con(iel,1:4), 1) + [v1; v2; v3; v4];
                ndcount(con(iel,1:4), 1) = ndcount(con(iel,1:4), 1) + 1 ;

            elseif elstat(iel,1) == 22

                xel = x(con(iel,:),:) ;
                vel = vn(con(iel,:),1) ;

                vecta = xel(8, :)' - xel(6, :)' ;
                vectb = xel(1, :)' - xel(6, :)' ;
                vectc = xel(2, :)' - xel(6, :)' ;
                vectd = xel(3, :)' - xel(6, :)' ;
                vecte = xel(4, :)' - xel(6, :)' ;

                vecti = vecta / norm(vecta) ;

                x6 = 0 ;
                x1 = x6 + dot(vecti, vectb) ;
                x2 = x6 + dot(vecti, vectc) ;
                x3 = x6 + dot(vecti, vectd) ;
                x4 = x6 + dot(vecti, vecte) ;
                x8 = x6 + norm(vecta) ;

                v1 = ( (vel(6,1) - vel(8,1)) / (x6 - x8) ) * ( x1 - x6 ) + vel(6,1) ;
                v2 = ( (vel(6,1) - vel(8,1)) / (x6 - x8) ) * ( x2 - x6 ) + vel(6,1) ;
                v3 = ( (vel(6,1) - vel(8,1)) / (x6 - x8) ) * ( x3 - x6 ) + vel(6,1) ;
                v4 = ( (vel(6,1) - vel(8,1)) / (x6 - x8) ) * ( x4 - x6 ) + vel(6,1) ;

                vn(con(iel,1:4), 1) = vn(con(iel,1:4), 1) + [v1; v2; v3; v4];
                ndcount(con(iel,1:4), 1) = ndcount(con(iel,1:4), 1) + 1 ;

            elseif elstat(iel,1) == 31

                xel = x(con(iel,:),:) ;
                vel = vn(con(iel,:),1) ;

                vecta = xel(8, :)' - xel(5, :)' ;
                vectb = xel(1, :)' - xel(5, :)' ;
                vectc = xel(2, :)' - xel(5, :)' ;
                vectd = xel(4, :)' - xel(5, :)' ;
                vecte = xel(7, :)' - xel(6, :)' ;
                vectf = xel(2, :)' - xel(6, :)' ;
                vectg = xel(3, :)' - xel(6, :)' ;
                vecth = xel(4, :)' - xel(6, :)' ;

                vecti1 = vecta / norm(vecta) ;
                vecti2 = vecte / norm(vecte) ;

                x5 = 0 ;
                x1 = x5 + dot(vecti1, vectb) ;
                x21 = x5 + dot(vecti1, vectc) ;
                x41 = x5 + dot(vecti1, vectd) ;
                x8 = x5 + norm(vecta) ;
                x6 = 0 ;
                x22 = x6 + dot(vecti2, vectf) ;
                x3 = x6 + dot(vecti2, vectg) ;
                x42 = x6 + dot(vecti2, vecth) ;
                x7 = x6 + norm(vecte) ;

                v1 = ( (vel(5,1) - vel(8,1)) / (x5 - x8) ) * ( x1 - x5 ) + vel(5,1) ;
                v2 = 0.5*( ( (vel(5,1) - vel(8,1)) / (x5 - x8) ) * ( x21 - x5 ) + vel(5,1) ) + 0.5*( ( (vel(6,1) - vel(7,1)) / (x6 - x7) ) * ( x22 - x6 ) + vel(6,1) ) ;
                v3 = ( (vel(6,1) - vel(7,1)) / (x6 - x7) ) * ( x3 - x6 ) + vel(6,1) ;
                v4 = 0.5*( ( (vel(5,1) - vel(8,1)) / (x5 - x8) ) * ( x41 - x5 ) + vel(5,1) ) + 0.5*( ( (vel(6,1) - vel(7,1)) / (x6 - x7) ) * ( x42 - x6 ) + vel(6,1) ) ;

                vn(con(iel,1:4), 1) = vn(con(iel,1:4), 1) + [v1; v2; v3; v4];
                ndcount(con(iel,1:4), 1) = ndcount(con(iel,1:4), 1) + 1 ;

            elseif elstat(iel,1) == 32

                xel = x(con(iel,:),:) ;
                vel = vn(con(iel,:),1) ;

                vecta = xel(6, :)' - xel(5, :)' ;
                vectb = xel(1, :)' - xel(5, :)' ;
                vectc = xel(2, :)' - xel(5, :)' ;
                vectd = xel(3, :)' - xel(5, :)' ;
                vecte = xel(8, :)' - xel(7, :)' ;
                vectf = xel(1, :)' - xel(7, :)' ;
                vectg = xel(3, :)' - xel(7, :)' ;
                vecth = xel(4, :)' - xel(7, :)' ;

                vecti1 = vecta / norm(vecta) ;
                vecti2 = vecte / norm(vecte) ;

                x5 = 0 ;
                x11 = x5 + dot(vecti1, vectb) ;
                x2 = x5 + dot(vecti1, vectc) ;
                x31 = x5 + dot(vecti1, vectd) ;
                x6 = x5 + norm(vecta) ;
                x7 = 0 ;
                x12 = x7 + dot(vecti2, vectf) ;
                x32 = x7 + dot(vecti2, vectg) ;
                x4 = x7 + dot(vecti2, vecth) ;
                x8 = x7 + norm(vecte) ;

                v1 = 0.5*( ( (vel(5,1) - vel(6,1)) / (x5 - x6) ) * ( x11 - x5 ) + vel(5,1) ) + 0.5*( ( (vel(7,1) - vel(8,1)) / (x7 - x8) ) * ( x12 - x7 ) + vel(7,1) ) ;
                v2 = ( (vel(5,1) - vel(6,1)) / (x5 - x6) ) * ( x2 - x5 ) + vel(5,1) ;
                v3 = 0.5*( ( (vel(5,1) - vel(6,1)) / (x5 - x6) ) * ( x31 - x5 ) + vel(5,1) ) + 0.5*( ( (vel(7,1) - vel(8,1)) / (x7 - x8) ) * ( x32 - x7 ) + vel(7,1) ) ;
                v4 = ( (vel(7,1) - vel(8,1)) / (x7 - x8) ) * ( x4 - x7 ) + vel(7,1) ;

                vn(con(iel,1:4), 1) = vn(con(iel,1:4), 1) + [v1; v2; v3; v4];
                ndcount(con(iel,1:4), 1) = ndcount(con(iel,1:4), 1) + 1 ;

            end

        end
        nds = find(ndcount);
        vn(nds,1) = vn(nds,1) ./ ndcount(nds,1) ;

    
    % EXTENSION OF THE VELOCITY FIELD TO THE ENTIRE DOMAIN
        dt = 1 * (h / max(abs(vn(:,1))) ) ;
        vr = vn ;
        els = setdiff((1:nel)',elsoutdom);
        for i = 1:5
            rows = zeros(ceil(0.01*nnode*nnode), 1) ;
            cols = rows ;
            vals = rows ;
            count = 1;
            f = zeros(nnode - nfloatnode, 1);
            for iel = els'
                if iel == 3094
                    pause(.1)
                end

                nnodesub = nnodel - nfloatnodel;
                nds = con(iel, 1:nnodesub) ;
                xel = x(nds,:) ;
                vel = vr(nds, :) ;
                phiel = phi(nds,:) ;
                %stifness matrix
                kel = ls_kquad4( nnodesub, ngptquad, xgpquad, wgpquad, ndim, nnodesub, xel, 0, c, h ) ;
                [irow,jcol,vval] = find(kel);
                rows(count:(count - 1 + size(irow,1)), 1) = nds(irow)';
                cols(count:(count - 1 + size(jcol,1)), 1) = nds(jcol)';
                vals(count:(count - 1 + size(vval,1)), 1) = vval ;
                count = count + size(irow,1);
                %force vector
                fel = lsforcequad4( ndim, nnodesub, ngptquad, xgpquad, wgpquad, nnodesub, xel, phiel, h, vel, dt ) ;
                f(nds,1) = f(nds,1) + fel ;

            end
            rows(count:end,:) = [];
            cols(count:end,:) = [];
            vals(count:end,:) = [];
            k = sparse(rows,cols,vals) ;

            %new vel
            ands = setdiff((1:(nnode - nfloatnode))',ndsoutdom);
            vr(ands,1) = k(ands, ands)\f(ands, 1);
        end   
        vn = vr;



    % NORMALISING THE VELOCITY FIELD
    mx = max(abs(vn(:,1)));
    vn = vn ./ mx ;
    
end


function [fel] = lsforcequad4( ndim, ndofel, ngpt, xgp, wgp, nnodel, xel, phiel, h, vel, dt )

    % **************************************************************************
    %
    %   INPUTS:
    %       ndim - number of spatial dimensions
    %       ndofel - number of total dof per element
    %       ngpt -  number of integration points
    %       xgp - matrix of coordinates for the integration points
    %       wgp - matrix of weights for the integration points
    %       nnodel - number of total nodes per element
    %       xel - matrix of nodal coordinates
    %       phiel - vector of nodal LS values
    %       h - minimum mesh size
    %       vel - vector of nodal velocities
    %       dt - time step
    %
    %
    %   OUTPUTS:
    %       fel - elemental force vector for velocity extension
    %
    % **************************************************************************

    % VARIABLE INITIALISATION
        fel = zeros(ndofel,1);
        
    % LOOP THROUGH INTEGRATION POINTS
        for igp = 1:ngpt
            
            % coordinates of integration points
                xi = xgp(igp,1);
                eta = xgp(igp,2);
                
            % shape functions
                [n, dn] = fem_sfquad4( ndim, nnodel, xi, eta ) ;
            
            % jacobian
                J = dn * xel ; detJ = det(J) ; dnxy = J\dn ;
            
            % values from previous iteration
                phi = n*phiel ; dphi = dnxy*phiel ; normdphi = norm(dphi) ; 
                vn = n*vel ; dvn = dnxy*vel ;
            
            % normal vector
                normal = -dphi / normdphi ;
            
            % s function
                s = phi / sqrt(phi^2 + h^2 * normdphi^2) ;
            
            % w function
                w = s * normal ;
            
            % velocity vector
                v = -w ;

            % elemental vector computation
                r3 = transpose(n) * vn * detJ * wgp(igp) ;
                r1 = transpose(n) * ( transpose(v) * dvn ) * detJ * wgp(igp) ;
                fel = fel - dt*r1 + r3;
        end
end
