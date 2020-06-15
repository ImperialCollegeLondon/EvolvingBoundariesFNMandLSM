% **************************************************************************
%  MODELLING EVOLVING BOUNDARIES USING FNM AND LSM (April 2020)
% **************************************************************************
%
%
%   fem_post.m
%
%
%   FUNCTIONS:
%       - fem_post - computes the stress and strain fields
%       - bmatrixquad - computes the b matrix for quadrilateral elements
%       - bmatrixtri - computes the b matrix for triangular elements
%
%
% **************************************************************************
%                   developed by Rui O.S.S. da Costa and Silvestre T. Pinho
%                                      Imperial College London, Aeronautics
% **************************************************************************


function [sig, eps] = fem_post( s_mesh, s_geo, s_param, s_sol)
    
    % **************************************************************************
    %
    %   INPUTS:
    %       s_mesh - structure arry with mesh data
    %       s_geo - structure arry with geometric data
    %       s_param - structure arry with generic parameters
    %       s_sol - structure array with solution data
    %
    %   OUTPUTS:
    %       sig - stress field matrix
    %       eps - strain field matrix
    %
    % **************************************************************************    

    % VARIABLE INITIALISATION
        D = s_geo.D ;
        nstr = s_mesh.nstr ;
        dofel = s_mesh.dofel ;
        nnode = s_mesh.nnode ;
        x = s_mesh.x;
        con = s_mesh.con;
        ndim = s_mesh.ndim ;
        ngptquad = s_mesh.ngptquad ;
        xgpquad = s_mesh.xgpquad ;
        ngpttri = s_mesh.ngpttri ;
        xgptri = s_mesh.xgptri ;
        u = s_sol.u ;
        elstat = s_mesh.elstat ;
        ndofnode = s_mesh.ndofnode ;
        p = s_param.p;
        sig = zeros(nnode, nstr);
        eps = zeros(nnode, nstr);
        ndcount = zeros(nnode, 1) ;
        elsin = s_mesh.elsin;
        elsbnd = s_mesh.elsbnd;
    
    % LOOP THROGH UN-PARTITIONED ELEMENTS INSIDE THE DOMAIN
        for iel = elsin'
            
            % local indices 
                scon = fnm_subcon( 0, 1 ) ;
                nnodesub = size(scon,2) ;
                ndofsub = nnodesub * ndofnode ;
                sdof = zeros(1, ndofsub);
                for ind = 1:nnodesub
                    for idof = 1:ndofnode
                        sdof(1, ndofnode * ind + idof - ndofnode) = ndofnode * scon(ind) + idof - ndofnode ;
                    end
                end

            % coordinates and displacement
                xel = x( con(iel, scon), :) ;
                uel = u( dofel(iel, sdof), :) ;

            % calculation of sig and eps at the integration points
                siggp = zeros(ngptquad, nstr);
                epsgp = zeros(ngptquad, nstr);
                for igp = 1:ngptquad
                    xi = xgpquad(igp,1);
                    eta = xgpquad(igp,2);
                    [n, dn] = fem_sfquad4( ndim, nnodesub, xi, eta ) ;
                    J = dn * xel ; dnxy = J\dn ;
                    b = bmatrixquad( nnodesub, dnxy, nstr, ndofsub) ;

                    epsgp(igp,:) = (b * uel)' ;
                    siggp(igp,:) = (D * epsgp(igp,:)' )' ;
                end

            % extrapolation to the nodes
                for ind = 1:nnodesub
                    xi = sqrt(3) * sign( xgpquad(ind,1) );
                    eta = sqrt(3) * sign( xgpquad(ind,2) );

                    [n, dn] = fem_sfquad4( ndim, nnodesub, xi, eta ) ;

                    for ist = 1:nstr
                        eps(con(iel,scon(ind)), ist) = eps(con(iel,scon(ind)), ist) + n * epsgp(:, ist) ;
                        sig(con(iel,scon(ind)), ist) = sig(con(iel,scon(ind)), ist) + n * siggp(:, ist) ;
                    end
                end

            % update nodal count for averaging
                ndcount(con(iel,scon), 1) = ndcount(con(iel,scon), 1) + 1 ;  
        end
    
        
    % LOOP THROUGH PARTITIONED ELEMENTS
        for iel = elsbnd'
            
            % loop through sub elements
                nsub = fnm_partitioninfo(elstat(iel,1)) ;
                for isub = 1:nsub

                    % local indices
                        scon = fnm_subcon( elstat(iel,1), isub ) ;
                        nnodesub = size(scon,2) ;
                        ndofsub = nnodesub * ndofnode ;
                        sdof = zeros(1, ndofsub);
                        for ind = 1:nnodesub
                            for idof = 1:ndofnode
                                sdof(1, ndofnode * ind + idof - ndofnode) = ndofnode * scon(ind) + idof - ndofnode ;
                            end
                        end

                    % coordinates, LS values and displacement
                        xel = x( con(iel, scon), :) ;
                        phiel = p( con(iel,scon), :) ;
                        uel = u( dofel(iel, sdof), :) ;
                        rho = ls_density(phiel); % pseudo density based on LS values

                    % compute only for solid sub-elements
                        if rho == 1

                            % quadrilateral sub-elements
                                if nnodesub == 4

                                    % calculation of sig and eps at the integration points
                                        siggp = zeros(ngptquad, nstr);
                                        epsgp = zeros(ngptquad, nstr);
                                        for igp = 1:ngptquad
                                            xi = xgpquad(igp,1);
                                            eta = xgpquad(igp,2);
                                            [n, dn] = fem_sfquad4( ndim, nnodesub, xi, eta ) ;
                                            J = dn * xel ; dnxy = J\dn ;
                                            b = bmatrixquad( nnodesub, dnxy, nstr, ndofsub) ;

                                            epsgp(igp,:) = (b * uel)' ;
                                            siggp(igp,:) = (D * epsgp(igp,:)' )' ;
                                        end

                                    % extrapolation to the nodes
                                        for ind = 1:nnodesub
                                            xi = sqrt(3) * sign( xgpquad(ind,1) );
                                            eta = sqrt(3) * sign( xgpquad(ind,2) );
                                            [n, dn] = fem_sfquad4( ndim, nnodesub, xi, eta ) ;

                                            for ist = 1:nstr
                                                eps(con(iel,scon(ind)), ist) = eps(con(iel,scon(ind)), ist) + n * epsgp(:, ist) ;
                                                sig(con(iel,scon(ind)), ist) = sig(con(iel,scon(ind)), ist) + n * siggp(:, ist) ;
                                            end
                                        end

                                elseif nnodesub == 3
                                    
                                    % calculation of sig and eps at the integration points
                                        siggp = zeros(ngpttri, nstr);
                                        epsgp = zeros(ngpttri, nstr);
                                        for igp = 1:ngpttri
                                            xi = xgptri(igp,1);
                                            eta = xgptri(igp,2);
                                            [n, dn] = fem_sftri3( ndim, nnodesub, xi, eta ) ;
                                            J = dn * xel ; dnxy = J\dn ;
                                            b = bmatrixtri( nnodesub, dnxy, nstr, ndofsub) ;

                                            epsgp(igp,:) = (b * uel)' ;
                                            siggp(igp,:) = (D * epsgp(igp,:)' )' ;
                                        end

                                    % extrapolation to the nodes
                                        for ind = 1:nnodesub
                                            xi = 2 * xgptri(ind,1) - 1/3 ;
                                            eta = 2 * xgptri(ind,2) - 1/3 ;
                                            [n, dn] = fem_sftri3( ndim, nnodesub, xi, eta ) ;

                                            for ist = 1:nstr
                                                eps(con(iel,scon(ind)), ist) = eps(con(iel,scon(ind)), ist) + n * epsgp(:, ist) ;
                                                sig(con(iel,scon(ind)), ist) = sig(con(iel,scon(ind)), ist) + n * siggp(:, ist) ;
                                            end
                                        end
                                end
                            
                            % update nodal count for averaging
                                ndcount(con(iel,scon), 1) = ndcount(con(iel,scon), 1) + 1 ;
                        end
                end
        end
    
    % SOLUTION AVERAGING AT THE NODES
        nds = find(ndcount);
        for ind = nds'
            sig(ind,:) = sig(ind,:) / ndcount(ind,1) ;
            eps(ind,:) = eps(ind,:) / ndcount(ind,1) ;
        end

end


function [b] = bmatrixquad( nnodel, dnxy, nstr, ndofel )

    % **************************************************************************
    %
    %   INPUTS:
    %       nnodel - number of nodes per element
    %       dnxy - matrix of shape function derivatives
    %       nstr - number os stress/strain variables
    %       ndofel - number of dof per element
    %
    %   OUTPUTS:
    %       b - strain derivative operator matrix
    %
    % **************************************************************************
    
    % VARIABLE INITIALISATION
        b = zeros(nstr, ndofel) ;
    
    % LOOP THROUGH NODES
        for ind = 1:nnodel
            b(:, (2*ind - 1):(2*ind) ) = [ dnxy(1,ind),           0 ;
                                                     0, dnxy(2,ind) ;
                                           dnxy(2,ind), dnxy(1,ind)]; 
        end
end


function [b] = bmatrixtri( nnodel, dnxy, nstr, ndofel )

    % **************************************************************************
    %
    %   INPUTS:
    %       nnodel - number of nodes per element
    %       dnxy - matrix of shape function derivatives
    %       nstr - number os stress/strain variables
    %       ndofel - number of dof per element
    %
    %   OUTPUTS:
    %       b - strain derivative operator matrix
    %
    % **************************************************************************

    % VARIABLE INITIALISATION
        b = zeros(nstr, ndofel) ;
    
    % LOOP THROUGH NODES
        for ind = 1:nnodel
            b(:, (2*ind - 1):(2*ind) ) = [ dnxy(1,ind),           0 ;
                                                     0, dnxy(2,ind) ;
                                           dnxy(2,ind), dnxy(1,ind)]; 
        end
end

