% **************************************************************************
%  MODELLING EVOLVING BOUNDARIES USING FNM AND LSM (April 2020)
% **************************************************************************
%
%
%   fem_k.m
%
%
%   FUNCTIONS:
%       - fem_k - creates the global stiffness matrix for the
%                 mechanical static analysis
%
%
% **************************************************************************
%                   developed by Rui O.S.S. da Costa and Silvestre T. Pinho
%                                      Imperial College London, Aeronautics
% **************************************************************************


function [k] = fem_k( s_mesh, s_geo, s_param, s_sol )

    % **************************************************************************
    %
    %   INPUTS:
    %       s_mesh - structure arry with mesh data
    %       s_geo - structure arry with geometric data
    %       s_param - structure arry with generic parameters
    %       s_sol - structure arry with solution data
    %
    %   OUTPUTS:
    %       k - global stiffness matrix
    %
    % **************************************************************************
    
    % VARIABLE INITIALISATION
        con = s_mesh.con;
        ndof = s_mesh.ndof ;
        ndofel = s_mesh.ndofel;
        thck = s_geo.thck;
        ngptquad = s_mesh.ngptquad;
        ngpttri = s_mesh.ngpttri;
        xgpquad = s_mesh.xgpquad;
        wgpquad = s_mesh.wgpquad;
        xgptri = s_mesh.xgptri;
        wgptri = s_mesh.wgptri;
        ndim = s_mesh.ndim;
        nnodel = s_mesh.nnodel;
        x = s_mesh.x;
        nstr = s_mesh.nstr ;
        D = s_geo.D;
        nfloatdofel = s_mesh.nfloatdofel; 
        nfloatnodel = s_mesh.nfloatnodel ;
        elstat = s_mesh.elstat ;
        ndofnode = s_mesh.ndofnode ;
        phi = s_param.p ;
        elsin = s_mesh.elsin;
        elsbnd = s_mesh.elsbnd;   
        rows = zeros(ceil(0.001*ndof*ndof), 1) ;
        cols = rows ;
        vals = rows ;
        count = 1;
    
    
    % ELEMENTAL STIFFNESS MATRIX
        % K for the non-paritioned quadrilateral element
            [kel] = fem_kquad4( ndofel - nfloatdofel, thck, ngptquad, xgpquad, wgpquad, ndim, nnodel - nfloatnodel, x(con(elsin(1),1:4), :), D, nstr ) ;
            [i,j,v] = find(kel);
            
        % loop through non-paritioned elements
            for iel = elsin'
                
                % local indices
                    scon = [1 2 3 4] ;
                    nnodesub = size(scon,2) ;
                    ndofsub = nnodesub * ndofnode ;

                % global indices
                    glbnds = con(iel, scon) ;
                    glbdfs = zeros(ndofsub,1) ;
                    for ind = 1:nnodesub
                        for idof = 1:ndofnode
                            glbdfs(ndofnode * ind + idof - ndofnode,1) = ndofnode * glbnds(ind) + idof - ndofnode ;
                        end
                    end
                    
                % assigning indices and values for assembly
                    rows(count:(count - 1 + size(i,1)), 1) = glbdfs(i,1);
                    cols(count:(count - 1 + size(j,1)), 1) = glbdfs(j,1);
                    vals(count:(count - 1 + size(v,1)), 1) = v ;
                    count = count + size(i,1);
            end
            
        % loop through partitioned elements   
            for iel = elsbnd'
                
                % loop through sub elements
                    nsub = fnm_partitioninfo(elstat(iel,1)) ;
                    for isub = 1:nsub
                        
                        % local indices
                            [scon] = fnm_subcon( elstat(iel,1), isub ) ;
                            nnodesub = size(scon,2) ;
                            ndofsub = nnodesub * ndofnode ;
                            
                        % coordinates, LS values and pseudo-density
                            xel = x(con(iel,scon), :) ;
                            phiel = phi(con(iel,scon), :) ;
                            rho = ls_density( phiel ) ;

                        % computation of K for solid elements only
                            if rho == 1
                                if nnodesub == 4 % quad sub element
                                    [kel] = fem_kquad4( ndofsub, thck, ngptquad, xgpquad, wgpquad, ndim, nnodesub, xel, D, nstr ) ;
                                elseif nnodesub == 3 % tri sub element
                                    [kel] = fem_ktri3( ndofsub, thck, ngpttri, xgptri, wgptri, ndim, nnodesub, xel, D, nstr ) ;
                                end
                                [i,j,v] = find(kel);

                                % global indices
                                    glbnds = con(iel, scon) ;
                                    glbdfs = zeros(ndofsub,1) ;
                                    for ind = 1:nnodesub
                                        for idof = 1:ndofnode
                                            glbdfs(ndofnode * ind + idof - ndofnode,1) = ndofnode * glbnds(ind) + idof - ndofnode ;
                                        end
                                    end
                                    
                                % assigning indices and values for assembly
                                    rows(count:(count - 1 + size(i,1)), 1) = glbdfs(i,1);
                                    cols(count:(count - 1 + size(j,1)), 1) = glbdfs(j,1);
                                    vals(count:(count - 1 + size(v,1)), 1) = v ;
                                    count = count + size(i,1);
                            end
                    end
            end
            
            
    % ASSEMBLY
        rows(count:end,:) = [];
        cols(count:end,:) = [];
        vals(count:end,:) = [];
        k = sparse(rows,cols,vals,ndof,ndof);
end