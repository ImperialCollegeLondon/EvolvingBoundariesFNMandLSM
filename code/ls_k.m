% **************************************************************************
%  MODELLING EVOLVING BOUNDARIES USING FNM AND LSM (April 2020)
% **************************************************************************
%
%
%   ls_k.m
%
%
%   FUNCTIONS:
%       - ls_k - creates the global stiffness matrix for the
%                level set evolution
%
%
% **************************************************************************
%                   developed by Rui O.S.S. da Costa and Silvestre T. Pinho
%                                      Imperial College London, Aeronautics
% **************************************************************************

function [k] = ls_k( s_mesh, s_geo, s_param, s_sol )
    
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
        nel = s_mesh.nel;
        con = s_mesh.con;
        ndof = s_mesh.ndof ;
        ndofel = s_mesh.ndofel;
        ngptquad = s_mesh.ngptquad;
        xgpquad = s_mesh.xgpquad;
        wgpquad = s_mesh.wgpquad;
        ndim = s_mesh.ndim;
        nnodel = s_mesh.nnodel;
        x = s_mesh.x;
        rinit = s_param.rinit;
        c = s_param.c;
        h = s_param.h;
        ndofnode = s_mesh.ndofnode;
        nfloatdofel = s_mesh.nfloatdofel; 
        nfloatnodel = s_mesh.nfloatnodel ;
        elsoutdom = s_mesh.elsoutdom;
        els = setdiff((1:nel)',elsoutdom);   
        rows = zeros(ceil(0.001*ndof*ndof), 1) ;
        cols = rows ;
        vals = rows ;
        count = 1;
    
    % ELEMENTAL STIFFNESS MATRIX
        % K for the non-paritioned quadrilateral element
            [kel] = ls_kquad4( ndofel - nfloatdofel, ngptquad, xgpquad, wgpquad, ndim, nnodel - nfloatnodel, x(con(1,1:4),:), rinit, c, h ) ;
    
        % loop all elements
            for iel = els'

                % global indices
                    glbnds = con(iel, 1:(nnodel - nfloatnodel)) ;
                    glbdfs = zeros(ndofel - nfloatdofel,1) ;
                    for ind = 1:(nnodel - nfloatnodel)
                        for idof = 1:ndofnode
                            glbdfs(ndofnode * ind + idof - ndofnode,1) = ndofnode * glbnds(ind) + idof - ndofnode ;
                        end
                    end
                    
                % assigning indices and values for assembly
                    [i,j,v] = find(kel);
                    rows(count:(count - 1 + size(i,1)), 1) = glbdfs(i,1);
                    cols(count:(count - 1 + size(j,1)), 1) = glbdfs(j,1);
                    vals(count:(count - 1 + size(v,1)), 1) = v ;
                    count = count + size(i,1);
            end

            
    % ASSEMBLY
        rows(count:end,:) = [];
        cols(count:end,:) = [];
        vals(count:end,:) = [];
        k = sparse(rows,cols,vals);   
end