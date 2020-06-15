% **************************************************************************
%  MODELLING EVOLVING BOUNDARIES USING FNM AND LSM (April 2020)
% **************************************************************************
%
%
%   ls_f.m
%
%
%   FUNCTIONS:
%       - ls_f - creates the global force vector for the
%                level set evolution
%
%
% **************************************************************************
%                   developed by Rui O.S.S. da Costa and Silvestre T. Pinho
%                                      Imperial College London, Aeronautics
% **************************************************************************

function [fe] = ls_f(s_mesh, s_geo, s_param, s_sol)

    % **************************************************************************
    %
    %   INPUTS:
    %       s_mesh - structure arry with mesh data
    %       s_geo - structure arry with geometric data
    %       s_param - structure arry with generic parameters
    %       s_sol - structure arry with solution data
    %
    %   OUTPUTS:
    %       fe - global LS force vector
    %
    % **************************************************************************

    % VARIABLE INITIALISATION
        ndof = s_mesh.ndof;
        nel = s_mesh.nel;
        con = s_mesh.con;
        u = s_sol.u;
        x = s_mesh.x;
        ndofel = s_mesh.ndofel;
        ngptquad = s_mesh.ngptquad;
        xgpquad = s_mesh.xgpquad;
        wgpquad = s_mesh.wgpquad;
        nnodel = s_mesh.nnodel;
        h = s_param.h;
        vn = s_param.vn;
        dt = s_param.dt;
        dtau = s_param.dtau;
        rinit = s_param.rinit;
        ndim = s_mesh.ndim;
        dofel = s_mesh.dofel;
        ndofnode = s_mesh.ndofnode;
        nfloatdofel = s_mesh.nfloatdofel; 
        nfloatnodel = s_mesh.nfloatnodel ;
        elsoutdom = s_mesh.elsoutdom;
        els = setdiff((1:nel)',elsoutdom);
        rows = zeros(ceil(0.1*ndof), 1) ;
        vals = rows ;
        count = 1;
        
    % ELEMENT LOOP
        for iel = els'
            
            % velocity, LS values and elemental coordinates
                vel = vn(dofel(iel,1:(ndofel - nfloatdofel)), 1);
                phiel = u(dofel(iel,1:(ndofel - nfloatdofel)), 1);
                xel = x(con(iel,1:(nnodel - nfloatnodel)), :);

            % elemental force vector
                fel = ls_fquad4( ndim, ndofel - nfloatdofel, ngptquad, xgpquad, wgpquad, nnodel - nfloatnodel, xel, phiel, h, vel, dt, dtau, rinit );

            % global indices
                glbnds = con(iel, 1:(nnodel - nfloatnodel)) ;
                glbdfs = zeros(ndofel - nfloatdofel,1) ;
                for ind = 1:(nnodel - nfloatnodel)
                    for idof = 1:ndofnode
                        glbdfs(ndofnode * ind + idof - ndofnode,1) = ndofnode * glbnds(ind) + idof - ndofnode ;
                    end
                end
                
            % assigning indices and values for assembly
                [i,j,v] = find(fel);
                rows(count:(count - 1 + size(i,1)), 1) = glbdfs(i,1);
                vals(count:(count - 1 + size(v,1)), 1) = v ;
                count = count + size(i,1);
        end
        
    % ASSEMBLY
        rows(count:end,:) = [];
        cols = 0*rows + 1 ;
        vals(count:end,:) = [];
        fe = sparse(rows,cols,vals);
    
    
end
