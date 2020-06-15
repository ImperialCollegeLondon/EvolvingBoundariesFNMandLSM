% **************************************************************************
%  MODELLING EVOLVING BOUNDARIES USING FNM AND LSM (April 2020)
% **************************************************************************
%
%
%   fem_f.m
%
%
%   FUNCTIONS:
%       - fem_f - creates the global force vector for the
%                 mechanical static analysis
%       - pointload - creates a force vector for the case of a point load
%
%
% **************************************************************************
%                   developed by Rui O.S.S. da Costa and Silvestre T. Pinho
%                                      Imperial College London, Aeronautics
% **************************************************************************

function [fe] = fem_f(s_mesh, s_geo, s_param, s_sol)
    
    % **************************************************************************
    %
    %   INPUTS:
    %       s_mesh - structure arry with mesh data
    %       s_geo - structure arry with geometric data
    %       s_param - structure arry with generic parameters
    %       s_sol - structure arry with solution data
    %
    %   OUTPUTS:
    %       fe - global force vector
    %
    % **************************************************************************

    % VARIABLE INITIALISATION
        ndof = s_mesh.ndof;
        x = s_mesh.x;
        Lx = s_geo.Lx;
        Ly = s_geo.Ly;
        dofnode = s_mesh.dofnode;
        F = s_param.F;
        ubc = s_param.ubc; 
        bndinfo = s_geo.bndinfo;
        anode = (1:(s_mesh.nnode - s_mesh.nfloatnode))';
        h = s_param.h;
    
    % POINT LOAD
        fe = pointload( ndof, x, Lx, Ly, h, dofnode, F, ubc, anode, bndinfo );
        fe = sparse(fe);
end

function [fp] = pointload( ndof, x, Lx, Ly, h, dofnode, F, ubc, anode, bndinfo )
    
    % **************************************************************************
    %
    %   INPUTS:
    %       ndof - total number of dof (real and floating)
    %       x - matrix of nodal coordinates
    %       Lx - length of the domain in the horizontal direction
    %       Ly - length of the domain in the vertical direction
    %       h - minimum element size
    %       dofnode - matrix of global node-dof connectivities
    %       F - force value
    %       ubc - boundary condition selector
    %       anode - vector of active nodes
    %       bndinfo - vector that lists some parameters regarding the fixed
    %                 boundaries due to the geometric restrictions
    %
    %   OUTPUTS:
    %       fp - force vector with point load
    %
    % **************************************************************************

    % VARIABLE INITIALISATION
        fp = zeros(ndof,1);

    % POINT LOAD VECTOR FOR EACH TYPE OF PROBLEM
        if ubc == 1 % cantilever beam

            nodeindx = find( x(:,1) == Lx );
            nodeindx2 = find( abs(x(nodeindx,2) - Ly/2) == min(abs(x(nodeindx,2) - Ly/2)) );
            nodeindx = nodeindx(nodeindx2);
            nodeindx = nodeindx(1);
            dofindx = dofnode(nodeindx,2);
            fp(dofindx,1) = F;

        elseif ubc == 2 % mbb beam

            nodeindx = find( x(:,1) == 0 & x(:,2) == Ly );
            nodeindx = nodeindx(1);
            dofindx = dofnode(nodeindx,2);
            fp(dofindx,1) = F;

        elseif ubc == 3 % L bracket

            d = abs(x(:,2) - Ly/5);
            nodeindx = find( d == min(d) & x(:,1) == Lx );
            nodeindx = nodeindx(1);
            dofindx = dofnode(nodeindx,2);
            fp(dofindx,1) = F;

        elseif ubc == 4 % bracket with 2 holes

            d = abs(x(anode,2) - (Ly/2 + bndinfo(4,1) + h/20));
            nodeindx1 = anode(find( d == min(d) & abs(x(anode,1) - bndinfo(2,1)) < h ),1) ;
            d = abs(x(anode,2) - (Ly/2 - bndinfo(4,1) - h/20));
            nodeindx2 = anode(find( d == min(d) & abs(x(anode,1) - bndinfo(2,1)) < h ),1) ;
            nodeindx = [nodeindx1; nodeindx2];
            dofindx = dofnode(nodeindx, 2) ;
            fp(dofindx,1) = F;

        end
        fp = sparse(fp);
    
end
