% **************************************************************************
%  MODELLING EVOLVING BOUNDARIES USING FNM AND LSM (April 2020)
% **************************************************************************
%
%
%   fem_mesh.m
%
%
%   FUNCTIONS:
%       - fem_mesh - creates a structured and regular grid of quadrilateral
%                    elements and their connectivities for the FEM analysis
%       - nodalcoord - creates the matrix of nodal coordinates
%       - connectivities - creates the matrix of global connectivities
%                          accounting only for the real nodes
%       - intgmesh - creates the matrices of coordinate and weights of the
%                    integration points for both quadrilateral and triangular
%                    isoparametric elements under the Gauss quadrature scheme
%
%
% **************************************************************************
%                   developed by Rui O.S.S. da Costa and Silvestre T. Pinho
%                                      Imperial College London, Aeronautics
% **************************************************************************


function [x, con, dofel, dofnode, xgpquad, wgpquad, xgptri, wgptri] = fem_mesh(s_mesh, s_geo, s_param)

    % **************************************************************************
    %
    %   INPUTS:
    %       s_mesh - structure arry with mesh data
    %       s_geo - structure arry with geometric data
    %       s_param - structure arry with generic parameters
    %
    %   OUTPUTS:
    %       x - matrix of nodal coordinates
    %       con - matrix of global connectivities 
    %             (accounting only for real nodes)
    %       dofel - matrix of global element-dof connectivities
    %               (accounting only for real dofs)
    %       dofnode - matrix of global node-dof connectivities
    %                 (accounting only for real nodes)
    %       xgpquad - matrix of coordinates for the integration points
    %                 for a quadrilateral element
    %       wgpquad - weight of the integration points
    %                 for a quadrilateral element
    %       xgptri - matrix of coordinates for the integration points
    %                for triangular elements
    %       wgptri - weight of the integration points
    %                for triangular elements
    %
    % **************************************************************************
    
    % VARIABLE INITIALISATION
        nnode = s_mesh.nnode;
        ndim = s_mesh.ndim;
        divx = s_mesh.divx;
        divy = s_mesh.divy;
        Lx = s_geo.Lx;
        Ly = s_geo.Ly;
        nel = s_mesh.nel;
        nnodel = s_mesh.nnodel;
        ndofel = s_mesh.ndofel;
        ndofnode = s_mesh.ndofnode;
        ngptquad = s_mesh.ngptquad;
        ngpttri = s_mesh.ngpttri;
        nfloatnode = s_mesh.nfloatnode ;
        nfloatnodel = s_mesh.nfloatnodel ;
        nfloatdofel = s_mesh.nfloatdofel ;
        
    % CREATE MESH
        [x] = nodalcoord( nnode, ndim, divx, divy, Lx, Ly ) ;
    
    % CREATE CONNECTIVITIES
        [con, dofel, dofnode] = connectivities( nel, nnode - nfloatnode, nnodel - nfloatnodel, ndofel - nfloatdofel, ndofnode, divx);
    
    % INTEGRATION POINTS AND WEIGHTS
        [xgpquad, wgpquad] = intgmesh( ngptquad, ndim, 'quad4' ) ;
        [xgptri, wgptri] = intgmesh( ngpttri, ndim, 'tri3' ) ;

end


function [x] = nodalcoord( nnode, ndim, divx, divy, Lx, Ly )

    % **************************************************************************
    %
    %   INPUTS:
    %       nnode - total number of nodes
    %       ndim - number of spatial dimensions
    %       divx - number of elements along the horizontal direction
    %       divy - number of elements along the vertical direction
    %       Lx - length of the domain in the horizontal direction
    %       Ly - length of the domain in the vertical direction
    %
    %   OUTPUTS:
    %       x - matrix of nodal coordinates
    %
    % **************************************************************************

    % VARIABLE INITIALISATION
        x = zeros(nnode, ndim);
        nodx = divx + 1;
        nody = divy + 1;
        dx = Lx / divx ;
        dy = Ly / divy ;

    % CREATION OF THE STRUCTURED GRID
        for iy = 1:nody
            for ix = 1:nodx
                x( nodx*(iy - 1) + ix , 1) = (ix - 1)*dx ;
                x( nodx*(iy - 1) + ix , 2) = (iy - 1)*dy ;
            end
        end
end


function [con, dofel, dofnode] = connectivities( nel, nnode, nnodel, ndofel, ndofnode, divx)

    % **************************************************************************
    %
    %   INPUTS:
    %       nel - number of elements
    %       nnode - total number of nodes
    %       nnodel - number of nodes per element
    %       ndofel - number of dof per element
    %       ndofnode - number of dof per node
    %       divx - number of elements along the horizontal direction
    %
    %   OUTPUTS:
    %       con - matrix of global connectivities 
    %             (accounting only for real nodes)
    %       dofel - matrix of global element-dof connectivities
    %               (accounting only for real dofs)
    %       dofnode - matrix of global node-dof connectivities
    %                 (accounting only for real nodes)
    %
    % **************************************************************************

    % VARIABLE INITIALISATION
        con = zeros(nel, nnodel);
        dofel = zeros(nel, ndofel); 
        dofnode = zeros(nnode, ndofnode);
        nodx = divx + 1;   
        yrow = 0;
        xcol = 1;
        
    % LOOP THROUGH ELEMENTS (BUILD CON AND DOFEL)
        for iel = 1:nel
            
            % elemental connectivities
                con(iel, 1) = iel + yrow;
                con(iel, 2) = con(iel,1) + 1;
                con(iel, 3) = con(iel,2) + nodx;
                con(iel, 4) = con(iel,1) + nodx;

            % element-dof connectivities
                for ind = 1:nnodel
                    for idf = 1:ndofnode
                        dofel(iel, 2*ind + idf - 2 ) = 2*con(iel,ind) + idf - 2 ;
                    end
                end
            
            % incrementing auxiliar variables
                if xcol >= divx
                    yrow = yrow + 1;
                    xcol = 1;
                else
                    xcol = xcol + 1;
                end
        end
    
    % LOOP THROUGH NODES (BUILD DOFNODE)
        for ind = 1:nnode
            for idf = 1:ndofnode
                dofnode(ind, idf) = 2*ind + idf - 2 ;
            end
        end

end


function [xgp, wgp] = intgmesh( ngpt, ndim, eltype )

    % **************************************************************************
    %
    %   INPUTS:
    %       ngpt - number of integration points
    %       ndim - number of spatial dimensions
    %       eltype - type of element (quad4, tri3, ...)
    %
    %   OUTPUTS:
    %       xgp - matrix of coordinates for the integration points
    %       wgp - weight of the integration points
    %
    % **************************************************************************

    % VARIABLE INITIALISATION
        xgp = zeros(ngpt,ndim);
        wgp = zeros(ngpt,1);
    
    % QUADRILATERAL WITH 4 NODES
        if strcmp(eltype,'quad4')

            if ngpt == 1 % 1 integration point
                xgp(1,:) = [0, 0];
                wgp(1) = 2;
            elseif ngpt == 4 % 2x2 rule, 4 integration points
                xgp(1,:) = [-1/sqrt(3), -1/sqrt(3)];
                xgp(2,:) = [ 1/sqrt(3), -1/sqrt(3)];
                xgp(3,:) = [ 1/sqrt(3),  1/sqrt(3)];
                xgp(4,:) = [-1/sqrt(3),  1/sqrt(3)];
                wgp(1) = 1;
                wgp(2) = 1;
                wgp(3) = 1;
                wgp(4) = 1;
            else % 2x2 rule, 4 integration points
                xgp(1,:) = [-1/sqrt(3), -1/sqrt(3)];
                xgp(2,:) = [ 1/sqrt(3), -1/sqrt(3)];
                xgp(3,:) = [ 1/sqrt(3),  1/sqrt(3)];
                xgp(4,:) = [-1/sqrt(3),  1/sqrt(3)];
                wgp(1) = 1;
                wgp(2) = 1;
                wgp(3) = 1;
                wgp(4) = 1;
            end
    
    % TRIANGLE WITH 3 NODES
        elseif strcmp(eltype,'tri3')

            if ngpt == 1 % 1 integration points
                xgp(1,:) = [1/3, 1/3];
                wgp(1,:) = 1/2;
            elseif ngpt == 3 % 3 integration points
                xgp(1,:) = [1/6, 1/6];
                xgp(2,:) = [2/3, 1/6];
                xgp(3,:) = [1/6, 2/3];
                wgp(1) = 1/6;
                wgp(2) = 1/6;
                wgp(3) = 1/6;
            else % 3 integration points
                xgp(1,:) = [1/6, 1/6];
                xgp(2,:) = [2/3, 1/6];
                xgp(3,:) = [1/6, 2/3];
                wgp(1) = 1/6;
                wgp(2) = 1/6;
                wgp(3) = 1/6;
            end
        end
end