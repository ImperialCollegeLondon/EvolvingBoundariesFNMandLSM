% **************************************************************************
%  MODELLING EVOLVING BOUNDARIES USING FNM AND LSM (April 2020)
% **************************************************************************
%
%
%   ls_mesh.m
%
%
%   FUNCTIONS:
%       - ls_mesh - creates a structured and regular grid of quadrilateral
%                    elements and their connectivities for the LSFEM analysis
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

function [x, con, dofel, dofnode, xgpquad, wgpquad, xgptri, wgptri, adof] = fem_mesh(s_mesh, s_geo, s_param)
    
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
    %       adof - vector of active DoFs
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
        ndofnode = s_mesh.ndofnode;
        ngptquad = s_mesh.ngptquad;
        ngpttri = s_mesh.ngpttri;
        nfloatnode = s_mesh.nfloatnode ;
        nfloatnodel = s_mesh.nfloatnodel ;
    

    % CREATE MESH
        [x] = nodalcoord( nnode, ndim, divx, divy, Lx, Ly ) ;
    
    % CREATE CONNECTIVITIES
        [con, dofel, dofnode] = connectivities( nel, nnode - nfloatnode, nnodel - nfloatnodel, ndofnode, divx);
    
    % INTEGRATION POINTS AND WEIGHTS
        [xgpquad, wgpquad] = intgmesh( ngptquad, ndim, 'quad4' ) ;
        [xgptri, wgptri] = intgmesh( ngpttri, ndim, 'tri3' ) ;
    
    % ACTIVE DOF
        adof = s_mesh.adof;

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


function [con, dofel, dofnode] = connectivities( nel, nnode, nnodel, ndofnode, divx)

    % **************************************************************************
    %
    %   INPUTS:
    %       nel - number of elements
    %       nnode - total number of nodes
    %       nnodel - number of nodes per element
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
            
            % incrementing auxiliar variables
                if xcol >= divx
                    yrow = yrow + 1;
                    xcol = 1;
                else
                    xcol = xcol + 1;
                end
        end
        
    % ELEMENT-DOF CONNECTIVITIES
        dofel = con;
        
    % LOOP THROUGH NODES (BUILD DOFNODE)
        for ind = 1:nnode
            for idf = 1:ndofnode
                dofnode(ind, idf) = ind ;
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
        if strcmp(eltype,'quad4') || strcmp(eltype,'lsquad4')

            if ngpt == 1 % 1 integration point
                xgp(1,:) = [0, 0];
                wgp(1) = 4;
            elseif ngpt == 4 % 2x2 rule, 4 integration points
                xgp(1,:) = [-1/sqrt(3), -1/sqrt(3)];
                xgp(2,:) = [ 1/sqrt(3), -1/sqrt(3)];
                xgp(3,:) = [ 1/sqrt(3),  1/sqrt(3)];
                xgp(4,:) = [-1/sqrt(3),  1/sqrt(3)];
                wgp(1) = 1;
                wgp(2) = 1;
                wgp(3) = 1;
                wgp(4) = 1;
            elseif ngpt == 9 % 3x3 rule, 9 integration points
                xgp = [-sqrt(3/5), -sqrt(3/5);
                                0, -sqrt(3/5);
                        sqrt(3/5), -sqrt(3/5);
                       -sqrt(3/5),          0;
                                0,          0;
                        sqrt(3/5),          0;
                       -sqrt(3/5),  sqrt(3/5);
                                0,  sqrt(3/5);
                        sqrt(3/5),  sqrt(3/5)] ;
                wgp = [25/81;
                       40/81;
                       25/81;
                       40/81;
                       64/81;
                       40/81;
                       25/81;
                       40/81;
                       25/81] ;
            elseif ngpt == 16 % 4x4 rule, 4 integration points
                xgp = [-0.8611363115940526, -0.8611363115940526;
                       -0.3399810435848563, -0.8611363115940526;
                        0.3399810435848563, -0.8611363115940526;
                        0.8611363115940526, -0.8611363115940526;
                       -0.8611363115940526, -0.3399810435848563;
                       -0.3399810435848563, -0.3399810435848563;
                        0.3399810435848563, -0.3399810435848563;
                        0.8611363115940526, -0.3399810435848563;
                       -0.8611363115940526,  0.3399810435848563;
                       -0.3399810435848563,  0.3399810435848563;
                        0.3399810435848563,  0.3399810435848563;
                        0.8611363115940526,  0.3399810435848563;
                       -0.8611363115940526,  0.8611363115940526;
                       -0.3399810435848563,  0.8611363115940526;
                        0.3399810435848563,  0.8611363115940526;
                        0.8611363115940526,  0.8611363115940526] ;
                 wgp = [0.3478548451374538 * 0.3478548451374538;
                        0.6521451548625461 * 0.3478548451374538;
                        0.6521451548625461 * 0.3478548451374538;
                        0.3478548451374538 * 0.3478548451374538;
                        0.3478548451374538 * 0.6521451548625461;
                        0.6521451548625461 * 0.6521451548625461;
                        0.6521451548625461 * 0.6521451548625461;
                        0.3478548451374538 * 0.6521451548625461;
                        0.3478548451374538 * 0.6521451548625461;
                        0.6521451548625461 * 0.6521451548625461;
                        0.6521451548625461 * 0.6521451548625461;
                        0.3478548451374538 * 0.6521451548625461;
                        0.3478548451374538 * 0.3478548451374538;
                        0.6521451548625461 * 0.3478548451374538;
                        0.6521451548625461 * 0.3478548451374538;
                        0.3478548451374538 * 0.3478548451374538] ;
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
        elseif strcmp(eltype,'tri3') || strcmp(eltype,'lstri3')

            if ngpt == 1 % 1 integration points
                xgp(1,:) = [1/3, 1/3];
                wgp(1,:) = 1;
            elseif ngpt == 3 % 3 integration points
                xgp(1,:) = [1/6, 1/6];
                xgp(2,:) = [2/3, 1/6];
                xgp(3,:) = [1/6, 2/3];
                wgp(1) = 1/3;
                wgp(2) = 1/3;
                wgp(3) = 1/3;
            elseif ngpt == 6 % 6 integration points
                xgp = [0.44594849091597, 0.44594849091597;
                       0.44594849091597, 0.10810301816807;
                       0.10810301816807, 0.44594849091597;
                       0.09157621350977, 0.09157621350977;
                       0.09157621350977, 0.81684757298046;
                       0.81684757298046, 0.09157621350977];
                wgp = [0.22338158967801;
                       0.22338158967801;
                       0.22338158967801;
                       0.10995174365532;
                       0.10995174365532;
                       0.10995174365532];
            elseif ngpt == 7 % 7 integration points
                xgp = [0.33333333333333, 0.33333333333333;
                       0.47014206410511, 0.47014206410511;
                       0.47014206410511, 0.05971587178977;
                       0.05971587178977, 0.47014206410511;
                       0.10128650732346, 0.10128650732346;
                       0.10128650732346, 0.79742698535309;
                       0.79742698535309, 0.10128650732346];
                wgp = [0.22500000000000;
                       0.13239415278851;
                       0.13239415278851;
                       0.13239415278851;
                       0.12593918054483;
                       0.12593918054483;
                       0.12593918054483];
            elseif ngpt == 12 % 12 integration points
                xgp = [0.24928674517091, 0.24928674517091;    
                       0.24928674517091, 0.50142650965818;
                       0.50142650965818, 0.24928674517091;
                       0.06308901449150, 0.06308901449150;
                       0.06308901449150, 0.87382197101700;
                       0.87382197101700, 0.06308901449150;
                       0.31035245103378, 0.63650249912140;
                       0.63650249912140, 0.05314504984482;
                       0.05314504984482, 0.31035245103378;
                       0.63650249912140, 0.31035245103378;
                       0.31035245103378, 0.05314504984482;
                       0.05314504984482, 0.63650249912140];
                wgp = [0.11678627572638;
                       0.11678627572638;
                       0.11678627572638;
                       0.05084490637021;
                       0.05084490637021;
                       0.05084490637021;
                       0.08285107561837;
                       0.08285107561837;
                       0.08285107561837;
                       0.08285107561837;
                       0.08285107561837;
                       0.08285107561837];

            else % 3 integration points
                xgp(1,:) = [1/6, 2/3];
                xgp(2,:) = [1/6, 1/6];
                xgp(3,:) = [2/3, 1/6];
                wgp(1) = 1/3;
                wgp(2) = 1/3;
                wgp(3) = 1/3;
            end

        end

end