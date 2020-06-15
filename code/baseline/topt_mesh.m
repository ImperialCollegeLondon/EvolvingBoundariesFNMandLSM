%% TOPOPT
%
%   CONSTRUCTION OF THE MESH

function [x, con, dofel, l_dofel, dofnode, l_dofnode, xgp, wgp, adof] = topt_mesh(s_mesh, s_geo, s_param)
    
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
    ngpt = s_mesh.ngpt;
    ubc = s_param.ubc;
    dofs = s_mesh.dofs;
    l_ndofel = s_mesh.l_ndofel;
    l_ndofnode = s_mesh.l_ndofnode;


    [x] = nodalcoord( nnode, ndim, divx, divy, Lx, Ly ) ;
    [con, dofel, dofnode, l_dofel, l_dofnode] = connectivities( nel, nnode, nnodel, ndofel, l_ndofel, ndofnode, l_ndofnode, divx ) ;
    [xgp, wgp] = intgmesh( ngpt, ndim ) ;
    [adof] = activequations( ubc, x, dofnode, ndofnode, dofs, Lx, Ly ) ;

end

%%
function [x] = nodalcoord( nnode, ndim, divx, divy, Lx, Ly )
% creates the vector of nodal coordinates
% divides lengths Lx and Ly into the number
% of divisions divx and divy specified

    x = zeros(nnode, ndim);
    nodx = divx + 1;
    nody = divy + 1;
    dx = Lx / divx ;
    dy = Ly / divy ;

    for iy = 1:nody
        for ix = 1:nodx
            x( nodx*(iy - 1) + ix , 1) = (ix - 1)*dx ;
            x( nodx*(iy - 1) + ix , 2) = (iy - 1)*dy ;
        end
    end

end

%%
function [con, dofel, dofnode, l_dofel, l_dofnode] = connectivities( nel, nnode, nnodel, ndofel, l_ndofel, ndofnode, l_ndofnode, divx )
% establishes nodal and elemental connectivities
% creates matrix con with nodal indices for each element
% matrix dofel with dof indices for each element
% and matrix dofnode with dof indices for each node

    con = zeros(nel, nnodel);
    dofel = zeros(nel, ndofel); 
    l_dofel = zeros(nel, l_ndofel); 
    dofnode = zeros(nnode, ndofnode);
    l_dofnode = zeros(nnode, l_ndofnode);
    nodx = divx + 1;
    
    yrow = 0;
    xcol = 1;
    for iel = 1:nel
        
        con(iel, 1) = iel + yrow;
        con(iel, 2) = con(iel,1) + 1;
        con(iel, 3) = con(iel,2) + nodx;
        con(iel, 4) = con(iel,1) + nodx;
        
        for ind = 1:nnodel
            for idf = 1:ndofnode
                dofel(iel, 2*ind + idf - 2 ) = 2*con(iel,ind) + idf - 2 ;
            end
        end
        
        l_dofel(iel,:) = con(iel,:);
        
        if xcol >= divx
            yrow = yrow + 1;
            xcol = 1;
        else
            xcol = xcol + 1;
        end
    end
    
    for ind = 1:nnode
        for idf = 1:ndofnode
            dofnode(ind, idf) = 2*ind + idf - 2 ;
        end
        for idf = 1:l_ndofnode
            l_dofnode(ind, idf) = ind ;
        end
    end
    
    

end

%%
function [xgp, wgp] = intgmesh( ngpt, ndim )
% creates the mesh of integration points
% According to the Gauss quadrature constructs vector xgp and wgp
% containing the coordinates and weights of gauss points

    xgp = zeros(ngpt,ndim);
    wgp = zeros(ngpt,ndim);
    if ngpt == 1
        xgp(1,:) = [0, 0];
        wgp(1,:) = [2, 2];
    elseif ngpt == 4
        xgp(1,:) = [-1/sqrt(3), -1/sqrt(3)];
        xgp(2,:) = [ 1/sqrt(3), -1/sqrt(3)];
        xgp(3,:) = [ 1/sqrt(3),  1/sqrt(3)];
        xgp(4,:) = [-1/sqrt(3),  1/sqrt(3)];
        wgp(1,:) = [1, 1];
        wgp(2,:) = [1, 1];
        wgp(3,:) = [1, 1];
        wgp(4,:) = [1, 1];
    elseif ngpt == 9
        xgp(1,:) = [-sqrt(3/5), -sqrt(3/5)];
        xgp(2,:) = [         0, -sqrt(3/5)];
        xgp(3,:) = [ sqrt(3/5), -sqrt(3/5)];
        xgp(4,:) = [-sqrt(3/5),          0];
        xgp(5,:) = [         0,          0];
        xgp(6,:) = [ sqrt(3/5),          0];
        xgp(7,:) = [-sqrt(3/5),  sqrt(3/5)];
        xgp(8,:) = [         0,  sqrt(3/5)];
        xgp(9,:) = [ sqrt(3/5),  sqrt(3/5)];
        wgp(1,:) = [5/9, 5/9];
        wgp(2,:) = [8/9, 5/9];
        wgp(3,:) = [5/9, 5/9];
        wgp(4,:) = [5/9, 8/9];
        wgp(5,:) = [8/9, 8/9];
        wgp(6,:) = [5/9, 8/9];
        wgp(7,:) = [5/9, 5/9];
        wgp(8,:) = [8/9, 5/9];
        wgp(9,:) = [5/9, 5/9];
    else
        xgp(1,:) = [-1/sqrt(3), -1/sqrt(3)];
        xgp(2,:) = [ 1/sqrt(3), -1/sqrt(3)];
        xgp(3,:) = [ 1/sqrt(3),  1/sqrt(3)];
        xgp(4,:) = [-1/sqrt(3),  1/sqrt(3)];
        wgp(1,:) = [1, 1];
        wgp(2,:) = [1, 1];
        wgp(3,:) = [1, 1];
        wgp(4,:) = [1, 1];
    end

end

%%
function [adof] = activequations( ubc, x, dofnode, ndofnode, dofs, Lx, Ly )
% determines the active equations(dof)
% based on the matrix of boundary conditions creates a vector activdof
% which contains the number of the active dof (non-constrained dof)

    if ubc == 1
        
        nodeindx = find( x(:,1) == 0 );
        dofindaux = dofnode(nodeindx, :) ;
        dofindaux = reshape(dofindaux, [size(dofindaux,1)*size(dofindaux,2), 1]);
        dofindx = unique(dofindaux);
        
        adof = dofs;
        adof(dofindx) = [];
        
    elseif ubc == 2
        
        nodeindx = find( x(:,1) == Lx );
        dofindaux = dofnode(nodeindx, 1) ;
        dofindaux = [dofindaux; 1; 2];
        dofindx = unique(dofindaux);
        
        adof = dofs;
        adof(dofindx) = [];
        
    elseif ubc == 3
        
        nodeindx = find( x(:,1) == Lx & x(:,2) == 0 );
        dofindaux = dofnode(nodeindx, 2) ;
        dofindaux = [dofindaux; 1; 2];
        dofindx = unique(dofindaux);
        
        adof = dofs;
        adof(dofindx) = [];
        
    elseif ubc == 4
        
        nodeindx = find( x(:,1) == 0 );
        dofindaux = dofnode(nodeindx, :) ;
        dofindaux = reshape(dofindaux, [size(dofindaux,1)*size(dofindaux,2), 1]);
        dofindx = unique(dofindaux);
        
        adof = dofs;
        adof(dofindx) = [];
        
        
    else
        
        nodeindx = find( x(:,1) == 0 );
        dofindaux = dofnode(nodeindx, :) ;
        dofindaux = reshape(dofindaux, [size(dofindaux,1)*size(dofindaux,2), 1]);
        dofindx = unique(dofindaux);
        
        adof = dofs;
        adof(dofindx) = [];
        
    end

end