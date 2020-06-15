% **************************************************************************
%  MODELLING EVOLVING BOUNDARIES USING FNM AND LSM (April 2020)
% **************************************************************************
%
%
%   fnm_confnm.m
%
%
%   FUNCTIONS:
%       - fnm_confnm - extension of the connectivity matrices to account
%                      for the addition of floating nodes to the mesh
%
%
% **************************************************************************
%                   developed by Rui O.S.S. da Costa and Silvestre T. Pinho
%                                      Imperial College London, Aeronautics
% **************************************************************************


function [con, dofel, dofnode] = fnm_confnm( s_mesh, s_geo, s_param, s_sol )

    % **************************************************************************
    %
    %   INPUTS:
    %       s_mesh - structure arry with mesh data
    %       s_geo - structure arry with geometric data
    %       s_param - structure arry with generic parameters
    %       s_sol - structure array with solution data
    %
    %   OUTPUTS:
    %       con - matrix of global connectivities 
    %             (accounting for real and floating nodes)
    %       dofel - matrix of global element-dof connectivities
    %               (accounting for real and floating dofs)
    %       dofnode - matrix of global node-dof connectivities
    %                 (accounting for real and floating nodes)
    %
    % **************************************************************************    

    % VARIABLE INITIALISATION
        divx = s_mesh.divx ;
        con = s_mesh.con ;
        nel = s_mesh.nel ;
        nnode = s_mesh.nnode ;
        nfloatnode = s_mesh.nfloatnode ;
        nfloatnodel = s_mesh.nfloatnodel ;
        dofel = s_mesh.dofel ;
        nfloatdofel = s_mesh.nfloatdofel ;
        nfloatdofnode = s_mesh.nfloatdofnode ;
        dofnode = s_mesh.dofnode ;
        confnm = zeros(nel, nfloatnodel) ;
        nnormalnd = nnode - nfloatnode ; 
        irow = 1;
        dofelfnm = zeros(nel, nfloatdofel) ;
        dofnodefnm = zeros(nfloatnode, nfloatdofnode) ;
    
    % LOOP THROUGH ELEMENTS (EXTEND CON)
        for iel = 1:nel
            ind = nnormalnd + con(iel,1) + (irow - 1)*divx ;
            confnm(iel, 1:4) = [ind, ind + divx + 1, ind + 2*divx + 1, ind + divx] ;
            if mod(iel,divx) == 0 
                irow = irow + 1;
            end
        end
        mxfnd = max(max(confnm)) ;
        confnm(:,end) = ( (mxfnd + 1):(mxfnd + nel) )' ; % list of floating nodes assigned to element areas
    
    % ASSIGNING THE EXTENSION OF ELEMENTAL CONNECTITIVIES
        aux = con ;
        con = [aux, confnm] ;
     
    % LOOP THROUGH ELEMENTS (EXTEND DOFEL)
        for iel = 1:nel
            for ind = 1:nfloatnodel
                for idof = 1:nfloatdofnode
                    dofelfnm(iel, nfloatdofnode * ind + idof - nfloatdofnode) = nfloatdofnode * confnm(iel, ind) + idof - nfloatdofnode ;
                end 
            end 
        end
    
    % ASSIGNING THE EXTENSION OF ELEMENT-DOF CONNECTITIVIES
        aux = dofel ;
        dofel = [aux, dofelfnm ];
    
    % LOOP THROUGH FLOATING NODES (EXTEND DOFNODE)
        for ind = (nnormalnd + 1):nnode
            for idof = 1:nfloatdofnode
                dofnodefnm(ind - nnormalnd, idof) = nfloatdofnode * ind + idof - nfloatdofnode ;
            end
        end
    
    % ASSIGNING THE EXTENSION OF NODE-DOF CONNECTITIVIES
        aux = dofnode ;
        dofnode = [aux; dofnodefnm ];
end