% **************************************************************************
%  MODELLING EVOLVING BOUNDARIES USING FNM AND LSM (April 2020)
% **************************************************************************
%
%
%   fem_kquad4.m
%
%
%   FUNCTIONS:
%       - fem_ktri3 - creates the elemental stiffness matrix of a
%                     triangular isoparametric element for the
%                     mechanical static analysis
%       - bmatrixtri - computes the B matrix for the triangular
%                      element
%
%
% **************************************************************************
%                   developed by Rui O.S.S. da Costa and Silvestre T. Pinho
%                                      Imperial College London, Aeronautics
% **************************************************************************


function [kel] = fem_ktri3( ndofel, thck, ngpt, xgp, wgp, ndim, nnodel, xel, D, nstr )

    % **************************************************************************
    %
    %   INPUTS:
    %       ndofel - number of total dof per element
    %       thck - thickness
    %       ngpt -  number of integration points
    %       xgp - matrix of coordinates for the integration points
    %       wgp - matrix of weights for the integration points
    %       ndim - number of spatial dimensions
    %       nnodel - number of total nodes per element
    %       xel - matrix of nodal coordinates
    %       D - elastic constitutive matrix
    %       nstr - number os stress/strain variables
    %
    %
    %   OUTPUTS:
    %       kel - elemental stiffness matrix
    %
    % **************************************************************************

    % VARIABLE INITIALISATION
        kel = sparse(zeros(ndofel,ndofel));
    
    % LOOP THROUGH INTEGRATION POINTS
        for igp = 1:ngpt
            % coordinates of integration points
                xi = xgp(igp,1);
                eta = xgp(igp,2);
            
            % shape functions
                [n, dn] = fem_sftri3( ndim, nnodel, xi, eta ) ;
            
            % jacobian
                J = dn * xel ; detJ = det(J) ;
            
            % b matrix
                dnxy = J\dn ;
                [b] = bmatrixtri( nnodel, dnxy, nstr, ndofel ) ;

            % mechanical part of the stiffness matrix
                kmech = thck * transpose(b) * D * b * detJ * wgp(igp,1) ;

            % atribution
                kel = kel + kmech;
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