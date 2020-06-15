% **************************************************************************
%  MODELLING EVOLVING BOUNDARIES USING FNM AND LSM (April 2020)
% **************************************************************************
%
%
%   ls_kquad4.m
%
%
%   FUNCTIONS:
%       - ls_kquad4 - creates the elemental stiffness matrix of a
%                     quadrilateral isoparametric element for the
%                     level set evolution
%
%
% **************************************************************************
%                   developed by Rui O.S.S. da Costa and Silvestre T. Pinho
%                                      Imperial College London, Aeronautics
% **************************************************************************


function [kel] = ls_kquad4( ndofel, ngpt, xgp, wgp, ndim, nnodel, xel, rinit, c, h )

    % **************************************************************************
    %
    %   INPUTS:
    %       ndofel - number of total dof per element
    %       ngpt -  number of integration points
    %       xgp - matrix of coordinates for the integration points
    %       wgp - matrix of weights for the integration points
    %       ndim - number of spatial dimensions
    %       nnodel - number of total nodes per element
    %       xel - matrix of nodal coordinates
    %       rinit - boolean value to toggle reinitialisation
    %       c - LSFEM parameter
    %       h - mesh size
    %
    %
    %   OUTPUTS:
    %       kel - elemental stiffness matrix
    %
    % **************************************************************************
    
    % VARIABLE INITIALISATION
        kel = zeros(ndofel, ndofel);
        
    % LOOP THROUGH INTEGRATION POINTS
        for igp = 1:ngpt
            % coordinates of integration points
                xi = xgp(igp,1);
                eta = xgp(igp,2);
                
            % shape functions
                [n, dn] = fem_sfquad4( ndim, nnodel, xi, eta ) ;
                
            % jacobian
                J = dn * xel ; detJ = det(J) ;
                dnxy = J\dn ;

            % stiffness matrix computation
                ka = transpose(n) * n * detJ * wgp(igp);
                
                % reinitialisation
                    if rinit
                        kb = transpose(dnxy) * dnxy * detJ * wgp(igp);
                        kphi = ka + c * h^2 * kb;
                    
                % regular
                    else
                        kphi = eye(nnodel) .* sum(ka,2); % lumped mass method
                    end

            % atribution
                kel = kel + kphi;
        end
end