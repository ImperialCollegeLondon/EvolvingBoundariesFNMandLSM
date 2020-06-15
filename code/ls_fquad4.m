% **************************************************************************
%  MODELLING EVOLVING BOUNDARIES USING FNM AND LSM (April 2020)
% **************************************************************************
%
%
%   ls_fquad4.m
%
%
%   FUNCTIONS:
%       - ls_fquad4 - creates the elemental force vector of a
%                     quadrilateral isoparametric element for the
%                     level set evolution
%
%
% **************************************************************************
%                   developed by Rui O.S.S. da Costa and Silvestre T. Pinho
%                                      Imperial College London, Aeronautics
% **************************************************************************


function [fel] = ls_fquad4( ndim, ndofel, ngpt, xgp, wgp, nnodel, xel, phiel, h, vel, dt, dtau, rinit )

    % **************************************************************************
    %
    %   INPUTS:
    %       ndim - number of spatial dimensions
    %       ndofel - number of total dof per element
    %       ngpt -  number of integration points
    %       xgp - matrix of coordinates for the integration points
    %       wgp - matrix of weights for the integration points
    %       nnodel - number of total nodes per element
    %       xel - matrix of nodal coordinates
    %       phiel - vector of nodal LS values
    %       h - minimum mesh size
    %       vel - vector of nodal velocities
    %       dt - time step
    %       dtau - reinitialisation time step
    %       rinit - boolean value to toggle reinitialisation
    %
    %
    %   OUTPUTS:
    %       fel - elemental force vector
    %
    % **************************************************************************

    % VARIABLE INITIALISATION
        fel = sparse(zeros(ndofel,1));
        
    % LOOP THROUGH INTEGRATION POINTS
        for igp = 1:ngpt
            % coordinates of integration points
                xi = xgp(igp,1);
                eta = xgp(igp,2);
                
            % shape functions
                [n, dn] = fem_sfquad4( ndim, nnodel, xi, eta ) ;
                
            % jacobian
                J = dn * xel ; detJ = det(J) ; dnxy = J\dn ;
                
            % values from previous iteration
                phi = n*phiel ; dphi = dnxy*phiel ; normdphi = norm(dphi) ; 
                
            % normal vector
                normal = -dphi / normdphi ;
                
            % s function
                s = phi / sqrt(phi^2 + h^2 * normdphi^2) ;
                
            % w function
                w = s * normal ;
                
            % velocity vector
                vn = n * vel;
                v = vn * normal ;

            % elemental vector computation
                r3 = transpose(n) * phi * detJ * wgp(igp) ;
                
                % regular
                    if rinit == 0
                        beta1 = max(0, 1/( 2 * sqrt( dt^(-2) + norm(J\v)^2 ) )) ;
                        W = n + beta1 * transpose(v) * dnxy ;

                        r1 = transpose(W) * vn * normdphi * detJ * wgp(igp) ;

                        fel = fel + dt*r1 + r3;
                    
                % reinitialisation
                    else
                        beta2 = max(0, 1/( 2 * sqrt( dtau^(-2) + norm(J\w)^2 ) )) ;
                        Wt = n + beta2 * transpose(w) * dnxy ;

                        r1 = transpose(Wt) * s * detJ * wgp(igp) ;
                        r2 = transpose(Wt) * (transpose(w) * dphi) * detJ * wgp(igp) ;

                        fel = fel + dtau*(r1 - r2) + r3;
                    end 
        end
end