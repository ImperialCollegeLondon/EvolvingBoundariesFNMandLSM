% **************************************************************************
%  MODELLING EVOLVING BOUNDARIES USING FNM AND LSM (April 2020)
% **************************************************************************
%
%
%   ls_cond.m
%
%
%   FUNCTIONS:
%       - ls_cond - calculates the LS conditioning number, a measure of
%                   the deviation from the signed distance function
%
%
% **************************************************************************
%                   developed by Rui O.S.S. da Costa and Silvestre T. Pinho
%                                      Imperial College London, Aeronautics
% **************************************************************************


function [cond] = ls_cond( s_mesh, s_geo, s_param, s_sol )

    % **************************************************************************
    %
    %   INPUTS:
    %       s_mesh - structure arry with mesh data
    %       s_geo - structure arry with geometric data
    %       s_param - structure arry with generic parameters
    %       s_sol - structure array with solution data
    %
    %   OUTPUTS:
    %       cond - conditioning number
    %
    % **************************************************************************

    % VARIABLE INITIALISATION
        nel = s_mesh.nel;
        con = s_mesh.con;
        x = s_mesh.x;
        ngptquad = s_mesh.ngptquad;
        ngpttri = s_mesh.ngpttri;
        xgpquad = s_mesh.xgpquad;
        xgptri = s_mesh.xgptri;
        ndim = s_mesh.ndim;
        u = s_sol.u;
        elstat = s_mesh.elstat ;
        gradv = zeros(20 * nel, 1);
        numgp = 1;
        cond = 0;
    
    % LOOP THROUGH ELEMENTS
        for iel = 1:nel

            % loop through sub-elements
                nsub = fnm_partitioninfo(elstat(iel,1)) ;
                for isub = 1:nsub
                    
                    % local indices
                        [scon] = fnm_subcon( elstat(iel,1), isub ) ;
                        nnodesub = size(scon,2) ;

                    % coordinates and LS values
                        xel = x(con(iel,scon), :) ;
                        phiel = u(con(iel,scon), 1) ;

                    if nnodesub == 4 % quad sub-element
                        
                        % loop through integration points
                            for igp = 1:ngptquad
                                
                                % coordinates
                                    xi = xgpquad(igp,1);
                                    eta = xgpquad(igp,2);
                                    
                                % shape functions
                                    [n, dn] = fem_sf( 'quad4', ndim, nnodesub, xi, eta ) ;
                                    
                                % jacobian
                                    J = dn * xel ; dnxy = J\dn;

                                % norm of the gradient of the LS field
                                    gradv(numgp,1) = norm( dnxy*phiel );

                                % increment counter
                                    numgp = numgp + 1;
                            end
                    elseif nnodesub == 3 % tri sub-element
                        
                        % loop through integration points
                            for igp = 1:ngpttri
                                
                                % coordinates
                                    xi = xgptri(igp,1);
                                    eta = xgptri(igp,2);
                                
                                % shape functions
                                    [n, dn] = fem_sf( 'tri3', ndim, nnodesub, xi, eta ) ;
                                
                                % jacobian
                                    J = dn * xel ; dnxy = J\dn;

                                % condition number of the level set
                                    gradv(numgp,1) = norm( dnxy*phiel );

                                % increment counter
                                    numgp = numgp + 1;
                            end
                    end
                end
        end
    
    % DROP VECTOR ENTRIES THAT ARE EMPTY
        if numgp < size(gradv,1)
            gradv(numgp:end,:) = [];
        end
    
    % COMPUTE THE CONDITIONING NUMBER
        for igp = 1:size(gradv,1)
            cond = cond + (gradv(igp,1) - 1)^2 ;
        end
        cond = sqrt(cond/size(gradv,1)) ;
end