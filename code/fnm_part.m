% **************************************************************************
%  MODELLING EVOLVING BOUNDARIES USING FNM AND LSM (April 2020)
% **************************************************************************
%
%
%   fnm_part.m
%
%
%   FUNCTIONS:
%       - fnm_part - creates the element partitioning based on the LS field
%                    and assigns the element status, activates the correct
%                    combination of floating and real nodes and the list of
%                    elements in and out of the domain is updated
%
%
% **************************************************************************
%                   developed by Rui O.S.S. da Costa and Silvestre T. Pinho
%                                      Imperial College London, Aeronautics
% **************************************************************************


function [fbnd, elstat, x, fadof, u, elsin, elsout, elsbnd] = fnm_part( s_mesh, s_geo, s_param, s_sol )

    % **************************************************************************
    %
    %   INPUTS:
    %       s_mesh - structure arry with mesh data
    %       s_geo - structure arry with geometric data
    %       s_param - structure arry with generic parameters
    %       s_sol - structure array with solution data
    %
    %   OUTPUTS:
    %       fbnd - list of global indicies of floating nodes on
    %              the evolving boundary
    %       elstat - vector containing the flag defining the
    %                element partitioning status
    %       x - matrix of nodal coordinates ( changed to account
    %           for the activation of floating nodes )
    %       fadof - vector that lists global indicies of all
    %               active floating DoFs
    %       u - global solution vector ( changed to account for the activation
    %           floating nodes in the element area )
    %       elsin - vector that lists all the elements inside domain
    %               \Omega at the current iteration
    %       elsout - vector that lists all the elements outside domain
    %                \Omega at the current iteration
    %       elsbnd - vector that lists all the elements close to the boundary
    %                that have been partitioned at the current iteration
    %
    % **************************************************************************

    % VARIABLE INITIALISATION
        nel = s_mesh.nel ;
        x = s_mesh.x;
        con = s_mesh.con;
        u = s_sol.u ;
        dofel = s_mesh.dofel;
        ndofel = s_mesh.ndofel;
        nfloatdofel = s_mesh.nfloatdofel ;
        h = s_param.h ;
        elstat = zeros(nel, 1);
        fbnd = [];
        fadof = [];
        dofnode = s_mesh.dofnode;
        elsoutdom = s_mesh.elsoutdom;
        els = setdiff((1:nel)',elsoutdom);
        tr = 0.05; % treshold for minimum edge length
        elsin = [];
        elsout = [];
        elsbnd = [];
    
    % LOOP THROUGH ELEMENTS INSIDE THE DOMAIN
        for iel = els'
            
            % elemental LS field
                phiel = u( dofel( iel, 1:(ndofel - nfloatdofel) ), 1);

            % check for elements in which LS values of different signs exist (intersected by the boundary)
                if any(phiel < 0) && any(phiel > 0)

                    % update element list
                        elsbnd = [elsbnd; iel];
                        
                    % check with local node indices are negative
                        indneg = find(phiel <= 0) ;
                        
                    % intersection case: bottom left corner 
                        if isequal(indneg, 1) || isequal(indneg, [2;3;4])

                            % coordinates and ls values
                                x_1 = x(con(iel,1), 1);
                                x_2 = x(con(iel,2), 1);
                                y_1 = x(con(iel,1), 2);
                                y_4 = x(con(iel,4), 2);
                                phi_1 = u(con(iel,1), 1);
                                phi_2 = u(con(iel,2), 1);
                                phi_4 = u(con(iel,4), 1);

                            % intersection coordinates (new floating node coordinates)
                                x_5 = x_1 - (x_2 - x_1) * ( phi_1 / (phi_2 - phi_1) ) ;
                                y_5 = y_1 ;

                                x_8 = x_1 ;
                                y_8 = y_1 - (y_4 - y_1) * ( phi_1 / (phi_4 - phi_1) ) ;

                                x_9 = x_1 + 0.5 * h ;
                                y_9 = y_1 + 0.5 * h ;

                            % calculating edge lengths
                                e1 = abs(x_5 - x_1)/h ;
                                e2 = abs(y_8 - y_1)/h ;

                            % check whether edges respect the threshold
                            % if not set the position of the floating node to 
                            % maximum possible within the threshold
                                if e1 < tr
                                    x_5 = x_1 + tr*h ;
                                elseif e1 > 1-tr
                                    x_5 = x_1 + (1-tr)*h ;
                                end
                                if e2 < tr
                                    y_8 = y_1 + tr*h ;
                                elseif e2 > 1-tr
                                    y_8 = y_1 + (1-tr)*h ;
                                end

                            % updating floating node coordinates and element status
                                elstat(iel,1) = 11 ;
                                x(con(iel,5), :) = [x_5 y_5];
                                x(con(iel,8), :) = [x_8 y_8];
                                x(con(iel,9), :) = [x_9 y_9];
                                
                            % updating floating node lists
                                fbnd = [fbnd; con(iel,5); con(iel,8)] ;
                                fadof = [fadof; dofnode(con(iel,5),1)' ; dofnode(con(iel,8),1)'; dofnode(con(iel,9),1)'];
                            
                            % updating LS value for the central floating node
                                u(con(iel,9), 1) = 1/2 * ( phi_4 + phi_2 ) ;

                                
                                
                    % intersection case: top left corner 
                        elseif isequal(indneg, 4) || isequal(indneg, [1;2;3])

                            % coordinates and ls values
                                x_4 = x(con(iel,4), 1);
                                x_3 = x(con(iel,3), 1);
                                y_4 = x(con(iel,4), 2);
                                y_1 = x(con(iel,1), 2);
                                phi_4 = u(con(iel,4), 1);
                                phi_3 = u(con(iel,3), 1);
                                phi_1 = u(con(iel,1), 1);

                            % intersection coordinates (new floating node coordinates)
                                x_7 = x_4 - (x_3 - x_4) * ( phi_4 / (phi_3 - phi_4) ) ;
                                y_7 = y_4 ;

                                x_8 = x_4 ;
                                y_8 = y_4 - (y_1 - y_4) * ( phi_4 / (phi_1 - phi_4) ) ;

                                x_9 = x_4 + 0.5 * h ;
                                y_9 = y_1 + 0.5 * h ;

                            % calculating edges lengths
                                e1 = abs(x_7 - x_4)/h ;
                                e2 = abs(y_8 - y_4)/h ;

                            % check whether edges respect the threshold
                            % if not set the position of the floating node to 
                            % maximum possible within the threshold
                                if e1 < tr
                                    x_7 = x_4 + tr*h ;
                                elseif e1 > 1-tr
                                    x_7 = x_4 + (1-tr)*h ;
                                end
                                if e2 < tr
                                    y_8 = y_4 - tr*h ;
                                elseif e2 > 1-tr
                                    y_8 = y_4 - (1-tr)*h ;
                                end
                            
                            % updating floating node coordinates and element status
                                elstat(iel,1) = 12 ;
                                x(con(iel,7), :) = [x_7 y_7];
                                x(con(iel,8), :) = [x_8 y_8];
                                x(con(iel,9), :) = [x_9 y_9];
                            
                            % updating floating node lists
                                fbnd = [fbnd; con(iel,7); con(iel,8)] ;
                                fadof = [fadof; dofnode(con(iel,7),1)'; dofnode(con(iel,8),1)'; dofnode(con(iel,9),1)'];

                            % updating LS value for the central floating node
                                u(con(iel,9), 1) = 1/2 * ( phi_1 + phi_3 ) ;

                                
                                
                    % intersection case: top right corner 
                        elseif isequal(indneg, 3) || isequal(indneg, [1;2;4])

                            % coordinates and ls values
                                x_3 = x(con(iel,3), 1);
                                x_4 = x(con(iel,4), 1);
                                y_3 = x(con(iel,3), 2);
                                y_2 = x(con(iel,2), 2);
                                y_4 = x(con(iel,4), 2);
                                phi_3 = u(con(iel,3), 1);
                                phi_4 = u(con(iel,4), 1);
                                phi_2 = u(con(iel,2), 1);

                            % intersection coordinates (new floating node coordinates)
                                x_7 = x_3 - (x_4 - x_3) * ( phi_3 / (phi_4 - phi_3) ) ;
                                y_7 = y_3 ;

                                x_6 = x_3 ;
                                y_6 = y_3 - (y_2 - y_3) * ( phi_3 / (phi_2 - phi_3) ) ;

                                x_9 = x_4 + 0.5 * h ;
                                y_9 = y_4 - 0.5 * h ;

                            % calculating edges lengths
                                e1 = abs(x_7 - x_3)/h ;
                                e2 = abs(y_6 - y_3)/h ;

                            % check whether edges respect the threshold
                            % if not set the position of the floating node to 
                            % maximum possible within the threshold
                                if e1 < tr
                                    x_7 = x_3 - tr*h ;
                                elseif e1 > 1-tr
                                    x_7 = x_3 - (1-tr)*h ;
                                end
                                if e2 < tr
                                    y_6 = y_3 - tr*h ;
                                elseif e2 > 1-tr
                                    y_6 = y_3 - (1-tr)*h ;
                                end

                            % updating floating node coordinates and element status
                                elstat(iel,1) = 13 ;
                                x(con(iel,6), :) = [x_6 y_6];
                                x(con(iel,7), :) = [x_7 y_7];
                                x(con(iel,9), :) = [x_9 y_9];
                            
                            % updating floating node lists
                                fbnd = [fbnd; con(iel,6); con(iel,7)] ;
                                fadof = [fadof; dofnode(con(iel,6),1)'; dofnode(con(iel,7),1)'; dofnode(con(iel,9),1)'];

                            % updating LS value for the central floating node
                                u(con(iel,9), 1) = 1/2 * ( phi_4 + phi_2 ) ;

                            
                            
                    % intersection case: bottom right corner 
                        elseif isequal(indneg, 2) || isequal(indneg, [1;3;4])

                            % coordinates and ls values
                                x_2 = x(con(iel,2), 1);
                                x_1 = x(con(iel,1), 1);
                                y_2 = x(con(iel,2), 2);
                                y_3 = x(con(iel,3), 2);
                                y_1 = x(con(iel,1), 2);
                                phi_2 = u(con(iel,2), 1);
                                phi_1 = u(con(iel,1), 1);
                                phi_3 = u(con(iel,3), 1);

                            % intersection coordinates (new floating node coordinates)
                                x_5 = x_2 - (x_1 - x_2) * ( phi_2 / (phi_1 - phi_2) ) ;
                                y_5 = y_2 ;

                                x_6 = x_2 ;
                                y_6 = y_2 - (y_3 - y_2) * ( phi_2 / (phi_3 - phi_2) ) ;

                                x_9 = x_1 + 0.5 * h ;
                                y_9 = y_1 + 0.5 * h ;

                            % calculating edges lengths
                                e1 = abs(x_5 - x_2)/h ;
                                e2 = abs(y_6 - y_2)/h ;

                            % check whether edges respect the threshold
                            % if not set the position of the floating node to 
                            % maximum possible within the threshold
                                if e1 < tr
                                    x_5 = x_2 - tr*h ;
                                elseif e1 > 1-tr
                                    x_5 = x_2 - (1-tr)*h ;
                                end
                                if e2 < tr
                                    y_6 = y_2 + tr*h ;
                                elseif e2 > 1-tr
                                    y_6 = y_2 + (1-tr)*h ;
                                end

                            % updating floating node coordinates and element status
                                elstat(iel,1) = 14 ;
                                x(con(iel,5), :) = [x_5 y_5];
                                x(con(iel,6), :) = [x_6 y_6];
                                x(con(iel,9), :) = [x_9 y_9];
                            
                            % updating floating node lists
                                fbnd = [fbnd; con(iel,5); con(iel,6)] ;
                                fadof = [fadof; dofnode(con(iel,5),1)'; dofnode(con(iel,6),1)'; dofnode(con(iel,9),1)'];
                            
                            % updating LS value for the central floating node
                                u(con(iel,9), 1) = 1/2 * ( phi_1 + phi_3 ) ;

                            
                            
                    % intersection case: through element vertical 
                        elseif isequal(indneg, [1;4]) || isequal(indneg, [2;3])

                            % coordinates and ls values
                                x_1 = x(con(iel,1), 1);
                                x_2 = x(con(iel,2), 1);
                                x_3 = x(con(iel,3), 1);
                                x_4 = x(con(iel,4), 1);
                                y_1 = x(con(iel,1), 2);
                                y_3 = x(con(iel,3), 2);
                                phi_1 = u(con(iel,1), 1);
                                phi_2 = u(con(iel,2), 1);
                                phi_3 = u(con(iel,3), 1);
                                phi_4 = u(con(iel,4), 1);

                            % intersection coordinates (new floating node coordinates)
                                x_5 = x_1 - (x_2 - x_1) * ( phi_1 / (phi_2 - phi_1) ) ;
                                y_5 = y_1 ;

                                x_7 = x_3 - (x_4 - x_3) * ( phi_3 / (phi_4 - phi_3) ) ;
                                y_7 = y_3 ;

                            % calculating edges lengths
                                e1 = abs(x_5 - x_1)/h ;
                                e2 = abs(x_7 - x_4)/h ;

                            % check whether edges respect the threshold
                            % if not set the position of the floating node to 
                            % maximum possible within the threshold
                                if e1 < tr
                                    x_5 = x_1 + tr*h ;
                                elseif e1 > 1-tr
                                    x_5 = x_1 + (1-tr)*h ;
                                end
                                if e2 < tr
                                    x_7 = x_4 + tr*h ;
                                elseif e2 > 1-tr
                                    x_7 = x_4 + (1-tr)*h ;
                                end

                            % updating floating node coordinates and element status
                                elstat(iel,1) = 21 ;
                                x(con(iel,5), :) = [x_5 y_5];
                                x(con(iel,7), :) = [x_7 y_7];
                            
                            % updating floating node lists
                                fbnd = [fbnd; con(iel,5); con(iel,7)] ;
                                fadof = [fadof; dofnode(con(iel,5),1)'; dofnode(con(iel,7),1)'];


                                                       
                    % intersection case: through element horizontal
                        elseif isequal(indneg, [1;2]) || isequal(indneg, [3;4])
                            
                            % coordinates and ls values
                                y_1 = x(con(iel,1), 2);
                                y_2 = x(con(iel,2), 2);
                                y_3 = x(con(iel,3), 2);
                                y_4 = x(con(iel,4), 2);
                                x_1 = x(con(iel,1), 1);
                                x_3 = x(con(iel,3), 1);
                                phi_1 = u(con(iel,1), 1);
                                phi_2 = u(con(iel,2), 1);
                                phi_3 = u(con(iel,3), 1);
                                phi_4 = u(con(iel,4), 1);

                            % intersection coordinates (new floating node coordinates)
                                x_6 = x_3 ;
                                y_6 = y_3 - (y_2 - y_3) * ( phi_3 / (phi_2 - phi_3) ) ;

                                x_8 = x_1 ;
                                y_8 = y_1 - (y_4 - y_1) * ( phi_1 / (phi_4 - phi_1) ) ;

                            % calculating edges lengths
                                e1 = abs(y_8 - y_1)/h ;
                                e2 = abs(y_6 - y_2)/h ;

                            % check whether edges respect the threshold
                            % if not set the position of the floating node to 
                            % maximum possible within the threshold
                                if e1 < tr
                                    y_8 = y_1 + tr*h ;
                                elseif e1 > 1-tr
                                    y_8 = y_1 + (1-tr)*h ;
                                end
                                if e2 < tr
                                    y_6 = y_2 + tr*h ;
                                elseif e2 > 1-tr
                                    y_6 = y_2 + (1-tr)*h ;
                                end
    
                            % updating floating node coordinates and element status
                                elstat(iel,1) = 22 ;
                                x(con(iel,6), :) = [x_6 y_6];
                                x(con(iel,8), :) = [x_8 y_8];
                            
                            % updating floating node lists
                                fbnd = [fbnd; con(iel,6); con(iel,8)] ;
                                fadof = [fadof; dofnode(con(iel,6),1)'; dofnode(con(iel,8),1)'];


                            
                    % intersection case: double corner, bottom left and top right
                        elseif isequal(indneg, [1;3]) 

                            % coordinates and ls values
                                x_1 = x(con(iel,1), 1);
                                x_2 = x(con(iel,2), 1);
                                x_3 = x(con(iel,3), 1);
                                x_4 = x(con(iel,4), 1);
                                y_1 = x(con(iel,1), 2);
                                y_2 = x(con(iel,2), 2);
                                y_3 = x(con(iel,3), 2);
                                y_4 = x(con(iel,4), 2);
                                phi_1 = u(con(iel,1), 1);
                                phi_2 = u(con(iel,2), 1);
                                phi_3 = u(con(iel,3), 1);
                                phi_4 = u(con(iel,4), 1);

                            % intersection coordinates (new floating node coordinates)
                                x_5 = x_1 - (x_2 - x_1) * ( phi_1 / (phi_2 - phi_1) ) ;
                                y_5 = y_1 ;

                                x_6 = x_3 ;
                                y_6 = y_3 - (y_2 - y_3) * ( phi_3 / (phi_2 - phi_3) ) ;

                                x_7 = x_3 - (x_4 - x_3) * ( phi_3 / (phi_4 - phi_3) ) ;
                                y_7 = y_3 ;

                                x_8 = x_1 ;
                                y_8 = y_1 - (y_4 - y_1) * ( phi_1 / (phi_4 - phi_1) ) ;

                                x_9 = x_1 + 0.5*h ;
                                y_9 = y_1 + 0.5*h ;

                            % calculating edges lengths
                                e1 = abs(x_5 - x_1)/h ;
                                e2 = abs(y_6 - y_2)/h ;
                                e3 = abs(x_7 - x_4)/h ;
                                e4 = abs(y_8 - y_1)/h ;

                            % check whether edges respect the threshold
                            % if not set the position of the floating node to 
                            % maximum possible within the threshold
                                if e1 < tr
                                    x_5 = x_1 + tr*h ;
                                elseif e1 > 1-tr
                                    x_5 = x_1 + (1-tr)*h ;
                                end
                                if e2 < tr
                                    y_6 = y_2 + tr*h ;
                                elseif e2 > 1-tr
                                    y_6 = y_2 + (1-tr)*h ;
                                end
                                if e3 < tr
                                    x_7 = x_4 + tr*h ;
                                elseif e3 > 1-tr
                                    x_7 = x_4 + (1-tr)*h ;
                                end
                                if e4 < tr
                                    y_8 = y_1 + tr*h ;
                                elseif e4 > 1-tr
                                    y_8 = y_1 + (1-tr)*h ;
                                end

                            % updating floating node coordinates and element status
                                elstat(iel,1) = 31 ;
                                x(con(iel,5), :) = [x_5 y_5];
                                x(con(iel,6), :) = [x_6 y_6];
                                x(con(iel,7), :) = [x_7 y_7];
                                x(con(iel,8), :) = [x_8 y_8];
                                x(con(iel,9), :) = [x_9 y_9];
                                
                            % updating floating node lists
                                fbnd = [fbnd; con(iel,5); con(iel,6); con(iel,7); con(iel,8)] ;
                                fadof = [fadof; dofnode(con(iel,5),1)'; dofnode(con(iel,6),1)'; dofnode(con(iel,7),1)'; dofnode(con(iel,8),1)'; dofnode(con(iel,9),1)'];

                            % updating LS value for the central floating node
                                u(con(iel,9), 1) = 1/2 * ( phi_2 + phi_4 ) ;

                            
                            
                    % intersection case: double corner, bottom right and top left
                        elseif isequal(indneg, [2;4])

                            % coordinates and ls values
                                x_1 = x(con(iel,1), 1);
                                x_2 = x(con(iel,2), 1);
                                x_3 = x(con(iel,3), 1);
                                x_4 = x(con(iel,4), 1);
                                y_1 = x(con(iel,1), 2);
                                y_2 = x(con(iel,2), 2);
                                y_3 = x(con(iel,3), 2);
                                y_4 = x(con(iel,4), 2);
                                phi_1 = u(con(iel,1), 1);
                                phi_2 = u(con(iel,2), 1);
                                phi_3 = u(con(iel,3), 1);
                                phi_4 = u(con(iel,4), 1);

                            % intersection coordinates (new floating node coordinates)
                                x_5 = x_1 - (x_2 - x_1) * ( phi_1 / (phi_2 - phi_1) ) ;
                                y_5 = y_1 ;

                                x_6 = x_3 ;
                                y_6 = y_3 - (y_2 - y_3) * ( phi_3 / (phi_2 - phi_3) ) ;

                                x_7 = x_3 - (x_4 - x_3) * ( phi_3 / (phi_4 - phi_3) ) ;
                                y_7 = y_3 ;

                                x_8 = x_1 ;
                                y_8 = y_1 - (y_4 - y_1) * ( phi_1 / (phi_4 - phi_1) ) ;

                                x_9 = x_1 + 0.5*h ;
                                y_9 = y_1 + 0.5*h ;

                            % calculating edges lengths
                                e1 = abs(x_5 - x_1)/h ;
                                e2 = abs(y_6 - y_2)/h ;
                                e3 = abs(x_7 - x_4)/h ;
                                e4 = abs(y_8 - y_1)/h ;

                            % check whether edges respect the threshold
                            % if not set the position of the floating node to
                            % maximum possible within the threshold
                                if e1 < tr
                                    x_5 = x_1 + tr*h ;
                                elseif e1 > 1-tr
                                    x_5 = x_1 + (1-tr)*h ;
                                end
                                if e2 < tr
                                    y_6 = y_2 + tr*h ;
                                elseif e2 > 1-tr
                                    y_6 = y_2 + (1-tr)*h ;
                                end
                                if e3 < tr
                                    x_7 = x_4 + tr*h ;
                                elseif e3 > 1-tr
                                    x_7 = x_4 + (1-tr)*h ;
                                end
                                if e4 < tr
                                    y_8 = y_1 + tr*h ;
                                elseif e4 > 1-tr
                                    y_8 = y_1 + (1-tr)*h ;
                                end

                            % updating floating node coordinates and element status
                                elstat(iel,1) = 32 ;
                                x(con(iel,5), :) = [x_5 y_5];
                                x(con(iel,6), :) = [x_6 y_6];
                                x(con(iel,7), :) = [x_7 y_7];
                                x(con(iel,8), :) = [x_8 y_8];
                                x(con(iel,9), :) = [x_9 y_9];
                            
                            % updating floating node lists
                                fbnd = [fbnd; con(iel,5); con(iel,6); con(iel,7); con(iel,8)] ;
                                fadof = [fadof; dofnode(con(iel,5),1)'; dofnode(con(iel,6),1)'; dofnode(con(iel,7),1)'; dofnode(con(iel,8),1)'; dofnode(con(iel,9),1)'];

                            % updating LS value for the central floating node
                                u(con(iel,9), 1) = 1/2 * ( phi_1 + phi_3 ) ;
                        end
                         
            % check for elements with only negative LS values (in the void region)
                elseif any(phiel < 0)

                    % update element list 
                        elsout = [elsout; iel];
            
            % the remainder are elements with only positive LS values (un-partitioned in the solid region)
                else

                    % update element list 
                        elsin = [elsin; iel];
                end
        end
    
    % REMOVING POSSIBLE REPEATED INDICES    
        fbnd = unique(fbnd);
        fadof = unique(fadof); 
end