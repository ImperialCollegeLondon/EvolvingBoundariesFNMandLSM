% **************************************************************************
%  MODELLING EVOLVING BOUNDARIES USING FNM AND LSM (April 2020)
% **************************************************************************
%
%
%   ls_reinit.m
%
%
%   FUNCTIONS:
%       - ls_reinit - reinitialises the LS field using the floating nodes
%
%
% **************************************************************************
%                   developed by Rui O.S.S. da Costa and Silvestre T. Pinho
%                                      Imperial College London, Aeronautics
% **************************************************************************

function [u] = ls_reinit( s_mesh, s_geo, s_param, s_sol )

    % **************************************************************************
    %
    %   INPUTS:
    %       s_mesh - structure arry with mesh data
    %       s_geo - structure arry with geometric data
    %       s_param - structure arry with generic parameters
    %       s_sol - structure array with solution data
    %
    %   OUTPUTS:
    %       u - reinitialised LS field
    %
    % **************************************************************************
    
    % VARIABLE INITIALISATION
        u = s_sol.u ;
        fbnd = s_param.fbnd ;
        nnode = s_mesh.nnode - s_mesh.nfloatnode ;
        x = s_mesh.x ;
        nel = s_mesh.nel;
        con = s_mesh.con;
        h = s_param.h;
        xbnd = x(fbnd, :);
    
    % LOOP THROUGH NODES TO UPDATE THE LS VALUE
        for ind = 1:nnode
            
            % nodal coordinates
                xind = x(ind, :);
                
            % distance to boundary nodes
                d = abs(xbnd - xind) ;
                d = d.^2 ;
                d = sum(d,2) ;
                d = sqrt(d) ;
            
            % obtain minimum distance
                d = min(d) ;
                d = d(1);

            % update LS value respecting the sign
                if u(ind,1) > 0
                    u(ind,1) = d ;
                else
                    u(ind,1) = -d ;
                end
        end
    
    % LOOP THROUGH ELEMENTS TO RANDOMISE LS VALUES
    % prevents the LSFEM analysis from reaching a singularity that happens
    % whenever a square element has the same LS value in all nodes
        for iel = 1:nel
            scon = [1 2 3 4];
            ndsr = con(iel,scon) ;
            phiel = u(ndsr,1);
            if abs(phiel(1) - phiel(2)) < h/100 && abs(phiel(3) - phiel(4)) < h/100 && abs(phiel(1)-phiel(3)) < h/100
                randv = (-h/4) + (h/2).*rand(4,1);
                u(ndsr,1) = phiel + randv;
            end
        end
end