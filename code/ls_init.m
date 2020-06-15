% **************************************************************************
%  MODELLING EVOLVING BOUNDARIES USING FNM AND LSM (April 2020)
% **************************************************************************
%
%
%   ls_init.m
%
%
%   FUNCTIONS:
%       - ls_init - initialises the LS field with a pattern of holes
%
%
% **************************************************************************
%                   developed by Rui O.S.S. da Costa and Silvestre T. Pinho
%                                      Imperial College London, Aeronautics
% **************************************************************************

function [p, u] = ls_init( s_mesh, s_geo, s_param, s_sol, nx, ny )

    % **************************************************************************
    %
    %   INPUTS:
    %       s_mesh - structure arry with mesh data
    %       s_geo - structure arry with geometric data
    %       s_param - structure arry with generic parameters
    %       s_sol - structure arry with solution data
    %       nx - number of holes in the x direction
    %       ny - number of holes in the y direction
    %
    %   OUTPUTS:
    %       p - level set field (for sync with mechanical analysis)
    %       u - level set field
    %
    % **************************************************************************

    % VARIABLE INITIALISATION
        Lx = s_geo.Lx;
        Ly = s_geo.Ly;
        nnode = s_mesh.nnode;
        nfloatnode = s_mesh.nfloatnode ;
        x = s_mesh.x;
        u = s_sol.u;
        dofs = s_mesh.dofs;
        ubc = s_param.ubc;
    
    % DEFINE NUMBER OF HOLES AND SIZE
        if nx == 0 && ny == 0
            nrows = 0;
            ncols = 0;
            radius = 0.2 ;
        else 
            nrows = ny;
            ncols = nx;
            if ubc == 1
                aa = 0.15;
            elseif ubc == 2
                aa = 0.1;
            elseif ubc == 3
                aa = 0.1;
            else
                aa = 0.09;
            end
            radius = aa * sqrt( Lx*Ly / pi ) ;
        end

    % CREATE THE HOLE PATTERN
        p = createholes( nrows, ncols, radius, Lx, Ly, nnode - nfloatnode, x, ubc ) ;
        u(dofs, 1) = p;

end


function [p] = createholes( nrows, ncols, radius, Lx, Ly, nnode, x, ubc )
    
    % **************************************************************************
    %
    %   INPUTS:
    %       nrows - number of rows (x)
    %       ncols - number of collumns (y)
    %       radius - radius of each hole
    %       Lx - size of domain in the x direction
    %       Ly - size of domain in the y direction
    %       nnode - number of nodes
    %       x - coordinate matrix
    %       ubc - boundary condition flag
    %
    %   OUTPUTS:
    %       p - level set field
    %
    % **************************************************************************
    
    % VARIABLE INITIALISATION
        p = zeros(nnode,1);
    
    % HOLE IN THE CENTRE
        if nrows == 0 && ncols == 0
            centers = [Lx/2 Ly/2];
            
    % GRID OF HOLES
        else
            
            % number of holes and their centres' coordinates
                nholes = nrows * ncols;
                xc = linspace(0,Lx,ncols);
                yc = linspace(0,Ly,nrows);
                centers = zeros(nholes, 2) ;

            % staggered hole pattern creation
                draw = 0;
                holes2elim = [];
                irow = 1;
                icol = 1;
                for i = 1:nholes

                    if draw == 0
                        holes2elim = [holes2elim i];

                        draw = 1; % draw next one
                    else
                        centers(i,:) = [xc(icol) yc(irow)];

                        draw = 0; % dont draw next one
                    end

                    if icol >= ncols
                        irow = irow + 1;
                        icol = 1;
                    else
                        icol = icol + 1;
                    end
                end
                centers(holes2elim,:) = []; 
        end

    % REMOVE ANY HOLES THAT MAY CONTAIN A LOAD POINT
        if ubc == 1
            dc = sqrt(sum((centers - [Lx, Ly/2]).^2, 2));
            cindx = find(dc < radius);
        elseif ubc == 2
            dc = sqrt(sum((centers - [0, Ly]).^2, 2));
            cindx = find(dc < radius);
        elseif ubc == 3
            dc = sqrt(sum((centers - [Lx, Ly/5]).^2, 2));
            cindx = find(dc < radius);
        elseif ubc == 4
            cindx = [];
        end
    centers(cindx,:) = [];
             
    
    % INITIALISE THE LS FIELD
        for ind = 1:nnode
            xnode = repmat(x(ind,:), size(centers,1), 1) ;
            d1 = xnode - centers ;
            d1 = d1.^2 ;
            d1 = sum(d1,2) ;
            d1 = sqrt(d1) ;
            d1 = min(d1);

            if d1 < radius 
                p(ind,1) = d1 - radius;
            elseif d1 == radius
                p(ind,1) = 0;
            else
                p(ind,1) = d1 - radius;
            end
        end

end