%% TOPOPT
%
%   INITIALIZATION OF THE LEVEL SET FIELD

function [p, l_u] = topt_initphi( s_mesh, s_geo, s_param, s_sol, nx, ny)

    Lx = s_geo.Lx;
    Ly = s_geo.Ly;
    divx = s_mesh.divx;
    divy = s_mesh.divy;
    nel = s_mesh.nel;
    nnode = s_mesh.nnode;
    x = s_mesh.x;
    l_u = s_sol.l_u;
    l_dofel = s_mesh.l_dofel;
    ndim = s_mesh.ndim;
    nnodel = s_mesh.nnodel;
    con = s_mesh.con;
    l_ndofel = s_mesh.l_ndofel;
    ngpt = s_mesh.ngpt;
    xgp = s_mesh.xgp;
    wgp = s_mesh.wgp;
    c = s_param.c;
    h = s_param.h;
    dtau = s_param.dtau;
    ubc = s_param.ubc;
    
    [p] = createholes( Lx, Ly, divx, divy, nel, nnode, x, ubc, nx, ny ) ;

    l_u(:, 1) = p;
end

%%
function [p] = createholes( Lx, Ly, divx, divy, nel, nnode, x, ubc, nhx, nhy )
% initializes the density function and level set field

    dx = Lx/divx ;
    dy = Ly/divy ;
    d = min([dx dy]) ;
    
    p = zeros(nnode, 1);
    
    % center points and radiuses for initial holes
    if nhx == 0 || nhy == 0
        holes = [Lx/2, Ly/2];
        r = 0.2;
    else
        nx = nhx ;
        ny = nhy ;

        xx = linspace(0, Lx, nx);
        yy = linspace(0, Ly, ny);
        count = 1;
        holes = zeros(17,3);
        r = 0.2 * sqrt( Lx*Ly / pi ) ;
        for iy = 1:ny
            for ix = 1:nx
                if mod(iy,2) == 0 && mod(ix,2) ~= 0
                    holes(count,:) = [xx(ix), yy(iy), r];
                    count = count + 1;
                elseif mod(iy,2) ~= 0 && mod(ix,2) == 0
                    holes(count,:) = [xx(ix), yy(iy), r];
                    count = count + 1;
                end
            end
        end
    end
    
    if ubc == 3
        indx = find( holes(:,1) == Lx/2 & holes(:,2) == 0 );
        holes(indx,:) = [];
    end
    
             
    for ind = 1:nnode
        xholes = holes(:,1:2) ;
        xnode = repmat(x(ind,:), size(xholes,1), 1) ;
        d1 = xnode - xholes ;
        d1 = d1.^2 ;
        d1 = sum(d1,2) ;
        d1 = d1.^0.5 ;
        d1 = min(d1);
        
        if d1 < r 
            p(ind,1) = d1-r;
        elseif any(d1 == r)
            p(ind,1) = 0;
        else
            p(ind,1) = d1-r;
        end
    end
end