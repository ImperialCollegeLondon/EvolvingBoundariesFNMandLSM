%% TOPOPT
%
%   PLOTTING

function [] = topt_plot(n, s_mesh, s_geo, s_param, s_sol, s_hist)

    nel = s_mesh.nel;
    con = s_mesh.con;
    x = s_mesh.x;
    Lx = s_geo.Lx;
    Ly = s_geo.Ly;
    divx = s_mesh.divx;
    divy = s_mesh.divy;
    u = s_sol.u ;
    l_u = s_sol.l_u ;
    dofel = s_mesh.dofel;
    l_dofel = s_mesh.l_dofel;
    nnodel = s_mesh.nnodel ;
    nnode = s_mesh.nnode ;
    p = s_param.p ;
    rho = s_param.rho;
    G = s_param.G;
    ndofel = s_mesh.ndofel;
    ngpt = s_mesh.ngpt;
    objhist = s_hist.objhist;
    volhist = s_hist.volhist;
    vn = s_sol.vn;
    xgp = s_mesh.xgp;
    ndim = s_mesh.ndim;
    elU = s_param.elU;
    ndU = s_param.ndU;
    sig = s_sol.sig;

    if n == 1 % mesh plotting
%         pltmesh( nel, con, x, divx, divy, Lx, Ly ) ;
    elseif n == 2 % deformed plotting
        pltdeformed ( nel, con, x, u, dofel, nnodel, Lx, Ly, divx, divy, nnode )
    elseif n == 3 % initial config plotting
        pltinitconfig( nel, con, x , divx, divy, Lx, Ly, p, rho, G )
    elseif n == 4 % level set plotting
        pltlevelset( nel, l_u, vn, l_dofel, con, x, Lx, Ly, divx, divy, xgp, ndim, nnodel, ngpt )
    elseif n == 6 % plot obj and vol
        pltobjvol( objhist, volhist )
    elseif n == 5 % plot vel
        pltvel( nel, vn, elU, ndU, con, x, Lx, Ly, divx, divy)
    elseif n == 7 
        pltstr( nel, con, x, sig, rho,  Lx, Ly, divx, divy )
    end

end

%%
function [] = pltmesh( nel, con, x , divx, divy, Lx, Ly )
% plots the mesh in the undeformed configuration

    figure(101) 
    clf
    hold on
    for iel = 1:nel
        nds = con(iel,:);
        xcoord = x(nds,1);
        ycoord = x(nds,2);
        patch(xcoord,ycoord,'blue','FaceAlpha',.2) ;
        plot(xcoord,ycoord,'ko') ;
    end
    hold off
    if Lx > Ly
        s = 1200*[1 Ly/Lx];
    else
        s = 1200*[Lx/Ly 1];
    end
    set(gcf, 'Position', [100 100 s]) 
    box off
    dx = Lx/divx;
    dy = Ly/divy;
    pbaspect([(Lx + 2*dx)/(Ly + 2*dy) 1 1]) ;
    xlim([-dx Lx+dx])
    ylim([-dy Ly+dy])

end

%%
function [] = pltdeformed( nel, con, x, u, dofel, nnodel, Lx, Ly, divx, divy, nnode )
% plots the mesh in the deformed configuration

    figure(201)
    clf
    dx = Lx/divx;
    dy = Ly/divy;
    delta = max( max( abs( u(:,1) ) ) );
    f = 0.8 * ( min([dx dy]) / delta ) ;
    hold on
    for iel = 1:nel
        nds = con(iel,:);
        xcoord = x(nds,1);
        ycoord = x(nds,2);
        dof = u( dofel(iel,:), 1);
        uel = zeros(4,1);
        vel = zeros(4,1);
        for ind = 1:nnodel
            uel(ind) = dof( 2*ind - 1, 1);
            vel(ind) = dof( 2*ind, 1);
        end
%         patch(xcoord,ycoord,'blue','FaceAlpha',0) ;
        patch(xcoord + f*uel,ycoord + f*vel, uel) ;
%         plot(xcoord,ycoord,'ko') ;
%         plot(xcoord + f*uel,ycoord + f*vel,'ko') ;
    end
    hold off
    if Lx > Ly
        s = 1200*[1 Ly/Lx];
    else
        s = 1200*[Lx/Ly 1];
    end
    set(gcf, 'Position', [500 100 s]) 
    box off
    uu = zeros(nnode, 1);
    vv = zeros(nnode, 1);
    for ind = 1:nnode
        uu(ind) = u( 2*ind - 1, 1) ;
        vv(ind) = u( 2*ind, 1) ;
    end
    cb = colorbar('XTick', linspace(min(uu(:,1)), max(uu(:,1)), 10) );
    title('u')
    pbaspect([(Lx + 2*dx)/(Ly + 2*dy) 1 1]) ;
    xlim([-dx Lx+dx])
    ylim([-dy Ly+dy])
    
    
    figure(202)
    clf
    dx = Lx/divx;
    dy = Ly/divy;
    delta = max( max( abs( u(:,1) ) ) );
    f = 0.8 * ( min([dx dy]) / delta ) ;
    hold on
    for iel = 1:nel
        nds = con(iel,:);
        xcoord = x(nds,1);
        ycoord = x(nds,2);
        dof = u( dofel(iel,:), 1);
        uel = zeros(4,1);
        vel = zeros(4,1);
        for ind = 1:nnodel
            uel(ind) = dof( 2*ind - 1, 1);
            vel(ind) = dof( 2*ind, 1);
        end
%         patch(xcoord,ycoord,'blue','FaceAlpha',0) ;
        patch(xcoord + f*uel,ycoord + f*vel, vel) ;
%         plot(xcoord,ycoord,'ko') ;
%         plot(xcoord + f*uel,ycoord + f*vel,'ko') ;
    end
    hold off
    if Lx > Ly
        s = 1200*[1 Ly/Lx];
    else
        s = 1200*[Lx/Ly 1];
    end
    set(gcf, 'Position', [1000 100 s]) 
    box off
    cb = colorbar('XTick', linspace(min(vv(:,1)), max(vv(:,1)), 10) );
    title('v')
    pbaspect([(Lx + 2*dx)/(Ly + 2*dy) 1 1]) ;
    xlim([-dx Lx+dx])
    ylim([-dy Ly+dy])
    


end

%%
function [] = pltinitconfig( nel, con, x , divx, divy, Lx, Ly, p, rho, G )
% plots the initial domain configuration

%     figure(301) 
%     clf
%     hold on
%     for iel = 1:nel
%         nds = con(iel,:);
%         xcoord = x(nds,1);
%         ycoord = x(nds,2);
%         patch(xcoord, ycoord, -p(nds,1), 'FaceAlpha', .5) ;
%     end
%     hold off
%     if Lx > Ly
%         s = 1200*[1 Ly/Lx];
%     else
%         s = 1200*[Lx/Ly 1];
%     end
%     set(gcf, 'Position', [100 300 s]) 
%     box off
%     colorMap = gray(2);
%     colormap(colorMap);
%     dx = Lx/divx;
%     dy = Ly/divy;
%     pbaspect([(Lx + 2*dx)/(Ly + 2*dy) 1 1]) ;
%     xlim([-dx Lx+dx])
%     ylim([-dy Ly+dy])
    
    
    figure(302)
    clf
    hold on
    for iel = 1:nel
        if rho(iel,1) == 1
            nds = con(iel,:);
            xcoord = x(nds,1);
            ycoord = x(nds,2);
            patch(xcoord, ycoord, 'black') ;
        end
    end
    hold off
    if Lx > Ly
        s = 800*[1 Ly/Lx];
    else
        s = 800*[Lx/Ly 1];
    end
    set(gcf, 'Position', [300 300 s]) 
    box off
    title('Numerical model')
    xlabel('x')
    ylabel('y')
    dx = Lx/divx;
    dy = Ly/divy;
    pbaspect([(Lx + 2*dx)/(Ly + 2*dy) 1 1]) ;
    xlim([-dx Lx+dx])
    ylim([-dy Ly+dy])
    
    
%     figure(303)
%     clf
%     hold on
%     for iel = 1:nel
%         nds = con(iel,:);
%         xcoord = x(nds,1);
%         ycoord = x(nds,2);
%         patch(xcoord, ycoord, 'black', 'FaceAlpha', 0.2) ;
%     end
%     for iint = 1:size(G,1)
%         xpt = G(iint,:) * x ;
%         plot( xpt(1,1), xpt(1,2), 'ro')
%     end
%     hold off
%     if Lx > Ly
%         s = 1200*[1 Ly/Lx];
%     else
%         s = 1200*[Lx/Ly 1];
%     end
%     set(gcf, 'Position', [600 300 s]) 
%     box off
%     dx = Lx/divx;
%     dy = Ly/divy;
%     pbaspect([(Lx + 2*dx)/(Ly + 2*dy) 1 1]) ;
%     xlim([-dx Lx+dx])
%     ylim([-dy Ly+dy])


end

%%
function [] = pltlevelset( nel, l_u, vn, l_dofel, con, x, Lx, Ly, divx, divy, xgp, ndim, nnodel, ngpt )
% plots the level set field

    figure(401) 
    clf
    hold on
%     for iel = 1:nel
%         phiel = l_u(l_dofel(iel,:), 1) ;
%         nds = con(iel,:);
%         xcoord = x(nds,1);
%         ycoord = x(nds,2);
%         patch(xcoord, ycoord, phiel) ;
%     end
    
    phiv = l_u(:,1);
    xun = unique(x(:,1));
    yun = unique(x(:,2));
    [X,Y] = meshgrid(xun,yun) ;
    Z = 0*X;
    nn = 1;
    for irow = 1:size(X,1)
        for icol = 1:size(X,2)
            Z(irow,icol) = phiv(nn,1) ;
            nn = nn + 1;
        end
    end
    
    contourf(X,Y,Z,[0 0], 'LineWidth', 3, 'LineColor', 'k')
    colormap('gray')
    
    hold off
    if Lx > Ly
        s = 800*[1 Ly/Lx];
    else
        s = 800*[Lx/Ly 1];
    end
        set(gcf, 'Position', [250 100 s])
    box off
%     cb = colorbar('XTick', linspace(min(phiv(:,1)), max(phiv(:,1)), 10) );
    dx = Lx/divx;
    dy = Ly/divy;
    title('Zero level set')
    xlabel('x')
    ylabel('y')
    pbaspect([(Lx + 2*dx)/(Ly + 2*dy) 1 1]) ;
    xlim([-dx Lx+dx])
    ylim([-dy Ly+dy])
    
    
    
%     figure(402)
%     clf
%     
%     hold on
%     numgp = 1;
%     vnn = 0.1*vn/max(abs(vn(:,1)));
%     for iel = 1:nel
%         xel = x(con(iel,:),:);
%         phiel = u(dofel(iel,(ndofelm+1):ndofel), 1) ;
%         for igp = 1:ngpt
%             xi = xgp(igp,1);
%             eta = xgp(igp,1);
%             [n, dn] = shapefuncquad( ndim, nnodel, xi, eta ) ;
%             xip = n*xel;
%             % jacobian
%             J = dn * xel ; dnxy = J\dn;
%             % values from previous iteration
%             dphi = dnxy*phiel ; normdphi = norm(dphi) ;
%             % normal vector
%             normal = dphi / normdphi ;
%             
%             v = vnn(numgp,1) * normal;
%             
%             quiver(xip(1),xip(2),v(1),v(2))
%             numgp = numgp + 1;
%         end
%     end
%     hold off
%     if Lx > Ly
%         s = 1200*[1 Ly/Lx];
%     else
%         s = 1200*[Lx/Ly 1];
%     end
%     set(gcf, 'Position', [1000 100 s]) 
%     box off
%     dx = Lx/divx;
%     dy = Ly/divy;
%     pbaspect([(Lx + 2*dx)/(Ly + 2*dy) 1 1]) ;
%     xlim([-dx Lx+dx])
%     ylim([-dy Ly+dy])
    
end

%%
function [] = pltvel( nel, vn, elU, ndU, con, x, Lx, Ly, divx, divy)
% plots the level set field

    figure(501) 
    clf
    hold on
    for iel = 1:nel
        vel = vn(con(iel,:), 1) ;
        nds = con(iel,:);
        xcoord = x(nds,1);
        ycoord = x(nds,2);
        patch(xcoord, ycoord, vel) ;
    end
    hold off
    if Lx > Ly
        s = 800*[1 Ly/Lx];
    else
        s = 800*[Lx/Ly 1];
    end
    set(gcf, 'Position', [100 400 s]) 
    box off
    cb = colorbar('XTick', linspace(min(vn(:,1)), max(vn(:,1)), 10) );
%     colormap(jet(12));
    dx = Lx/divx;
    dy = Ly/divy;
    title('V_n')
    xlabel('x')
    ylabel('y')
    pbaspect([(Lx + 2*dx)/(Ly + 2*dy) 1 1]) ;
    xlim([-dx Lx+dx])
    ylim([-dy Ly+dy])
    
    
%     figure(502) 
%     clf
%     hold on
%     for iel = 1:nel
%         vel = ndU(con(iel,:), 1) ;
%         nds = con(iel,:);
%         xcoord = x(nds,1);
%         ycoord = x(nds,2);
%         patch(xcoord, ycoord, vel) ;
%     end
%     hold off
%     if Lx > Ly
%         s = 1200*[1 Ly/Lx];
%     else
%         s = 1200*[Lx/Ly 1];
%     end
%     set(gcf, 'Position', [-1600 400 s]) 
%     box off
%     cb = colorbar('XTick', linspace(min(vn(:,1)), max(vn(:,1)), 10) );
%     dx = Lx/divx;
%     dy = Ly/divy;
%     title('ndU')
%     pbaspect([(Lx + 2*dx)/(Ly + 2*dy) 1 1]) ;
%     xlim([-dx Lx+dx])
%     ylim([-dy Ly+dy])
    
    
    
%     figure(503) 
%     clf
%     hold on
%     for iel = 1:nel
%         vel = vn(con(iel,:), 1) ;
%         vel = 0*vel + elU(iel,1);
%         nds = con(iel,:);
%         xcoord = x(nds,1);
%         ycoord = x(nds,2);
%         patch(xcoord, ycoord, vel) ;
%     end
%     hold off
%     if Lx > Ly
%         s = 1200*[1 Ly/Lx];
%     else
%         s = 1200*[Lx/Ly 1];
%     end
%     set(gcf, 'Position', [-1600 400 s]) 
%     box off
%     cb = colorbar('XTick', linspace(min(vn(:,1)), max(vn(:,1)), 10) );
%     dx = Lx/divx;
%     dy = Ly/divy;
%     title('elU')
%     pbaspect([(Lx + 2*dx)/(Ly + 2*dy) 1 1]) ;
%     xlim([-dx Lx+dx])
%     ylim([-dy Ly+dy])
%         
    
end

%%
function [] = pltobjvol( objhist, volhist)
% plots objective function and volume
    
    figure(601)
    clf
    plot(1:size(objhist,2), objhist)
    title('Objective function')
    ylabel('Compliance')
    xlabel('Iterations')
    set(gcf, 'Position', [100 200 500 500])
    pbaspect([1 1 1]) ;
    
    
    figure(602)
    clf
    plot(1:size(volhist,2), volhist)
    title('Constraint')
    ylabel('Volume')
    xlabel('Iterations')
    set(gcf, 'Position', [300 200 500 500])
    pbaspect([1 1 1]) ;
    
end

%%
function [] = pltstr( nel, con, x, sig, rho, Lx, Ly, divx, divy )
% plots the mesh in the deformed configuration

    figure(601)
    clf
    hold on
    ndss = [];
    for iel = 1:nel      
        nds = con(iel,:) ;
        xcoord = x(nds,1);
        ycoord = x(nds,2);
        sxx = sig(nds, 1) ;
        if rho(iel,1) == 1
            patch(xcoord,ycoord, sxx, 'EdgeAlpha', 0) ;
            ndss = [ndss; nds];
        end
    end
    hold off
    colormap(jet(16))
    ndss = unique(ndss);
    cb = colorbar('XTick', linspace(min(sig(ndss,1)), max(sig(ndss,1)), 17) );
    title('sxx')
    dx = Lx/divx;
    dy = Ly/divy;
    pbaspect([(Lx + 2*dx)/(Ly + 2*dy) 1 1]) ;
end











