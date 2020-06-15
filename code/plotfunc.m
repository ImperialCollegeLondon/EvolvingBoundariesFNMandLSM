% **************************************************************************
%  MODELLING EVOLVING BOUNDARIES USING FNM AND LSM (April 2020)
% **************************************************************************
%
%
%   plotfunc.m
%
%
% **************************************************************************
%                   developed by Rui O.S.S. da Costa and Silvestre T. Pinho
%                                      Imperial College London, Aeronautics
% **************************************************************************

function [] = plotfunc(n, s_mesh, s_geo, s_param, s_sol)


    % VARIABLE INITIALISATION
        nel = s_mesh.nel;
        con = s_mesh.con;
        x = s_mesh.x;
        Lx = s_geo.Lx;
        Ly = s_geo.Ly;
        divx = s_mesh.divx;
        divy = s_mesh.divy;
        u = s_sol.u ;
        dofel = s_mesh.dofel;
        nnodel = s_mesh.nnodel ;
        ndofel = s_mesh.ndofel ;
        nnode = s_mesh.nnode ;
        elstat = s_mesh.elstat;
        nfloatnodel = s_mesh.nfloatnodel;
        dofs = s_mesh.dofs;
        fbnd = s_param.fbnd;
        adof = s_mesh.adof;
        sig = s_sol.sig ;
        eps = s_sol.eps ;
        dofnode = s_mesh.dofnode ;
        ndofnode = s_mesh.ndofnode ;
        vn = s_param.vn ;
        elsin = s_mesh.elsin ;
        elsbnd = s_mesh.elsbnd ;
        fadof = s_mesh.fadof ;
        p = s_param.p ;
        elsoutdom = s_mesh.elsoutdom;
        
        
    % FUNCTION CALLS
        if n == 1 % mesh plotting
            pltmesh( nel, con, x, divx, divy, Lx, Ly ) 
        elseif n == 2 % deformed plotting
            pltdof( nel, con, x, u, dofel, ndofel, nnodel, Lx, Ly, divx, divy, nnode, elstat, dofnode, ndofnode, elsin, elsbnd, fadof, p )
        elseif n == 3 % level set plotting
            pltlevelset( u(dofs,1), x(dofs,:), Lx, Ly, divx, divy )
        elseif n == 4 % plotting conforming mesh
            pltmeshfnm( nel, con, x, divx, divy, Lx, Ly, elstat, nnodel, nfloatnodel, fbnd, u, elsoutdom )
        elseif n == 5 % plot density
            pltdensity( nel, con, x, u , divx, divy, Lx, Ly, elstat, nnodel, nfloatnodel, adof )
        elseif n == 6 % stress and strain fields
            pltstr( nel, con, x, sig, eps, nnodel, nfloatnodel, Lx, Ly, divx, divy, elstat, p )
        elseif n == 7 % velocity field
            pltvel( nel, con, x, vn, Lx, Ly, divx, divy, elstat )
        end
end


function [] = pltmesh( nel, con, x , divx, divy, Lx, Ly)
% plots the mesh in the undeformed configuration

    figure(101) 
    clf
    hold on
    for iel = 1:nel
        nds = con(iel,:);
        xcoord = x(nds,1);
        ycoord = x(nds,2);
        patch(xcoord,ycoord,'blue','FaceAlpha',.2) ;
%         plot(xcoord,ycoord,'ko') ;
    end
    hold off
    if Lx > Ly
        s = 800*[1 Ly/Lx];
    else
        s = 800*[Lx/Ly 1];
    end
    set(gcf, 'Position', [100 100 s]) 
    box off
    dx = Lx/divx;
    dy = Ly/divy;
    pbaspect([(Lx + 2*dx)/(Ly + 2*dy) 1 1]) ;
    xlim([-dx Lx+dx])
    ylim([-dy Ly+dy])

end


function [] = pltdof( nel, con, x, u, dofel, ndofel, nnodel, Lx, Ly, divx, divy, nnode, elstat, dofnode, ndofnode, elsin, elsbnd, fadof, p )
% plots the mesh in the deformed configuration

    figure(201)
    clf
    dx = Lx/divx;
    dy = Ly/divy;
    delta = max( max( abs( u(:,1) ) ) );
    f = 0.8 * ( min([dx dy]) / delta ) ;
    hold on
    for iel = elsin'
        scon = [1 2 3 4] ;
        nnodesub = size(scon,2);
        ndofsub = ndofnode * nnodesub ;
        nds = con(iel,scon) ;
        dfs = zeros(ndofsub,1) ;
        for ind = 1:nnodesub
            for idof = 1:ndofnode
                dfs(ndofnode * ind + idof - ndofnode,1) = ndofnode * nds(ind) + idof - ndofnode ;
            end
        end
        xcoord = x(nds,1);
        ycoord = x(nds,2);
        dof = u( dfs, 1);
        uel = zeros(nnodesub,1);
        for ind = 1:nnodesub
            uel(ind) = dof( 2*ind - 1, 1);
        end
        patch(xcoord,ycoord, uel) ;
    end
    for iel = elsbnd'
        nsub = fnm_partitioninfo( elstat(iel) ) ;
        for isub = 1:nsub
            scon = fnm_subcon( elstat(iel), isub ) ;
            nnodesub = size(scon,2);
            ndofsub = ndofnode * nnodesub ;
            nds = con(iel,scon) ;
            dfs = zeros(ndofsub,1) ;
            for ind = 1:nnodesub
                for idof = 1:ndofnode
                    dfs(ndofnode * ind + idof - ndofnode,1) = ndofnode * nds(ind) + idof - ndofnode ;
                end
            end
            xcoord = x(nds,1);
            ycoord = x(nds,2);
            dof = u( dfs, 1);
            pel = p( nds, 1);
            rho = ls_density(pel);
            if rho == 1
                uel = zeros(nnodesub,1);
                for ind = 1:nnodesub
                    uel(ind) = dof( 2*ind - 1, 1);
                end
                patch(xcoord,ycoord, uel) ;
            end
        end
    end
    for idf = fadof'
        ind = find(dofnode(:,1) == idf | dofnode(:,2) == idf) ;
        xcoord = x(ind,1);
        ycoord = x(ind,2);
        plot(xcoord,ycoord,'ro','MarkerSize',2,'MarkerFaceColor','red') ;
    end
    for iel = 1:nel
        scon = [1 2 3 4];
        nnodesub = size(scon,2);
        ndofsub = ndofnode * nnodesub ;
        ndsr = con(iel,scon) ;
        dfsr = zeros(ndofsub,1) ;
        for ind = 1:nnodesub
            for idof = 1:ndofnode
                dfsr(ndofnode * ind + idof - ndofnode,1) = ndofnode * ndsr(ind) + idof - ndofnode ;
            end
        end
        xcoord = x(ndsr,1);
        ycoord = x(ndsr,2);
        dof = u( dfsr, 1);
        uelr = zeros(nnodesub,1);
        velr = zeros(nnodesub,1);
        for ind = 1:nnodesub
            uelr(ind) = dof( 2*ind - 1, 1);
            velr(ind) = dof( 2*ind, 1);
        end
        dfs = dofel(iel,:);
        uel = u(dfs,1);
        if any(isnan(uel))
            patch(xcoord,ycoord,'red','FaceAlpha',0,'EdgeAlpha',1,'EdgeColor','red','LineWidth',0.75) ;
        end
        text( min(xcoord) + ( max(xcoord) - min(xcoord) )/2 , ...
              min(ycoord) + ( max(ycoord) - min(ycoord) )/2 , ...
              num2str(elstat(iel)), 'FontSize', 6)
    end
    hold off
    if Lx > Ly
        s = 1650*[1 Ly/Lx];
    else
        s = 800*[Lx/Ly 1];
    end
    set(gcf, 'Position', [10 50 s]) 
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
    f = 2.5 * ( min([dx dy]) / delta ) ;
    hold on
    for iel = [elsin;elsbnd]'
        nsub = fnm_partitioninfo( elstat(iel) ) ;
        for isub = 1:nsub
            scon = fnm_subcon( elstat(iel), isub ) ;
            nnodesub = size(scon,2);
            ndofsub = ndofnode * nnodesub ;
            nds = con(iel,scon) ;
            dfs = zeros(ndofsub,1) ;
            for ind = 1:nnodesub
                for idof = 1:ndofnode
                    dfs(ndofnode * ind + idof - ndofnode,1) = ndofnode * nds(ind) + idof - ndofnode ;
                end
            end
            xcoord = x(nds,1);
            ycoord = x(nds,2);
            dof = u( dfs, 1);
            pel = p( nds, 1);
            rho = ls_density(pel);
            if rho == 1
                uel = zeros(nnodesub,1);
                vel = zeros(nnodesub,1);
                for ind = 1:nnodesub
                    uel(ind) = dof( 2*ind - 1, 1);
                    vel(ind) = dof( 2*ind, 1);
                end
%                 patch(xcoord + f*uel,ycoord + f*vel, vel) ;
                patch(xcoord,ycoord,'white','FaceAlpha',0,'EdgeAlpha',0.5) ;
            end
        end
    end
    for iel = [elsin;elsbnd]'
        nsub = fnm_partitioninfo( elstat(iel) ) ;
        for isub = 1:nsub
            scon = fnm_subcon( elstat(iel), isub ) ;
            nnodesub = size(scon,2);
            ndofsub = ndofnode * nnodesub ;
            nds = con(iel,scon) ;
            dfs = zeros(ndofsub,1) ;
            for ind = 1:nnodesub
                for idof = 1:ndofnode
                    dfs(ndofnode * ind + idof - ndofnode,1) = ndofnode * nds(ind) + idof - ndofnode ;
                end
            end
            xcoord = x(nds,1);
            ycoord = x(nds,2);
            dof = u( dfs, 1);
            pel = p( nds, 1);
            rho = ls_density(pel);
            if rho == 1
                uel = zeros(nnodesub,1);
                vel = zeros(nnodesub,1);
                for ind = 1:nnodesub
                    uel(ind) = dof( 2*ind - 1, 1);
                    vel(ind) = dof( 2*ind, 1);
                end
                patch(xcoord + f*uel,ycoord + f*vel, vel) ;
                %             patch(xcoord,ycoord, vel) ;
            end
        end
    end
    hold off
    if Lx > Ly
        s = 800*[1 Ly/Lx];
    else
        s = 800*[Lx/Ly 1];
    end
    set(gcf, 'Position', [300 50 s]) 
    box off
    cb = colorbar('XTick', linspace(min(vv(:,1)), max(vv(:,1)), 10) );
    title('v')
    pbaspect([(Lx + 6*dx)/(Ly + 6*dy) 1 1]) ;
    xlim([-3*dx Lx+3*dx])
    ylim([-3*dy Ly+3*dy])
    


end


function [] = pltlevelset( u, x, Lx, Ly, divx, divy )
% plots the level set field

    figure(301) 
    clf
    
    phiv = u(:,1);
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
    Mcont = contourf(X,Y,Z, [0 0], 'LineStyle', 'none') ;
%     save('contour','Mcont')
    gr2 = gray(1);
    colormap(gr2)
    set(gca, 'visible','off')
    set(gca, 'xtick',[])
    if Lx > Ly
        s = 2000*[1 Ly/Lx];
    else
        s = 1300*[Lx/Ly 1];
    end
    set(gcf, 'Position', [50 50 s])
    box off
    dx = Lx/divx;
    dy = Ly/divy;
    title('Zero level set')
    xlabel('x')
    ylabel('y')
    pbaspect([(Lx + 2*dx)/(Ly + 2*dy) 1 1]) ;
    xlim([-dx Lx+dx])
    ylim([-dy Ly+dy])
%     grid on
    

end


function [] = pltmeshfnm( nel, con, x , divx, divy, Lx, Ly, elstat, nnodel, nfloatnodel, fbnd, phi, elsoutdom )
% plots the mesh in the undeformed configuration

    figure(401) 
    clf
    hold on
    els = setdiff((1:nel)',elsoutdom);
    for iel = els'
        nsub = fnm_partitioninfo( elstat(iel) ) ;
        for isub = 1:nsub
            scon = fnm_subcon( elstat(iel), isub ) ;
            nds = con(iel,scon);
            xcoord = x(nds,1);
            ycoord = x(nds,2);
            if elstat(iel) == 0
                patch(xcoord,ycoord,'black','FaceAlpha',.2) ;
            elseif elstat(iel) == 11 || elstat(iel) == 12 || elstat(iel) == 13 || elstat(iel) == 14
                patch(xcoord,ycoord,'blue','FaceAlpha',.2) ;
            elseif elstat(iel) == 21 || elstat(iel) == 22
                patch(xcoord,ycoord,'green','FaceAlpha',.2) ;
            elseif elstat(iel) == 31 || elstat(iel) == 32
                patch(xcoord,ycoord,'red','FaceAlpha',.2) ;
            elseif elstat(iel) == 41 || elstat(iel) == 42 || elstat(iel) == 43 || elstat(iel) == 44
                patch(xcoord,ycoord,'yellow','FaceAlpha',.2) ;
            elseif elstat(iel) == 51 || elstat(iel) == 52
                patch(xcoord,ycoord,'magenta','FaceAlpha',.2) ;
            else
                patch(xcoord,ycoord,'red','FaceAlpha',.2) ;
            end
            %                 plot(xcoord,ycoord,'ko') ;
        end
    end
    for ind = fbnd'
        xcoord = x(ind,1);
        ycoord = x(ind,2);
        plot(xcoord,ycoord,'ro') ;
    end
    for iel = 1:nel
        scon = [1 2 3 4];
        ndsr = con(iel,scon) ;
        phiel = phi(ndsr,1);
        xcoord = x(ndsr,1);
        ycoord = x(ndsr,2);
        if phiel(1)==phiel(2) && phiel(3)==phiel(4) && phiel(1)==phiel(3)
            patch(xcoord,ycoord,'red','FaceAlpha',0,'EdgeAlpha',1,'EdgeColor','red','LineWidth',0.75) ;
        end
    end
    hold off
    if Lx > Ly
        s = 800*[1 Ly/Lx];
    else
        s = 800*[Lx/Ly 1];
    end
    set(gcf, 'Position', [100 100 s]) 
    box off
    dx = Lx/divx;
    dy = Ly/divy;
    pbaspect([(Lx + 2*dx)/(Ly + 2*dy) 1 1]) ;
    xlim([-dx Lx+dx])
    ylim([-dy Ly+dy])

end


function [] = pltdensity( nel, con, x, u , divx, divy, Lx, Ly, elstat, nnodel, nfloatnodel, adof )
% plots the mesh in the undeformed configuration

    figure(501) 
    clf
    hold on
    for iel = 1:nel
        phiel = u(con(iel,:), 1) ;
        nsub = fnm_partitioninfo( elstat(iel) ) ;
        for isub = 1:nsub
            scon = fnm_subcon( elstat(iel), isub ) ;
            rho = ls_density( phiel(scon,1) ) ;
            if rho == 1
                nds = con(iel,scon);
                xcoord = x(nds,1);
                ycoord = x(nds,2);
                patch(xcoord,ycoord,'black','FaceAlpha',.5) ;
            end
        end
    end
%     for ind = adof
%         xcoord = x(ind,1);
%         ycoord = x(ind,2);
%         plot(xcoord,ycoord,'ro') ;
%     end
    hold off
    if Lx > Ly
        s = 800*[1 Ly/Lx];
    else
        s = 800*[Lx/Ly 1];
    end
    set(gcf, 'Position', [100 100 s]) 
    box off
    dx = Lx/divx;
    dy = Ly/divy;
    pbaspect([(Lx + 2*dx)/(Ly + 2*dy) 1 1]) ;
    xlim([-dx Lx+dx])
    ylim([-dy Ly+dy])

end


function [] = pltstr( nel, con, x, sig, eps, nnodel, nfloatnodel, Lx, Ly, divx, divy, elstat, phi )
% plots the mesh in the deformed configuration

    figure(601)
    clf
    dx = Lx/divx;
    dy = Ly/divy;
    hold on
    for iel = 1:nel      
        nsub = fnm_partitioninfo( elstat(iel) ) ;
        for isub = 1:nsub
            scon = fnm_subcon( elstat(iel), isub ) ;
            nds = con(iel,scon) ;
            xcoord = x(nds,1);
            ycoord = x(nds,2);
            rho = ls_density(phi(nds,1));
            sxx = sig(nds, 1) ;
            if rho == 1
                patch(xcoord,ycoord, sxx, 'EdgeAlpha', 0) ;
            end
        end
    end
    hold off
    if Lx > Ly
        s = 800*[1 Ly/Lx];
    else
        s = 800*[Lx/Ly 1];
    end
    set(gcf, 'Position', [50 50 s]) 
    box off
    colormap(parula(16))
    cb = colorbar('XTick', linspace(min(sig(:,1)), max(sig(:,1)), 10) );
    title('sxx')
    pbaspect([(Lx + 2*dx)/(Ly + 2*dy) 1 1]) ;
    xlim([-dx Lx+dx])
    ylim([-dy Ly+dy])
    
    
    figure(602)
    clf
    dx = Lx/divx;
    dy = Ly/divy;
    hold on
    for iel = 1:nel
        nsub = fnm_partitioninfo( elstat(iel) ) ;
        for isub = 1:nsub
            scon = fnm_subcon( elstat(iel), isub ) ;
            nds = con(iel,scon) ;
            xcoord = x(nds,1);
            ycoord = x(nds,2);
            exx = eps(nds, 1) ;
            patch(xcoord,ycoord, exx) ;
        end
    end
    hold off
    if Lx > Ly
        s = 800*[1 Ly/Lx];
    else
        s = 800*[Lx/Ly 1];
    end
    set(gcf, 'Position', [500 50 s]) 
    box off
    cb = colorbar('XTick', linspace(min(eps(:,1)), max(eps(:,1)), 11) );
    title('exx')
    pbaspect([(Lx + 2*dx)/(Ly + 2*dy) 1 1]) ;
    xlim([-dx Lx+dx])
    ylim([-dy Ly+dy])


end


function [] = pltvel( nel, con, x, vn, Lx, Ly, divx, divy, elstat )

    figure(701)
    clf
    dx = Lx/divx;
    dy = Ly/divy;
    hold on
    for iel = 1:nel
        scon = [1 2 3 4] ;
        nds = con(iel,scon) ;
        xcoord = x(nds,1);
        ycoord = x(nds,2);
        vel = vn(nds,1);
        patch(xcoord, ycoord, vel, 'FaceAlpha', 1, 'EdgeAlpha',0) ;
    end
    hold off
    if Lx > Ly
        s = 800*[1 Ly/Lx];
    else
        s = 800*[Lx/Ly 1];
    end
    set(gcf, 'Position', [1600 50 s]) 
    box off
    colormap(gray)
    cb = colorbar('XTick', linspace(min(vn(:,1)), max(vn(:,1)), 11) );
    title('velocity')
    pbaspect([(Lx + 2*dx)/(Ly + 2*dy) 1 1]) ;
    xlim([-dx Lx+dx])
    ylim([-dy Ly+dy])


end