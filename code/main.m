% **************************************************************************
%  MODELLING EVOLVING BOUNDARIES USING FNM AND LSM (April 2020)
% **************************************************************************
%
%
%   main.m
%
%
%   FUNCTIONS:
%       - main - handles calls the secondary functions and procedures
%                and controls the iterations, plotting and output actions
%       - syncvars - handles the syncing of l_* (LSFEM part) and 
%                    f_* (FEM part) structure variables
%       - activdof - determines the active DOF based on the testcase
%
%
% **************************************************************************
%                   developed by Rui O.S.S. da Costa and Silvestre T. Pinho
%                                      Imperial College London, Aeronautics
% **************************************************************************
%
%     Manuscript available at: https://doi.org/10.1016/j.cma.2020.113077
%
% **************************************************************************

function [hobj, hvol, file] = main(s_testcase, s_initLx, s_initLy, s_initdx, s_initdy, s_initnhx, s_initnhy, s_initvfrac, s_initalpha, s_maxiter, s_plot)

    % **************************************************************************
    %
    %   INPUTS:
    %       s_testcase - problem choice (cantilever, MBB, ...)
    %       s_initLx - length in the horizontal axis
    %       s_initLy - length in the vertical axis
    %       s_initdx - number of elements in the horizontal direction
    %       s_initdy - number of elements in the vertical direction
    %       s_initnhx - number of initial holes in the horizontal direction
    %       s_initnhy - number of initial holes in the vertical direction
    %       s_initvfrac - volume fraction desired for the optimisation
    %       s_initalpha - time stepping parameter
    %       s_maxiter - number of maximum iterations
    %
    %
    %   OUTPUTS:
    %       hobj - vector of the objective function at all iterations
    %       hvol - vector of the volume of the domain at all iterations
    %       file - string containing a timestamp used to save .mat files
    %
    % **************************************************************************

    
    % STARTING
        clearvars -except s_*; close all; clc
        if s_testcase==1,tcas='Cantilever beam'; elseif s_testcase==2,tcas='MBB beam'; elseif s_testcase==3,tcas='L-bracket'; elseif s_testcase==4,tcas='Bracket with 2 holes'; end
        fprintf('** STARTED %s ** \n\nCase           : %s \nDomain         : %dx%d [mm] \nMesh           : %dx%d [elements] \nInitialization : %dx%d [holes] \nTime step      : %2.2f [proportion] \nMax iterations : %d \n\n',datestr(now,'dd/mm/yyyy HH:MM:SS'), tcas, s_initLx, s_initLy, s_initdx, s_initdy, s_initnhx, s_initnhy, s_initalpha, s_maxiter)
        
        tic  

        
        
    % INITIALIZATION OF DATA & VARIABLES
        [f_mesh, f_geo, f_param, f_sol] = datafile('fem', s_initLx, s_initLy, s_initdx, s_initdy, s_testcase, s_initvfrac) ;
        [l_mesh, l_geo, l_param, l_sol] = datafile('ls', s_initLx, s_initLy, s_initdx, s_initdy, s_testcase, s_initvfrac) ;
        file = strcat(datestr(now,'yyyymmdd_HHMMSS'),'_',tcas,'_','LxLy_',num2str(s_initLx),'x',num2str(s_initLy),'_dxdy_',num2str(s_initdx),'x',num2str(s_initdy),'_nhxnhy_',num2str(s_initnhx),'x',num2str(s_initnhy),'_ts_',num2str(s_initalpha),'_maxit_',num2str(s_maxiter));
        hobj = zeros(s_maxiter,1); hvol = zeros(s_maxiter,1);
        

        
    % CREATION OF THE MESH
        [f_mesh.x, f_mesh.con, f_mesh.dofel, f_mesh.dofnode, f_mesh.xgpquad, f_mesh.wgpquad, f_mesh.xgptri, f_mesh.wgptri] = fem_mesh(f_mesh, f_geo, f_param) ;
        [l_mesh.x, l_mesh.con, l_mesh.dofel, l_mesh.dofnode, l_mesh.xgpquad, l_mesh.wgpquad, l_mesh.xgptri, l_mesh.wgptri] = ls_mesh(l_mesh, l_geo, l_param) ;

        if s_plot.mesh, plotfunc(1, f_mesh, f_geo, f_param, f_sol); end

        % modifying the connectivities to include the floating nodes
        [f_mesh.con, f_mesh.dofel, f_mesh.dofnode] = fnm_confnm( f_mesh, f_geo, f_param, f_sol ) ;
        [l_mesh.con, l_mesh.dofel, l_mesh.dofnode] = fnm_confnm( l_mesh, l_geo, l_param, l_sol ) ;

        
        
    % LEVEL SET INITIALIZATION
        [l_param.p, l_sol.u] = ls_init( l_mesh, l_geo, l_param, l_sol, s_initnhx, s_initnhy) ;
        
        % enforcing domain restrictions
        [l_sol.u, l_mesh.ndsoutdom, l_mesh.elsoutdom, l_mesh.elsfixbnd, l_geo.bndinfo] = ls_domainrestrict( l_mesh, l_geo, l_param, l_sol );
        l_mesh.adof = setdiff(l_mesh.adof, l_mesh.ndsoutdom);


        
        
    % FIRST ELEMENT PARTITIONING
        [l_param.fbnd, l_mesh.elstat, l_mesh.x, l_mesh.fadof, l_sol.u, l_mesh.elsin, l_mesh.elsout, l_mesh.elsbnd] = fnm_part( l_mesh, l_geo, l_param, l_sol ) ;
        [l_sol.u] = ls_reinit( l_mesh, l_geo, l_param, l_sol ) ;
        [f_mesh, f_param, f_geo] = syncvars(f_mesh,f_param,f_geo,l_mesh,l_param,l_sol,l_geo) ;


        
    % TIME MARCHING
        drawnow
        iter = 1;
        while iter <= s_maxiter
            
            % ACTIVE DOF
                f_mesh.fadof = activdof(f_mesh, f_geo, f_param);
                
                

            % MECHANICAL PART
                % force vector
                    [f_sol.fe] = fem_f(f_mesh, f_geo, f_param, f_sol) ;

                % stiffness matrix
                    [f_sol.k] = fem_k(f_mesh, f_geo, f_param, f_sol) ;

                % calculation of displacements
                    f_sol.u = sparse(f_mesh.ndof,1);
                    f_sol.u(f_mesh.fadof,1) = f_sol.k(f_mesh.fadof,f_mesh.fadof) \ f_sol.fe(f_mesh.fadof,1) ;

                % stress and strain computaiton
                    [f_sol.sig, f_sol.eps] = fem_post( f_mesh, f_geo, f_param, f_sol) ;

                % optimisation velocity field and timestep
                    [f_param.vn, f_sol.obj, f_sol.vol, f_param.lmult, f_param.pen] = ls_velocity( f_mesh, f_geo, f_param, f_sol ) ;
                    f_param.dt = s_initalpha * ( f_param.h / max(abs(f_param.vn)) ) ;
                    l_param.dt = f_param.dt; l_param.vn = f_param.vn;

                    
                    
                    
            % PLOTTING
                if s_plot.ls, plotfunc(3, l_mesh, l_geo, l_param, l_sol); end
                if s_plot.fnm, plotfunc(4, l_mesh, l_geo, l_param, l_sol); end
                if s_plot.dens, plotfunc(5, l_mesh, l_geo, l_param, l_sol); end
                if s_plot.dof, plotfunc(2, f_mesh, f_geo, f_param, f_sol); end
                if s_plot.str, plotfunc(6, f_mesh, f_geo, f_param, f_sol); end
                if s_plot.vel, plotfunc(7, f_mesh, f_geo, f_param, f_sol); end
                drawnow

            % PRINTING CURRENT OPTIMISATION STATE
                fprintf('%4d : obj = %15f | vol = %15f \n', iter, f_sol.obj, f_sol.vol)

            % SAVING ITERATION DATA
                if s_plot.save, save(sprintf('runs/%s_iter%d',file,iter),'l_mesh', 'l_geo','l_param','l_sol', 'f_mesh', 'f_geo', 'f_param', 'f_sol', '-v7.3'); end
                hobj(iter,1) = f_sol.obj; hvol(iter,1) = f_sol.vol;
                
                    
            % LEVEL SET PART
                % force vector
                    [l_sol.fe] = ls_f(l_mesh, l_geo, l_param, l_sol) ;

                % stiffness matrix
                    [l_sol.k] = ls_k(l_mesh, l_geo, l_param, l_sol) ;

                % calculation of level set field
                    l_sol.u(l_mesh.adof,1) = l_sol.k(l_mesh.adof,l_mesh.adof) \ l_sol.fe(l_mesh.adof,1) ;
                    f_param.p = l_sol.u;
                    
                % enforcing domain restrictions
                    l_sol.u(l_mesh.floatdofs,1) = 0 ; l_mesh.x(l_mesh.floatdofs,:) = 0 ;
                    [l_sol.u, l_mesh.ndsoutdom, l_mesh.elsoutdom, l_mesh.elsfixbnd, l_geo.bndinfo] = ls_domainrestrict( l_mesh, l_geo, l_param, l_sol );
                    

                % partitioning
                    [l_param.fbnd, l_mesh.elstat, l_mesh.x, l_mesh.fadof, l_sol.u, l_mesh.elsin, l_mesh.elsout, l_mesh.elsbnd] = fnm_part( l_mesh, l_geo, l_param, l_sol ) ;
                    [l_sol.u] = ls_reinit( l_mesh, l_geo, l_param, l_sol ) ;
                    [f_mesh, f_param,f_geo] = syncvars(f_mesh,f_param,f_geo,l_mesh,l_param,l_sol,l_geo) ;
                    

                    
            % UPDATING ITERATOR
                iter = iter + 1 ;
        end
        

    % FINISHING
        t = toc;
        fprintf('\n** FINALISED IN %6.2f SECONDS ** \n\n', t)
        if s_plot.save, save(sprintf('runs/%s_end',file),'t', '-v7.3'), end
end

function [f_mesh,f_param,f_geo] = syncvars(f_mesh,f_param,f_geo,l_mesh,l_param,l_sol,l_geo)
    
    % **************************************************************************
    %
    %   INPUTS:
    %       f_mesh - FEM mesh variables
    %       f_param - FEM parameter variables
    %       f_geo - FEM geometry variables
    %       l_mesh - LSFEM mesh variables
    %       l_param - LSFEM parameter variables
    %       l_sol - LSFEM solution variables
    %       l_geo - LSFEM geometry variables
    %
    %
    %   OUTPUTS:
    %       f_mesh - synced FEM mesh structure variable
    %       f_param - synced FEM parameter structure variable
    %       f_geo - synced FEM geometry structure variable
    %
    % **************************************************************************

    f_mesh.elstat = l_mesh.elstat ;
    f_mesh.x = l_mesh.x ;
    f_param.fbnd = l_param.fbnd ;
    f_param.p = l_sol.u;
    f_mesh.elsin = l_mesh.elsin;
    f_mesh.elsout = l_mesh.elsout;
    f_mesh.elsbnd = l_mesh.elsbnd;
    f_param.fanode = l_mesh.fadof;
    f_mesh.ndsoutdom = l_mesh.ndsoutdom;
    f_mesh.elsoutdom = l_mesh.elsoutdom;
    f_mesh.elsfixbnd = l_mesh.elsfixbnd;
    f_geo.bndinfo = l_geo.bndinfo;

end

function [fadof] = activdof(f_mesh, f_geo, f_param)

    % **************************************************************************
    %
    %   INPUTS:
    %       f_mesh - structure arry with mesh data
    %       f_geo - structure arry with geometric data
    %       f_param - structure arry with generic parameters
    %
    %   OUTPUTS:
    %       fadof - array with active DOF indices
    %
    % **************************************************************************

    % VARIABLE INITIALISATION
        ubc = f_param.ubc;
        x = f_mesh.x;
        dofnode = f_mesh.dofnode;
        p = f_param.p;
        Lx = f_geo.Lx;
        Ly = f_geo.Ly;
        h = f_param.h;
        nnode = f_mesh.nnode;
        nfloatnode = f_mesh.nfloatnode;
        ndsoutdom = f_mesh.ndsoutdom;
        fanode = f_param.fanode;
        bndinfo = f_geo.bndinfo;
        anode = setdiff(unique([(1:(nnode-nfloatnode))'; fanode]), ndsoutdom);
        fadof = unique(reshape(dofnode(anode,:), [size(dofnode(anode,:),1)*size(dofnode(anode,:),2), 1])) ;
        
        
    % CANTILEVER BEAM
        if ubc == 1 
            
            nodeindx = anode(find( x(anode,1) == 0 ), 1) ;
            dofindaux = dofnode(nodeindx, :) ;
            dofindaux = reshape(dofindaux, [size(dofindaux,1)*size(dofindaux,2), 1]) ;
            dofindaux = unique(dofindaux) ;

            nodneg = anode(find(p(anode,1) < 0) ,1) ;
            dofneg = dofnode(nodneg, :) ;
            dofneg = reshape(dofneg, [size(dofneg,1)*size(dofneg,2), 1]) ;
            dofneg = unique(dofneg) ;

            dofx = unique( [dofindaux; dofneg] ) ;

            
    % MBB BEAM
        elseif ubc == 2 

            nodeindx = anode(find( x(anode,1) == 0 ), 1);
            dofindaux = dofnode(nodeindx, 1) ;
            dofindaux = reshape(dofindaux, [size(dofindaux,1)*size(dofindaux,2), 1]);
            dofindx1 = unique(dofindaux);

            nodeindx = anode(find( x(anode,1) == Lx & x(anode,2) == 0 ), 1);
            dofindaux = dofnode(nodeindx, 2) ;
            dofindaux = reshape(dofindaux, [size(dofindaux,1)*size(dofindaux,2), 1]);
            dofindx2 = unique(dofindaux);

            dofindaux = unique([dofindx1; dofindx2]);

            nodneg = anode(find(p(anode,1) < 0),1) ;
            dofneg = dofnode(nodneg, :) ;
            dofneg = reshape(dofneg, [size(dofneg,1)*size(dofneg,2), 1]) ;
            dofneg = unique(dofneg) ;

            dofx = unique( [dofindaux; dofneg] ) ;

            
    % L BRACKET
        elseif ubc == 3 

            nodeindx = anode(find( x(anode,2) == Ly & x(anode,1) >= 0 & x(anode,1) <= 2*Lx/5+h), 1);
            dofindaux = dofnode(nodeindx, :) ;
            dofindaux = reshape(dofindaux, [size(dofindaux,1)*size(dofindaux,2), 1]);
            dofindaux = unique(dofindaux);

            nodneg = anode(find(p(anode,1) < 0), 1) ;
            dofneg = dofnode(nodneg, :) ;
            dofneg = reshape(dofneg, [size(dofneg,1)*size(dofneg,2), 1]) ;
            dofneg = unique(dofneg) ;

            dofx = unique( [dofindaux; dofneg] ) ;

            
    % BRACKET WITH TWO HOLES
        elseif ubc == 4 

            d = sqrt(sum( ( x(fanode,:) - [bndinfo(1,1) Ly/2] ).^2 ,2));
            nodeindx = fanode( find( d > bndinfo(4,1) - h/20 & d < bndinfo(4,1) + h/20 ) ,1);
            dofindaux = dofnode(nodeindx, :) ;
            dofindaux = reshape(dofindaux, [size(dofindaux,1)*size(dofindaux,2), 1]);
            dofindaux = unique(dofindaux);

            nodneg = anode(find(p(anode,1) < 0),1) ;
            dofneg = dofnode(nodneg, :) ;
            dofneg = reshape(dofneg, [size(dofneg,1)*size(dofneg,2), 1]) ;
            dofneg = unique(dofneg) ;

            dofx = unique( [dofindaux; dofneg] ) ;

        end
    
    % OBTAINING ARRAY OF ACTIVE DOF
        fadof = setdiff(fadof,dofx);
end