%% TOPOPT 17/04/2019
%
%   2D topology optimisation using the level set integrated approach
%
%   MAIN FILE

function [] = topopt(file, initLx, initLy, initdx, initdy, initnhx, initnhy, initalpha, maxiter)

        tic

    %% INITIALIZATION OF VARIABLES

        [s_mesh, s_geo, s_param, s_sol, s_hist] = topt_datafile(initLx, initLy, initdx, initdy) ;


    %% CREATION OF THE MESH

        [s_mesh.x, s_mesh.con, s_mesh.dofel, s_mesh.l_dofel, s_mesh.dofnode, s_mesh.l_dofnode, s_mesh.xgp, s_mesh.wgp, s_mesh.adof] = topt_mesh(s_mesh, s_geo, s_param) ;
    %     topt_plot(1, s_mesh, s_geo, s_param, s_sol, s_hist) ;


    %% INITIALIZATION OF THE LEVEL SET FIELD

        [s_param.p, s_sol.l_u] = topt_initphi( s_mesh, s_geo, s_param, s_sol, initnhx, initnhy ) ;
        [s_param.rho, s_param.G] = topt_intersect( s_mesh, s_geo, s_param, s_sol ) ;
    %     topt_plot(3, s_mesh, s_geo, s_param, s_sol, s_hist) ;
    %     topt_plot(4, s_mesh, s_geo, s_param, s_sol, s_hist) ;


    %% ITERATIONS

        drawnow
        iter = 1;
        while iter <= maxiter

            % ==== MECHANICAL
            % force vector
            [s_sol.fe] = topt_fext(1, s_mesh, s_geo, s_param, s_sol) ;

            % stiffness matrix
            [s_sol.k] = topt_k(1, s_mesh, s_geo, s_param, s_sol) ;

            % calculate displacements
            s_sol.u(s_mesh.adof,1) = s_sol.k(s_mesh.adof,s_mesh.adof) \ s_sol.fe(s_mesh.adof,1) ;

            [s_sol.sig, s_sol.eps] = topt_post( s_mesh, s_geo, s_param, s_sol);
    %         topt_plot(7, s_mesh, s_geo, s_param, s_sol, s_hist);

            % velocity and time step
            [s_sol.vn, s_param.vol, s_sol.obj, s_param.lmult, s_param.rpen, s_param.elU, s_param.ndU] = topt_velocity( s_mesh, s_geo, s_param, s_sol ) ;
            s_param.dt = initalpha * (s_param.h / max(abs(s_sol.vn(:,1))) ) ;
            s_hist.objhist = [s_hist.objhist s_sol.obj];
            s_hist.volhist = [s_hist.volhist s_param.vol];

    %         topt_plot(3, s_mesh, s_geo, s_param, s_sol, s_hist)
    %         topt_plot(5, s_mesh, s_geo, s_param, s_sol, s_hist)

            % output
            fprintf('%3d: obj = %d  |  vol = %d  \n',iter, s_sol.obj, s_param.vol)

            save(sprintf('runs/%s/iter%d',file,iter),'s_mesh', 's_geo','s_param','s_sol', '-v7.3')
    %         save(sprintf('runs/iter%d',iter),'s_mesh', 's_geo','s_param','s_sol', '-v7.3')

            % ==== LEVEL SET
            % external forces
            [s_sol.l_fe] = topt_fext(2, s_mesh, s_geo, s_param, s_sol) ;

            % stiffness matrix
            [s_sol.l_k] = topt_k(2, s_mesh, s_geo, s_param, s_sol) ;

            % displacement increment
            s_sol.l_u = s_sol.l_k \ s_sol.l_fe ;

            % update intersections matrix and density field
            [s_param.rho, s_param.G] = topt_intersect( s_mesh, s_geo, s_param, s_sol ) ;
    %         topt_plot(3, s_mesh, s_geo, s_param, s_sol, s_hist)

            % calculate conditioning number
            [s_param.cond] = topt_cond( s_mesh, s_geo, s_param, s_sol );
%             if s_param.cond > 0.5
                s_param.rinit = 1;

                % reinitialization
                reiter = 1 ;
                fprintf('\n')
                while reiter <= 3

                    % external forces
                    [s_sol.l_fe] = topt_fext(2, s_mesh, s_geo, s_param, s_sol) ;

                    % stiffness matrix
                    [s_sol.l_k] = topt_k(2, s_mesh, s_geo, s_param, s_sol) ;

                    % displacement
                    s_sol.l_u = s_sol.l_k \ s_sol.l_fe ;

                    % new conditioning number
                    [s_param.cond] = topt_cond( s_mesh, s_geo, s_param, s_sol );

                    % output
                    fprintf('%d:  cond = %d  (reinit) \n',iter, s_param.cond)
                    reiter = reiter + 1;
                end
                fprintf('\n')

                s_param.rinit = 0;
%             end

            % plots
    %         topt_plot(4, s_mesh, s_geo, s_param, s_sol, s_hist) ;
    %         topt_plot(6, s_mesh, s_geo, s_param, s_sol, s_hist) ;

            % update iterator
            iter = iter + 1;
            drawnow
        end

    %     topt_plot(2, s_mesh, s_geo, s_param, s_sol, s_hist)

    t = toc;
    save(sprintf('runs/%s/end',file),'t')

end