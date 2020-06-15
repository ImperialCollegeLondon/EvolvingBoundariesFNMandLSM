% **************************************************************************
%  MODELLING EVOLVING BOUNDARIES USING FNM AND LSM (April 2020)
% **************************************************************************
%
%
%   runscript.m
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


%% ------
% OPTIONS
% -------

    % RUNTIME PLOTTING OPTIONS
    % 0 = no plot  |  1 = plot

        s_plot.mesh = 0; % initial mesh
        s_plot.ls = 1; % ls field
        s_plot.fnm = 0; % mesh partition
        s_plot.dens = 1; % meshed domain
        s_plot.dof = 0; % displacement solution
        s_plot.str = 0; % stress and strain fields
        s_plot.vel = 0; % velocity field


    % RUNTIME OUTPUT OPTIONS
    % 0 = no output  |  1 = output

        s_plot.save = 1; % save .mat file for each iteration
        
        
        
%% ----------
% RUN COMMAND
% -----------
%
% main(testcase, Lx, Ly, dx, dy, nhx, nhy, volfrac, alpha, maxiter, s_plot)


    % CANTILEVER BEAM EXAMPLE
        [obj, vol, file] = main(1, 2, 1, 100, 50, 7, 5, 0.5, 1, 10, s_plot);
        
    % MBB BEAM EXAMPLE
        %[obj, vol, file] = main(2, 3, 2, 60, 40, 7, 3, 0.5, 1, 500, s_plot);
        
    % L-BRACKET EXAMPLE
        %[obj, vol, file] = main(3, 2, 2, 80, 80, 9, 9, 0.25, 1, 500, s_plot);
        
    % BRACKET WITH 2 HOLES
        %[obj, vol, file] = main(4, 3, 1, 300, 100, 17, 5, 0.15, 1, 500, s_plot);

        
    % SAVING OUTPUT VARIABLES
        save(sprintf('data/%s_output',file),'obj', 'vol', '-v7.3')