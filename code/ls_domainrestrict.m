% **************************************************************************
%  MODELLING EVOLVING BOUNDARIES USING FNM AND LSM (April 2020)
% **************************************************************************
%
%
%   ls_domainrestrict.m
%
%
%   FUNCTIONS:
%       - ls_domainrestric - imposes fixed domain restrictions by modifying
%                            the LS field
%
%
% **************************************************************************
%                   developed by Rui O.S.S. da Costa and Silvestre T. Pinho
%                                      Imperial College London, Aeronautics
% **************************************************************************


function [u, ndsoutdom, elsoutdom, elsfixbnd, bndinfo] = ls_domainrestrict( s_mesh, s_geo, s_param, s_sol )

    % **************************************************************************
    %
    %   INPUTS:
    %       s_mesh - structure arry with mesh data
    %       s_geo - structure arry with geometric data
    %       s_param - structure arry with generic parameters
    %       s_sol - structure array with solution data
    %
    %   OUTPUTS:
    %       u - LS field
    %       ndsoutdom - vector containing the nodes outside the fixed domain to
    %                   be deactivated for analysis
    %       elsoutdom - vector containing the elements fully outside the fixed
    %                   domain, used to save computation time
    %       elsfixbnd - vector containing the elements close to the fix boundary
    %                   used in the paritioning
    %       bndinfo - vector containing information about the fixed boundary
    %
    % **************************************************************************

    % VARIABLE INITIALISATION
        u = s_sol.u ;
        nnode = s_mesh.nnode - s_mesh.nfloatnode ;
        x = s_mesh.x ;
        nel = s_mesh.nel;
        con = s_mesh.con;
        h = s_param.h;
        ubc = s_param.ubc;
        Lx = s_geo.Lx;
        Ly = s_geo.Ly;
        divx = s_mesh.divx;
        elsoutdom = [];
        ndsoutdom = [];
        elsfixbnd = [];
        bndinfo = [];

    % CHECK IF DOMAIN RESTRICTIONS ARE ALREADY DEFINED
        aux = s_geo.bndinfo;
        if isempty(aux)
            start = 1;
        else
            start = 0;
        end
    

    % CONDITIONAL CHECK FOR DIFFERENT BOUNDARY CONDITIONS
        % cantilever beam
            if ubc == 1 

                % do nothing
                
        % MBB beam
            elseif ubc == 2 

                % do nothing
                
        % L bracket - remove nodes in upper right corner square
            elseif ubc == 3 

                xs = ((0+h/2):h:Lx)';
                ys = ((0+h/2):h:Ly)';
                dx = abs(xs - (2*Lx/5 + h/100));
                xs = xs(find(dx == min(dx)),1);
                dy = abs(ys - (2*Ly/5 + h/100));
                ys = ys(find(dy == min(dy)),1);

                % impose square 0.6x0.6 in the upper right corner
                ndindx = find( x(1:nnode,1) > xs + h/100 & x(1:nnode,2) > ys + h/100) ;
                u(ndindx,1) = -1 ;

                ndsv1 = setdiff((1:nnode)', ndindx);
                for ind = ndsv1'
                    xind = x(ind,:);
                    if u(ind,1) > 0

                        if xind(1,1) < xs - h /100 && xind(1,2) < ys - h/100
                            d = sqrt(sum((xind - [xs ys]).^2));
                        elseif xind(1,1) > xs + h/100 && xind(1,2) < ys - h/100
                            d = sqrt(sum((xind - [xind(1,1) ys]).^2));
                        elseif xind(1,1) < xs - h/100 && xind(1,2) > ys + h/100
                            d = sqrt(sum((xind - [xs xind(1,2)]).^2));
                        elseif xind(1,1) > xs - h/100 && xind(1,2) > ys + h/100
                            d = min(sqrt(sum((xind - [xind(1,1) ys]).^2)), sqrt(sum((xind - [xs xind(1,2)]).^2)));
                        end

                        if u(ind,1) > d
                            u(ind,1) = d;
                        end
                    end     
                end
                for ind = ndindx'
                    xind = x(ind,:);
                    d = min(sqrt(sum((xind - [xind(1,1) ys]).^2)), sqrt(sum((xind - [xs xind(1,2)]).^2)));
                    u(ind,1) = -d;
                end

                ndsoutdom = ndindx;
                els = [];
                for ind = ndsoutdom'
                    els = [els; find(con(:,1) == ind);...
                                find(con(:,2) == ind);...
                                find(con(:,3) == ind);...
                                find(con(:,4) == ind)];
                end
                elsoutdom = unique(els);
                ndindx = find( (x(1:nnode,1) > xs - h & x(1:nnode,1) < xs + h & x(1:nnode,2) > ys - h) | ...
                               (x(1:nnode,2) > ys - h & x(1:nnode,2) < ys + h & x(1:nnode,1) > xs - h) ) ;
                els = [];
                for ind = ndindx'
                    els = [els; find(con(:,1) == ind);...
                                find(con(:,2) == ind);...
                                find(con(:,3) == ind);...
                                find(con(:,4) == ind)];
                end
                els = unique(els);
                elsoutdom = setdiff(elsoutdom, els);

        % Bracket with 2 holes - remove two circles and two half-circle rings
            elseif ubc == 4 

                xs = (0:h:Lx)';
                ys = ((0+h/2):h:Ly)';
                dx1 = abs(xs - (5*Lx/6));
                xright = xs(find(dx1 == min(dx1)),1);
                dx2 = abs(xs - (Lx/6));
                xleft = xs(find(dx2 == min(dx2)),1);
                dy = abs(ys - (4*Ly/5));
                ys = ys(find(dy == min(dy)),1);
                rads = ys(1) - Ly/2;
                radb = Ly/2 + h/100;
                radm = rads + 2*h;

                % find nodes outside big left circle
                dleft = sqrt(sum(( x(1:nnode,:) - [xleft Ly/2] ).^2,2));
                ndsbl = find( dleft > radb + h/100 & x(1:nnode,1) < Lx/6 ) ;
                % nodes inside small left circle
                ndssl = find( dleft < rads - h/100 ) ;
                % nodes outside big right circle
                dright = sqrt(sum(( x(1:nnode,:) - [xright Ly/2] ).^2,2));
                ndsbr = find( dright > radb + h/100 & x(1:nnode,1) > 5*Lx/6 ) ;
                % nodes inside small right circle
                ndssr = find( dright < rads - h/100 ) ;
                ndindx = unique([ndsbl; ndssl; ndsbr; ndssr]);


                uinit = u;
                umed = 0*u;
                usma= 0*u;
                ubig = 0*u;
                ubig2 = 0*u;
                uinit_u_med = 0*u;
                uinitrmed_i_cbig = 0*u;
                unitirmedicbig_i_sma = 0*u;
                for ind = 1:nnode
                    xind = x(ind,:);

                    dl = sqrt(sum(( xind - [xleft Ly/2] ).^2,2));
                    dr = sqrt(sum(( xind - [xright Ly/2] ).^2,2));
                    d = min(dl,dr);
                    d = d(1);

                    umed(ind,1) = - (d - radm);
                    usma(ind,1) = d - rads;
                    ubig(ind,1) = - (d - radb) ;
                    ubig2(ind,1) = (Lx/2-xleft) - abs(xind(1,1) - Lx/2);
                end
                for ind = 1:nnode
                    ubig(ind,1) = ubig(ind,1) + ubig2(ind,1) + sqrt( ubig(ind,1)^2 + ubig2(ind,1)^2 );
                    uinit_u_med(ind,1) = uinit(ind,1) + umed(ind,1) + sqrt( uinit(ind,1)^2 + umed(ind,1)^2 );
                    uinitrmed_i_cbig(ind,1) = uinit_u_med(ind,1) + ubig(ind,1) - sqrt( uinit_u_med(ind,1)^2 + ubig(ind,1)^2 );
                    unitirmedicbig_i_sma(ind,1) = uinitrmed_i_cbig(ind,1) + usma(ind,1) - sqrt( uinitrmed_i_cbig(ind,1)^2 + usma(ind,1)^2 );
                end
                u = unitirmedicbig_i_sma;


                ndsoutdom = ndindx;
                els = [];
                for ind = ndsoutdom'
                    els = [els; find(con(:,1) == ind);...
                                find(con(:,2) == ind);...
                                find(con(:,3) == ind);...
                                find(con(:,4) == ind)];
                end
                elsoutdom = unique(els);


                % find nodes around big left circle
                dleft = sqrt(sum(( x(1:nnode,:) - [xleft Ly/2] ).^2,2));
                ndsbl = find( dleft > radb - h & dleft < radb + h & x(1:nnode,1) < Lx/6 ) ;
                % nodes around small left circle
                ndssl = find( dleft > rads - h & dleft < rads + h ) ;
                % nodes around big right circle
                dright = sqrt(sum(( x(1:nnode,:) - [xright Ly/2] ).^2,2));
                ndsbr = find( dright > radb - h & dright < radb + h & x(1:nnode,1) > 5*Lx/6 ) ;
                % nodes around small left circle
                ndssr = find( dright > rads - h & dright < rads + h ) ;
                ndindx = unique([ndsbl; ndssl; ndsbr; ndssr]);

                els = [];
                for ind = ndindx'
                    els = [els; find(con(:,1) == ind);...
                                find(con(:,2) == ind);...
                                find(con(:,3) == ind);...
                                find(con(:,4) == ind)];
                end
                elsfixbnd = unique(els);
                elsoutdom = setdiff(elsoutdom, elsfixbnd);

                bndinfo = [xleft;xright;radb;rads];

            end
end