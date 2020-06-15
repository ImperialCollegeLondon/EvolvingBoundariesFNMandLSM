%% TOPOPT
%
%   DATA FILE

function [s_mesh, s_geo, s_param, s_sol, s_hist] = topt_datafile(lenx, leny, dx, dy)

    % input data
    F = -1; %[N]
    E = 1; %[MPa]
    nu = 0.3;
    Lx = lenx; %[mm]
    Ly = leny; %[mm]
    thck = 1; %[mm]

    % mesh parameters
    divx = dx;
    divy = dy;
    nel = divx * divy;
    Ael = (Lx * Ly) / nel;
    nnodel = 4;
    ndofnode = 2; % u, v
    nnode = (divx + 1) * (divy + 1);
    ndofel = ndofnode * nnodel;
    ndof = ndofnode * nnode;
    l_nnodel = 4;
    l_ndofnode = 1; % phi
    l_ndofel = l_ndofnode * l_nnodel;
    l_ndof = l_ndofnode * nnode;
    ndim = 2;
    ngpt = 4;
    nstr = 3; % number of stress/strain components
    h = min([ Lx/divx, Ly/divy]) ;
      
    % initialization of vectors and matrices
    x = zeros(nnode, ndim);
    con = zeros(nel, nnodel);
    dofel = zeros(nel, ndofel);
    l_dofel = zeros(nel, l_ndofel);
    dofnode = zeros(nnode, ndofnode);
    l_dofnode = zeros(nnode, l_ndofnode);
    xgp = zeros(ngpt, ndim);
    wgp = zeros(ngpt, ndim);
    dofs = (1:ndof)';
    adof = dofs;
    
    u = sparse(ndof, 1);
    k = sparse(ndof, ndof);
    fe = sparse(ndof, 1);
    l_u = sparse(l_ndof, 1);
    l_du = sparse(l_ndof, 1);
    l_k = sparse(l_ndof, l_ndof);
    l_fe = sparse(l_ndof, 1);
    l_fi = sparse(l_ndof, 1);

    vn = zeros(nnode, 1);
    dt = 0;
    dtau = h/2;
    cond = 0;
    vol = 0;
    obj = 0;
    lamb = 0;
    lmult = 1;
    rpen = 1;
    
    
    rho = zeros(nel, 1);
    p = zeros(nnode, 1);
    rinit = 0;
    c = 0.1;
    G = [];
    
    objhist = [];
    volhist = [];
    lhist = [];
 
    % boundary conditions
    ubc = 1; % cantilever load middle
%     ubc = 2; % half bridge
%     ubc = 3; % complete bridge
%     ubc = 4; % cantilever load bottom

    
    % definition of the structure variables
    s_mesh.divx = divx;
    s_mesh.divy = divy;
    s_mesh.nel = nel;
    s_mesh.nnodel = nnodel;
    s_mesh.ndofnode = ndofnode;
    s_mesh.nnode = nnode;
    s_mesh.ndofel = ndofel;
    s_mesh.ndof = ndof;
    s_mesh.ndim = ndim;
    s_mesh.ngpt = ngpt;
    s_mesh.x = x;
    s_mesh.con = con;
    s_mesh.dofel = dofel;
    s_mesh.dofnode = dofnode;
    s_mesh.xgp = xgp;
    s_mesh.wgp = wgp;
    s_mesh.dofs = dofs;
    s_mesh.adof = adof;
    s_mesh.nstr = nstr ;
    s_mesh.l_nnodel = l_nnodel;
    s_mesh.l_ndofnode = l_ndofnode;
    s_mesh.l_ndofel = l_ndofel;
    s_mesh.l_ndof = l_ndof;
    s_mesh.l_dofel = l_dofel;
    s_mesh.l_dofnode = l_dofnode;
    
    s_geo.E = E;
    s_geo.nu = nu;
    s_geo.Lx = Lx;
    s_geo.Ly = Ly;
    s_geo.thck = thck;
    
    s_param.F = F;
    s_param.ubc = ubc;
    s_param.rho = rho;
    s_param.p = p;
    s_param.rinit = rinit;
    s_param.c = c ;
    s_param.h = h ;
    s_param.G = G ;
    s_param.dtau = dtau;
    s_param.dt = dt;
    s_param.cond = cond;
    s_param.vol = vol;
    s_param.lamb = lamb;
    s_param.lmult = lmult;
    s_param.rpen = rpen;
    s_param.Ael = Ael;
    s_param.elU = [];
    s_param.ndU = [];
    
    s_sol.u = u;
    s_sol.k = k;
    s_sol.fe = fe;
    s_sol.vn = vn;
    s_sol.obj = obj;
    s_sol.l_u = l_u;
    s_sol.l_du = l_du;
    s_sol.l_k = l_k;
    s_sol.l_fe = l_fe;
    s_sol.l_fi = l_fi;

    s_hist.objhist = objhist;
    s_hist.volhist = volhist;
    s_hist.lhist = lhist;

end

