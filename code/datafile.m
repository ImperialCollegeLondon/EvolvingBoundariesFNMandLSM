% **************************************************************************
%  MODELLING EVOLVING BOUNDARIES USING FNM AND LSM (April 2020)
% **************************************************************************
%
%
%   datafile.m
%
%
%   FUNCTIONS:
%       - datafile - creates the structure array variables with data 
%                    needed for the numerical method (including the  
%                    input data), and initializes vectors and matrices
%
%
% **************************************************************************
%                   developed by Rui O.S.S. da Costa and Silvestre T. Pinho
%                                      Imperial College London, Aeronautics
% **************************************************************************


function [s_mesh, s_geo, s_param, s_sol] = datafile(analysis, lenx, leny, dx, dy, testcase, volfrac)

    % **************************************************************************
    %
    %   INPUTS:
    %       analysis - string used to toggle between mechanical FEM and LSFEM 
    %                  analyis types
    %       lenx - length of the domain in the horizontal direction
    %       leny - length of the domain in the vertical direction
    %       dx - number of elements along the horizontal direction
    %       dy - number of elements along the vertical direction
    %       testcase - type of problem to be analysed
    %       volfrac - volume fraction
    %
    %   OUTPUTS:
    %       s_mesh - structure array used for mesh data
    %       s_geo - structure array used for geometrical data
    %       s_param - structure array used for generic analysis parameters
    %       s_sol - structure array for the the solution data
    %
    % **************************************************************************


    % ASSIGING INPUT AND INITIAL DATA
        F = -1; %[N]
        E = 1; %[MPa]
        nu = 0.3;
        Lx = lenx; %[mm]
        Ly = leny; %[mm]
        divx = dx;
        divy = dy;
        thck = 1; %[mm]

        
    % CONSTITUTIVE MATRIX
        pstress = 1;
        if pstress
            c = E/(1 - nu^2);
            D = c*[ 1, nu,          0 ;
                nu,  1,          0 ;
                0,  0, 0.5*(1-nu)];
        end
        D = sparse(D);
        Ael = (Lx/divx) * (Ly/divy) ;

        
    % MESH PARAMETERS
        if strcmp(analysis, 'fem')
            nel = divx * divy;
            nnodel = 4;
            ndofnode = 2;
            ngptquad = 4;
            ngpttri = 3;
        elseif strcmp(analysis, 'ls')
            nel = divx * divy;
            nnodel = 4;
            ndofnode = 1;
            ngptquad = 4;
            ngpttri = 3;
        end
        nnode = (divx + 1) * (divy + 1);
        ndofel = ndofnode * nnodel;
        ndof = ndofnode * nnode;
        ndim = 2;
        nstr = 3; 
        h = min([ Lx/divx, Ly/divy]) ;
        dofs = (1:ndof)';
        adof = dofs;
    
    
    % FNM MESH PARAMETERS
        if ~strcmp(analysis,'fem')
            nfloatnodel = 5;
            nfloatdofnode = 1;
            nfloatnode = divx*(divy + 1) + (divx + 1)*divy + nel;
        else
            nfloatnodel = 5;
            nfloatdofnode = 2;
            nfloatnode = divx*(divy + 1) + (divx + 1)*divy + nel;
        end
        nfloatdofel = nfloatdofnode * nfloatnodel ;
        nfloatdof = nfloatdofnode * nfloatnode ;
        floatdofs = (ndof + 1):(ndof + nfloatdof);
        if nfloatdof > 0
            ndof = ndof + nfloatdof ;
            nnode = nnode + nfloatnode ;
            nnodel = nnodel + nfloatnodel ;
            ndofel = ndofel + nfloatdofel ;
        end
        elstat = zeros(nel, 1);
        
        
    % LEVEL SET AND OPTIMISATION VARIABLES
        vn = zeros(nnode, 1);
        dt = 0;
        dtau = h/2;
        rho = zeros(nel,1);
        p = zeros(nnode - nfloatnode,1);
        rinit = 0;
        cond = 0;
        c = 0.1;
        fbnd = zeros(nnode,1);
        lmult = 1 ;
        pen = 1 ;
        vfrac = volfrac;
        obj = 0;
        vol = 0;
 
        
    % TESTCASE AND BOUNDARY CONDITIONS
        ubc = testcase; % 1 - cantilever beam;
                        % 2 - L-bracket;
                        % 3 - MBB beam;
                        % 4 - bracket with holes;
    
     
    % INITIALIZATION OF VECTORS AND MATRICES
        x = zeros(nnode, ndim);
        con = zeros(nel, nnodel);
        dofel = zeros(nel, ndofel);
        dofnode = zeros(nnode, ndofnode);
        xgpquad = zeros(ngptquad, ndim);
        wgpquad = zeros(ngptquad, ndim);
        xgptri = zeros(ngpttri, ndim);
        wgptri = zeros(ngpttri, ndim);

        u = sparse(ndof, 1);
        k = sparse(ndof, ndof);
        fe = sparse(ndof, 1);
        sig = zeros(nnode, nstr);
        eps = zeros(nnode, nstr);
    
    
    

    
    % ASSIGINING THE OUTPUT STRUCTURE ARRAY VARIABLES
        s_mesh.divx = divx; % number of elements along the horizontal direction
        s_mesh.divy = divy; % number of elements along the vertical direction
        s_mesh.nel = nel; % total number of elements
        s_mesh.nnodel = nnodel; % number of total nodes per element (real and floating)
        s_mesh.ndofnode = ndofnode; % number of dofs per node
        s_mesh.nnode = nnode; % total number of nodes (real and floating)
        s_mesh.ndofel = ndofel; % number of total dof per element (real and floating)
        s_mesh.ndof = ndof; % total number of dof (real and floating)
        s_mesh.ndim = ndim; % number of spatial dimensions (1 = 1D, 2 = 2D, 3 = 3D)
        s_mesh.ngptquad = ngptquad; % number of integration points for the quadrilateral elements
        s_mesh.ngpttri = ngpttri; % number of integration points for the triangular elements
        s_mesh.x = x; % matrix of nodal coordinates
        s_mesh.con = con; % matrix of global elemental connectivities
        s_mesh.dofel = dofel; % matrix of global element-dof connectivities
        s_mesh.dofnode = dofnode; % matrix of global node-dof connectivities
        s_mesh.xgpquad = xgpquad; % matrix of coordinates for the integration points for quadrilateral elements
        s_mesh.wgpquad = wgpquad; % matrix of weights for the integration points for quadrilateral elements
        s_mesh.xgptri = xgptri; % matrix of coordinates for the integration points for triangular elements
        s_mesh.wgptri = wgptri; % matrix of weights for the integration points for triangular elements
        s_mesh.dofs = dofs; % vector that lists global indices of all DoFs
        s_mesh.adof = adof; % vector that lists global indicies of all active DoFs
        s_mesh.fadof = []; % vector that lists global indicies of all active floating DoFs
        s_mesh.nstr = nstr ; % number os stress/strain variables
        s_mesh.nfloatnode = nfloatnode ; % total number of floating nodes
        s_mesh.nfloatnodel = nfloatnodel; % number floating nodes per element
        s_mesh.nfloatdofnode = nfloatdofnode; % number of flating dof per floating node
        s_mesh.nfloatdofel = nfloatdofel; % number of floating dof per element  
        s_mesh.nfloatdof = nfloatdof; % total number of floating dof
        s_mesh.floatdofs = floatdofs; % vector that lists all floating dofs
        s_mesh.elstat = elstat; % vector containing the flag defining the element partitioning status
        s_mesh.elsin = []; % vector that lists all the elements inside domain \Omega at the current iteration
        s_mesh.elsout = []; % vector that lists all the elements outside domain \Omega at the current iteration
        s_mesh.elsbnd = []; % vector that lists all the elements close to the boundary that have been partitioned at the current iteration
        s_mesh.elsoutdom = []; % vector that lists all the elements permanently out of the domain due to geometric restrictions
        s_mesh.ndsoutdom = []; % vector that lists the global indices of all the nodes out of the domain due to geometric restrictions
        s_mesh.elsfixbnd = []; % vector that lists all the elements close to a fixed boundary resulting from geometric restrictions
        s_mesh.ndsfixbnd = []; % vector that lists the global indicies of all the nodes close to a fixed boundary resulting from geometric restrictions
 
        s_geo.E = E; % Young's modulus
        s_geo.nu = nu; % Poisson ratio
        s_geo.Lx = Lx; % length of the domain in the horizontal direction
        s_geo.Ly = Ly; % length of the domain in the vertical direction
        s_geo.thck = thck; % thickness
        s_geo.D = D; % elastic constitutive matrix
        s_geo.Ael = Ael ; % Area of the quadrilateral element
        s_geo.bndinfo = []; % vector that lists some parameters regarding the fixed boundaries due to the geometric restrictions
  
        s_param.F = F; % force value
        s_param.ubc = ubc; % boundary condition selector
        s_param.h = h ; % minimum element size
        s_param.analysis = analysis; % type of analysis (mechanical FEM, LSFEM)
        s_param.vn = vn; % vector of nodal normal velocities from the optimisation sensitivities
        s_param.dt = dt; % time step
        s_param.dtau = dtau; % pseudo time step used for velocity extension and reinitialization
        s_param.rho = rho; % pseudo density (0.001 for void, 1 for solid)
        s_param.p = p; % vector for the pseudo density field
        s_param.rinit = rinit; % variable that toggles the reinitialization
        s_param.cond = cond; % measure of quality of the LS field
        s_param.c = c; % paratemer used to stabilise the LSFEM solution
        s_param.fbnd = fbnd; % list of global indicies of floating nodes on the evolving boundary
        s_param.lmult = lmult; % lagrange multiplier 
        s_param.pen = pen; % penalty parameter of the constraint
        s_param.vfrac = vfrac; % volume fraction

        s_sol.u = u; % global solution vector
        s_sol.k = k; % global stiffness matrix
        s_sol.fe = fe; % global force vector
        s_sol.sig = sig; % stress matrix
        s_sol.eps = eps; % strain matrix
        s_sol.uold = u; % global solution vector of the previous iteration
        s_sol.obj = obj; % value of the objective function at the current iteration
        s_sol.vol = vol; % value of the volume of the domain at the current iteration

end

