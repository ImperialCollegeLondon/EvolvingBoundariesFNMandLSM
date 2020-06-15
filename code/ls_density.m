% **************************************************************************
%  MODELLING EVOLVING BOUNDARIES USING FNM AND LSM (April 2020)
% **************************************************************************
%
%
%   ls_density.m
%
%
%   FUNCTIONS:
%       - ls_density - calculates the \rho, the pseudo-density field used
%                      to select solid elements
%
%
% **************************************************************************
%                   developed by Rui O.S.S. da Costa and Silvestre T. Pinho
%                                      Imperial College London, Aeronautics
% **************************************************************************

function rho = ls_density( phiel )

    % **************************************************************************
    %
    %   INPUTS:
    %       phiel - vector of the elemental LS values
    %
    %   OUTPUTS:
    %       rho - pseudo-density value
    %
    % **************************************************************************

    % COMPUTE RHO BASED ON THE SIGN OF THE AVERAGE OF THE LS VALUES
        if sum(phiel,1)/size(phiel,1) < 0
            rho = 0.0001 ; % void
        else
            rho = 1 ; % solid
        end
end