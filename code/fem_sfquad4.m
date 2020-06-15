% **************************************************************************
%  MODELLING EVOLVING BOUNDARIES USING FNM AND LSM (April 2020)
% **************************************************************************
%
%
%   fem_sfquad4.m
%
%
%   FUNCTIONS:
%       - fem_sfquad4 - calulates de shape functions and its derivatives
%                       for a quadrilateral isoparametric element with 4 nodes
%
%
% **************************************************************************
%                   developed by Rui O.S.S. da Costa and Silvestre T. Pinho
%                                      Imperial College London, Aeronautics
% **************************************************************************


function [n, dn] = fem_sfquad4( ndim, nnodel, xi, eta)

    % **************************************************************************
    %
    %   INPUTS:
    %
    %       ndim - number of spatial dimensions
    %       nnodel - number of total nodes per element
    %       xi - \xi coordinate of the natural domain
    %       eta - \eta coordinate of the natural domain
    %
    %
    %   OUTPUTS:
    %       n - vector of shape functions
    %       dn - matrix of shape functions derivatives
    %
    % **************************************************************************

    % VARIABLE INITIALISATION
        n = zeros(1, nnodel) ;
        dn = zeros(ndim, nnodel) ;
        
    % SHAPE FUNCTIONS
        n(1,1) = 0.25 * (1 - xi) * (1 - eta) ;
        n(1,2) = 0.25 * (1 + xi) * (1 - eta) ;
        n(1,3) = 0.25 * (1 + xi) * (1 + eta) ;
        n(1,4) = 0.25 * (1 - xi) * (1 + eta) ;
    
    % SHAPE FUNCTION DERIVATIVES
        dn(1,1) = -0.25 * (1 - eta) ;
        dn(1,2) = +0.25 * (1 - eta) ;
        dn(1,3) = +0.25 * (1 + eta) ;
        dn(1,4) = -0.25 * (1 + eta) ;
        dn(2,1) = -0.25 * (1 - xi) ;
        dn(2,2) = -0.25 * (1 + xi) ;
        dn(2,3) = +0.25 * (1 + xi) ;
        dn(2,4) = +0.25 * (1 - xi) ;
end