% **************************************************************************
%  MODELLING EVOLVING BOUNDARIES USING FNM AND LSM (April 2020)
% **************************************************************************
%
%
%   fem_sftri3.m
%
%
%   FUNCTIONS:
%       - fem_sftri3 - calulates de shape functions and its derivatives
%                      for a trianglular isoparametric element with 3 nodes
%
%
% **************************************************************************
%                   developed by Rui O.S.S. da Costa and Silvestre T. Pinho
%                                      Imperial College London, Aeronautics
% **************************************************************************


function [n, dn] = fem_sftri3( ndim, nnodel, xi, eta)

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
        n(1,1) = 1 - xi - eta;
        n(1,2) = xi;
        n(1,3) = eta;
    
    % SHAPE FUNCTION DERIVATIVES
        dn(1,1) = -1;
        dn(1,2) = 1;
        dn(1,3) = 0;
        dn(2,1) = -1;
        dn(2,2) = 0;
        dn(2,3) = 1;

end