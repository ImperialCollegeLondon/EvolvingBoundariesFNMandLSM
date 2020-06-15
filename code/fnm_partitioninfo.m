% **************************************************************************
%  MODELLING EVOLVING BOUNDARIES USING FNM AND LSM (April 2020)
% **************************************************************************
%
%
%   fnm_partitioninfo.m
%
%
%   FUNCTIONS:
%       - fnm_partitioninfo - outputs the number of sub elements for a given
%                             partitionioning scheme
%
%
% **************************************************************************
%                   developed by Rui O.S.S. da Costa and Silvestre T. Pinho
%                                      Imperial College London, Aeronautics
% **************************************************************************


function nsub = fnm_partitioninfo( elstatid )

    % **************************************************************************
    %
    %   INPUTS:
    %       elstatid - flag that identifies the type of partitioning
    %
    %   OUTPUTS:
    %       nsub - number of sub elements
    %
    % **************************************************************************

    % NO INTERSECTION
        if elstatid == 0
            nsub = 1;
            
    % INTERSECTION CASE: SINGLE CORNER
        elseif elstatid == 11 || elstatid == 12 || elstatid == 13 || elstatid == 14
            nsub = 6;
            
    % INTERSECTION CASE: THROUGH ELEMENT
        elseif elstatid == 21 || elstatid == 22
            nsub = 2;
            
    % INTERSECTION CASE: DOUBLE CORNER
        elseif elstatid == 31 || elstatid == 32
            nsub = 8;
        end
end