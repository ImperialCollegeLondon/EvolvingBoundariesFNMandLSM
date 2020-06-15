% **************************************************************************
%  MODELLING EVOLVING BOUNDARIES USING FNM AND LSM (April 2020)
% **************************************************************************
%
%
%   fnm_subcon.m
%
%
%   FUNCTIONS:
%       - fnm_subcon - outputs the local sub element connectivities 
%                      for a given partitionioning scheme
%
%
% **************************************************************************
%                   developed by Rui O.S.S. da Costa and Silvestre T. Pinho
%                                      Imperial College London, Aeronautics
% **************************************************************************


function [scon] = fnm_subcon( elstatid, isub )

    % **************************************************************************
    %
    %   INPUTS:
    %       elstatid - flag that identifies the type of partitioning
    %       isub - index identifying a sub element
    %
    %   OUTPUTS:
    %       scon - vector of local connectivities for a given sub element
    %
    % **************************************************************************

    % INTERSECTION CASE: SINGLE CORNER
        if elstatid == 0

            scon = [1 2 3 4];

    % INTERSECTION CASE: BOTTOM LEFT CORNER
        elseif elstatid == 11

            if isub == 1
                scon = [1 5 8];
            elseif isub == 2
                scon = [5 9 8];
            elseif isub == 3
                scon = [5 2 9];
            elseif isub == 4
                scon = [2 3 9];
            elseif isub == 5
                scon = [8 9 4];
            elseif isub == 6
                scon = [9 3 4];
            end
    
    % INTERSECTION CASE: TOP LEFT CORNER
        elseif elstatid == 12

            if isub == 1
                scon = [1 9 8];
            elseif isub == 2
                scon = [1 2 9];
            elseif isub == 3
                scon = [2 3 9];
            elseif isub == 4
                scon = [8 9 7];
            elseif isub == 5
                scon = [8 7 4];
            elseif isub == 6
                scon = [9 3 7];
            end

    % INTERSECTION CASE: TOP RIGHT CORNER
        elseif elstatid == 13

            if isub == 1
                scon = [1 9 4];
            elseif isub == 2
                scon = [1 2 9];
            elseif isub == 3
                scon = [2 6 9];
            elseif isub == 4
                scon = [9 6 7];
            elseif isub == 5
                scon = [6 3 7];
            elseif isub == 6
                scon = [9 7 4];
            end

    % INTERSECTION CASE: BOTTOM RIGHT CORNER
        elseif elstatid == 14

            if isub == 1
                scon = [1 5 9];
            elseif isub == 2
                scon = [5 6 9];
            elseif isub == 3
                scon = [5 2 6];
            elseif isub == 4
                scon = [6 3 9];
            elseif isub == 5
                scon = [9 3 4];
            elseif isub == 6
                scon = [1 9 4];
            end

    % INTERSECTION CASE: THROUGH ELEMENT VERTICAL
        elseif elstatid == 21

            if isub == 1
                scon = [1 5 7 4];
            elseif isub == 2
                scon = [5 2 3 7];
            end

    % INTERSECTION CASE: THROUGH ELEMENT HORIZONTAL
        elseif elstatid == 22

            if isub == 1
                scon = [1 2 6 8];
            elseif isub == 2
                scon = [8 6 3 4];
            end
    
    % INTERSECTION CASE: DOUBLE CORNER, BOTTOM LEFT TOP RIGHT
        elseif elstatid == 31

            if isub == 1
                scon = [1 5 8];
            elseif isub == 2
                scon = [5 9 8];
            elseif isub == 3
                scon = [5 2 9];
            elseif isub == 4
                scon = [2 6 9];
            elseif isub == 5
                scon = [8 9 4];
            elseif isub == 6
                scon = [4 9 7];
            elseif isub == 7
                scon = [9 6 7];
            elseif isub == 8
                scon = [6 3 7];
            end   

    % INTERSECTION CASE: DOUBLE CORNER, TOP LEFT BOTTOM RIGHT
        elseif elstatid == 32

            if isub == 1
                scon = [1 9 8];
            elseif isub == 2
                scon = [1 5 9];
            elseif isub == 3
                scon = [5 6 9];
            elseif isub == 4
                scon = [5 2 6];
            elseif isub == 5
                scon = [8 7 4];
            elseif isub == 6
                scon = [8 9 7];
            elseif isub == 7
                scon = [9 3 7];
            elseif isub == 8
                scon = [9 6 3];
            end
        end
end 