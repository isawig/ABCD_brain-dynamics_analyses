%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Code by Selen Atasoy
%
% projects the data to CH basis
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [P] = surfFMRI_projectCH(data, CH, functions, timepoints, type)

if ~exist('functions', 'var')|| isempty(timepoints)
    functions = 1:size(CH,2);
end

if ~exist('timepoints', 'var')|| isempty(timepoints)
    timepoints = 1:size(data,2);
end

if ~exist('type', 'var')
    type = 'cos';
end

switch type
    
    %% compute the cosine distance (this doesnt take amplitude into account)
    case 'cos'
        P = acosd(1 - pdist2(data(:, timepoints)',CH(:,functions)', 'cosine'));

        % CH can be inverted, so it is the absolute value of the distances that matters
        P = min(P, 180-P);

        % normalize
        P = (90-P)./90;
        
    %% compute the dot product
    case 'dot'
        P = zeros(length(timepoints), length(functions));
        for i=1:length(timepoints)
            pattern = data(:, timepoints(i));
            P(i,:) = connProject2Basis(pattern, CH(:,functions));
        end
        
    %% do the dot product efficiently
    case 'dotefficient'
        
        P = data' * CH'.';
        
end

