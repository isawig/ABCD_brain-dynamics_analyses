%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Pearson Lab UNSW
% Code modified by Selen Atasoy
%
% original code by Michail Belkin, Laplacian Eigenmaps
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% given the distance matrix computes the adjacency matrix by finding the
% k-nearest neighbours 
%
%
% INPUT:    
%           D           ~   distance matrix
%           TYPE        ~   type of the neighbourhood system, k-nn or eps
%                           balls
%           PARAM       ~   value of the previous parameter:
%                           number of nearest neighbours to find or the
%                           epsilon threshold
%           sigma       ~   is a the edges will be weighted
%                           0           - no weighting - binary adjacency
%                           'weight'    - edges weighted with the inverse
%                                       of the distances
%                           'var'       - Gaussian kernel with sigma = 3*variance
%                                       of the data
%                           otherwise - Gaussian kernel with the given sigma
%
% OUTPUT:   W           ~   adjacency matrix
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [A] = computeAdjacencyFastEfficient(dist, type, param, sigma)

%% check the input and set the default values
if ~exist('type', 'var') || isempty(type)
    display('Using k-nn neighbourhood as default!');
    type = 'nn';
end

if ~exist('param', 'var') || isempty(param)
    display('Setting k-nn = 6 as default!');
    param = 6;
end

if ~exist('sigma', 'var') || isempty(param)
    display('Setting sigma = 1 as default!');
    sigma = 1;
end

%% 
n = ceil(sqrt(2*length(dist)));
fprintf(['DATA: %d points in %d dimensional space.', n, n, '\n']);

switch type
    case {'nn'}
        fprintf(['Creating the adjacency matrix. Nearest neighbors, N=%d.', param, '\n']); 
    case{'eps', 'epsballs'} 
        fprintf(['Creating the adjacency matrix. Epsilon balls, eps=%f.', param, '\n']);
end;
  
A = sparse(n,n);

%%
if (strcmp(type,'nn'))
    
    % find the minimum value of z
    %[Z,I] = sort(D,2);
    
    for i=1:n
        
        %i
        el_array = setdiff(1:n, i);
        [idx, D] = pdistIndex(dist, i, el_array);
        [Z, I] = sort(D);
        
        for j=1:param
                 
            A( i, el_array(I(j)) ) = Z(j); 
            A( el_array(I(j)), i ) = Z(j); 
%             A(i,I(i,j))= Z(i,j); 
%             A(I(i,j),i)= Z(i,j); 
        end; 
        
    end;
    
else
     display('Efficient method with eps neighbourhood is not implemented yet');
     return;
end;

% [A_i, A_j, A_v] = find(A);


%%
switch sigma
    case 0
        display('Computing the adjacency matrix for combinatorial Laplacian..');
        %val = 1;
        A = double(A>0);

    case 'weight'
        display('Computing the adjacency matrix with the actual weights.. ');
        display('Normalizing the weights first.. ');
        
%         [A_i, A_j, A_v] = find(A);
%         val = max(A_v)./A_v;
        % efficient indexing %%
%         A = zeros(size(A));
%         indexMatrix             = [A_i, A_j];
%         indexCell               = num2cell(indexMatrix,1);
%         linearIndexMatrix       = sub2ind(size(A),indexCell{:});
%         A(linearIndexMatrix)    = val;
%         A = sparse(A);
        
        A(A>0) = max(max(A))./ A(A>0);
        
        
    otherwise
        display('Computing the adjacency matrix with exp(W_ij^2/sigma).. ');
%         [A_i, A_j, A_v] = find(A);
%         val = exp(-A_v.^2/sigma); 
        % efficient indexing %%
%         A = zeros(size(A));
%         indexMatrix             = [A_i, A_j];
%         indexCell               = num2cell(indexMatrix,1);
%         linearIndexMatrix       = sub2ind(size(A),indexCell{:});
%         A(linearIndexMatrix)    = val;
%         A = sparse(A);
        
    A = A.*(A>0);  
    A = -A.^2/sigma;
    A = exp(A);
    
    %A = exp(-A.^2/sigma).*(A>0);

end;



