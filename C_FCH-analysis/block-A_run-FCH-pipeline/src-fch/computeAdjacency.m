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
%           PARAM       ~   value of the previous parametre:
%                           number of nearest neighbours to find or the
%                           epsilon threshold
% OUTPUT:   W           ~   adjacency matrix
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [W] = computeAdjacency(D, TYPE, PARAM, sigma)

if (nargin < 3) | (strcmp(TYPE,'nn') & strcmp(TYPE,'epsballs')) | ~isreal(PARAM)
  
  disp(sprintf('ERROR: Too few arguments given or incorrect arguments.\n'));
  disp(sprintf('USAGE:\n A = laplacian(DATA, TYPE, PARAM)'));
  disp(sprintf('DATA - the data matrix. Data points are rows.'));
  disp(sprintf('Nearest neigbors: TYPE =''nn''    PARAM = number of nearest neigbors')); 
  disp(sprintf('Epsilon balls: TYPE =''epsballs''    PARAM = redius of the ball\n'));
  return;
end

n = size(D,1);
disp (sprintf ('DATA: %d points in %d dimensional space.',n,size (D,2)));

switch TYPE
 case {'nn'}
  disp(sprintf('Creating the adjacency matrix. Nearest neighbors, N=%d.', PARAM)); 
 case{'eps', 'epsballs'} 
  disp(sprintf('Creating the adjacency matrix. Epsilon balls, eps=%f.', PARAM));
end;

  
A = sparse(n,n);
% A = ones(n,n);
% A = A * max(max(D));

if (strcmp(TYPE,'nn'))
    
    %%% set all zero entries to maximum value %%%
%     [d_i, d_j, d_val] = find(D);
%     min_val = min(d_val)    
%     max_val = max(d_val)
%     z_matrix = max_val*2*(D==0);
%     D = D + z_matrix;
    
    %% find the minimum value of z
    [Z,I] = sort ( D,2);
    disp('finding the minimum value of z');
    for i=1:n
        for j=2:PARAM+1
            
            %% eliminate zeros - set to minimum
            if (Z(i,j)==0)
                display(['Setting Z(', int2str(i), ',', int2str(j), ') = ', int2str( 0 ), ')' ]);
%                 Z(i,j) = min_val;
            end
            
            %% original version 
            A(i,I(i,j))= Z(i,j); 
            A(I(i,j),i)= Z(i,j); 
        end; 
    end;
    
else
     A = D.*(D<PARAM);
end;

[A_i, A_j, A_v] = find(A);

if (sigma==0)

    for i = 1: size(A_i)  
        W(A_i(i), A_j(i)) = 1;
    end;
else

    if (strcmp(sigma, 'var'))
        display(['Variance of the dataset is ', num2str(var(A_v))]);
        for i = 1: size(A_i)  
            W(A_i(i), A_j(i)) = exp(-A_v(i)^2/(3*var(A_v))); %% my version
        end;

    else
        for i = 1: size(A_i)  
            W(A_i(i), A_j(i)) = exp(-A_v(i)^2/sigma); %% my version
         end;
    end;

end;



 

