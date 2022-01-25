# nsub
nsub    nonsingular submatrix

Assuming that A is an m-by-n matrix and rank(A)>=r, [R] = nsub(A,r) returns a vector R with r elements in the range (1:m), such that the rows of A(R,:) are linear independent.

[R,C] = nsub(A,r) returns a vector R with r elements in the range (1:m) and a vector C with r elements in the range (1:n), such that A(R,C) is a nonsingular submatrix of A.

Note: If r>rank(A) nsub returns an error message.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Example:

r = 3;

A = [-1 -1 1 1 -5; -1 -1 1 0 -5; 0 0 0 1 0; 2 1 1 1 1];

% Note: row3 = row1 - row2

[R,C] = nsub(A,r);

% Result:

   % R = [2 4 1];
   
   % C = [5 1 4];
   
   % A(R,C) = [-5 1 0; 1 2 1; -5 -1 1];
   
   % rank(A(R,C)) = 3;
   
   
Author of the code:  Gabriel Ponte. Collaborative work with:  Marcia Fampa, Jon Lee and Luze Xu.  
Related work can be found in: Fampa, M., Lee, J., Ponte, G. et al. Experimental analysis of local searches for sparse reflexive generalized inverses. J Glob Optim 81, 1057–1093 (2021). https://doi.org/10.1007/s10898-021-01087-y

