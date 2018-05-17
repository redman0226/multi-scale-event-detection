function A = tmg_tdm(A)
% This code is extracted (and simplified) from the TMG MATLAB toolbox
% For computing the term-document matrix using TF-IDF and cosine distance
% The input A is a word-document matrix
% Rows of A: words
% Columns of A: documents

%==================================================================    
%construct term-document matrix
%==================================================================    

%compute global weights
global_weights=log2(size(A, 2)./((A~=0)*ones(size(A, 2), 1)));

%compute local weights
A=double(A);

%form final term-document matrix
A=spdiags(global_weights, 0, size(A, 1), size(A, 1))*A;
A(isnan(A))=0;

%compute normalization factors and normalize term-document matrix
normalization_factors=ones(size(A, 2), 1);
for i=1:size(A, 2), 
    normalization_factors(i)=norm(A(:, i));
end
A=A*spdiags(1./normalization_factors, 0, size(A, 2), size(A, 2));