function T = tensor_reshape(T,idx_all,varargin)
%{
FUNCTION TENSOR_RESHAPE v1.0 18/JAU/2017

Article
Jing Chen, Song Cheng, Haidong Xie, Lei Wang, and Tao Xiang
On the Equivalence of Restricted Boltzmann Machines and Tensor Network States
[arXiv.17799849 ](https://arxiv.org/submit/1779849)
A Matlab code for rbm2mps on GitHub
https://github.com/yzcj105/rbm2mps

Authors
  Jing Chen     E-mail: yzcj105@126.com
  Song Cheng
  Haidong Xie
  The Institute of Physics, Chinese Academy of Sciences

A Matlab code for permute and reshape of tensors. 
For example a 4 index tensor A can be transformed into a matrix by ¡°permute¡± and ¡°reshape¡± 
B( ac,bd ) = A(a,b,c,d); 
MATLAB code: 
B = tensor\_reshape( A,'abcd','ac','bd')

Parameters:
    T (Input)       Input tensor
    idx_all (Input) The indices for all the bonds of T
    varargin(Input) All the other inputs are the indices after the reshape.
    T (output)      The tenosr after permutation and reshape of indices


Remarks:
    The function is often used together with tensor_product



%}

sizeT = size(T);
sizeT(end+1:length(idx_all)) = 1;
sz_part = zeros(1,nargin-2);
%idx_part = cell(1,nargin-2);
reshape_sz = ones(1,max(nargin-2,2));
reshape_idx = [];
for k = 1:nargin-2
    part = varargin{k};
    sz_part(k) = length(part);
    temp = zeros(1,sz_part(k));
    for s = 1:sz_part(k)
        [~,temp(s)] = find(idx_all == part(s),1);
    end
    %idx_part{k} = temp;
    reshape_sz(k) = prod(sizeT(temp));
    reshape_idx = [reshape_idx,temp];
end
if any(reshape_idx ~= 1:length(idx_all))
    T = permute(T,reshape_idx);
end
T = reshape(T,reshape_sz);

