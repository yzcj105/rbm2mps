function [ A1, L, Coef] = MPS_Canonicalize( A )
% finite size mps canonicalize
% Writen by Jing Chen, Haidong Xie, Song Chen
% Matrix Product State Representations
%D. Perez-Garcials, F. Verstraete, M.M. Wolf, J.I. Cirac
%  	arXiv:quant-ph/0608197
% INPUT: A : cell of 3 leg tensor A
% OUTPUT: A1 :cell of 3 leg tensor
%           The MPS is the tensor product of A1
%          L : cell of vector,
%             the entanglement spectrum on each bond
%        Coef: The coefficient factor for the MPS

%------QR-from right to left---------------%
len = numel( A );
A1 = cell( len,1 );
L = cell( len,1 );
R = 1; % A{len} is D * 1 * d ;
for k = len : -1 :1
    [ R, A1{k} ] = MPS_qr(A{k},R);
end
A1{1} = A1{1} * sign( R );
Coef = abs(R); % number
%----------svd from left to right----------%
L{1} = 1;
V = 1;
for k = 1:len-1
    [ A1{k}, L{k+1}, V ] = MPS_svd(L{k},V,A1{k});
end
[ A1{len} ]  = MPS_svd(L{len},V,A1{len});
end
function [ R1, A1 ] = MPS_qr( A,R )
[ dm ] = size( A,3 );
Tmp1 = tensor_product( 'AB1',A,'Ab1',R,'Bb');
dB1 = size( Tmp1,2 );
Tmp1 = tensor_reshape( Tmp1,'ABm','Bm','A' );
[ Q, R1 ] = qr( Tmp1,0 );
dA1 = numel(Q)/dm/dB1;
Tmp2 = reshape( Q, [dB1,dm,dA1 ] );
A1 = permute( Tmp2 ,[ 3 1 2 ]);
end
function [ A1, L1, V1 ] = MPS_svd( L,V,A )
LV = diag(L)*V';
[ LVA ] = tensor_product('aB1',LV,'aA',A,'AB1');
[ LVA ] = tensor_reshape( LVA, 'aB1', 'a1','B');
[ ~, L1, V1 ] = svd(LVA,'econ');
L1 = diag(L1);
% smaller than the machine error has been truncated
idx = find( L1 > eps,1,'last' );
L1 = L1(1:idx);
V1 = V1(:,1:idx);
%A1 = tensor_n_product('ab1',V,'Aa',A,'AB1',V1,'Bb');
Tmp = tensor_product('Ab1',A,'AB1',V1,'Bb');
A1 = tensor_product('ab1',Tmp,'Ab1',V,'Aa');
end






