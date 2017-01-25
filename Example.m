function Example()
%{
FUNCTION EXAMPLE v1.0 18/JAU/2017

This code construct the MPS in different ways Fig.2 and Fig.4 in Sec.II.
And verify their equivalence under canonical form.  

Article
Jing Chen, Song Cheng, Haidong Xie, Lei Wang, and Tao Xiang
On the Equivalence of Restricted Boltzmann Machines and Tensor Network States
[arXiv:1701.04831](http://arxiv.org/abs/1701.04831)

A Matlab code for rbm2mps on GitHub
https://github.com/yzcj105/rbm2mps


Authors
  Song Cheng
  Haidong Xie
  Jing Chen     E-mail: yzcj105@126.com
  The Institute of Physics, Chinese Academy of Sciences


Parameters:( See Eq.1 )
    a (Input) Double precision array of length Nv. 
            a(k) means the bias of v_i
    b (Input) Double precision array of length Nh.
            b(k) means the bias of h_i
    W (Input) Double precision matrix by Nv * Nh. 
            W(i,j) means the weight in the connection between  v_i and h_j
    bpos: 1*nh(int) (bias_position),tell which piece the hidden units belongs
             to in the cut, see Fig.2(b) bpos = [ 2 3 4 5];
    mps (output) Cell of length Nv, Each element of cell mps{k} is a tensor
            of mps rep... with three bands: left,right,and phy.
    Nv The number of visible units.
    Nh The number of hidden uits.


Remarks:
    



%}
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% generate an RBM as Fig.1(a) with 6 visible units and 4 hidden units
a=rand(1,6);
b=rand(1,4);
W = [  0.43498            0            0            0
    -0.020515      0.49548            0            0
     -0.26821      0.46243     0.080192            0
     -0.10371     0.035067     0.030964     -0.48333
            0            0      0.40121      0.30092
            0            0            0     -0.35749];
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Construct MPS as Fig.(2,3)
% put each 
bpos = [ 2 3 4 5 ]; 
mps1 = rbm2mps(W,b,a,bpos);
% show the bond dimension
D1 = MPS_D(mps1);
fprintf('The bond dimensions of MPS in Fig.2 are\n');
disp(D1);
% Construct MPS as Fig.4
mps2 = rbm2mps2(W,b,a);
D2 = MPS_D(mps2);
fprintf('The bond dimensions of MPS in Fig.4 are\n');
disp(D2);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Canonical mps1 and mps2 and compare them.
% The redundency of MPS can be removed by Canonicalization,
% which is a important property when simplifing an RBM
% See Sec.V Fig.8 
[ mps1_cano L1 Coef1 ]= MPS_Canonicalize( mps1 );
[ mps2_cano L2 Coef2 ] = MPS_Canonicalize( mps2 );

ErrCoef = abs(Coef1-Coef2)/Coef2;
[ ErrL ErrMPS ] = CompareMPS( mps1_cano,L1,mps2_cano,L2);
fprintf('The difference between MPS by Fig.2 and Fig.4\n');
fprintf('The difference of the coefficient:\t%g\n',ErrCoef);
fprintf('The difference of the entanglement spectrum in each bond:\n');
disp(ErrL);
%If the entanglment is nearly degenerate, the difference may be large. 
% This comes from the gauge freedom in the nearly degenerate space.
fprintf('The difference of the MPS tensors:\n');
disp(ErrMPS);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
D1_cano = MPS_D(mps1_cano);
fprintf('The bond dimensions of canonicalized MPS are\n');
disp(D1_cano);
% Which is the same as Fig.4
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% In Fig.2(b) If we permute the order of the hidden units
 bpos = [ 2 5 3 4 ];
 mps3 = rbm2mps(W,b,a,bpos);
 [ mps3_cano, L3, Coef3 ] = MPS_Canonicalize(mps3);
 [ ErrL13 ErrMPS13 ] = CompareMPS( mps1_cano,L1,mps3_cano,L3);
 fprintf('Difference of MPS if we permute hidden units in Fig.2(b):\t%g\n',...
     max([ErrL13,ErrMPS13]));
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end
function Ds = MPS_D(MPS)
Ds = zeros(1,length(MPS)-1);
% collect the bond dimension of a MPS
for k = 1:length(MPS)-1
    Ds(k) = size(MPS{k},2);
end
end
function [ ErrL ErrMPS ] = CompareMPS( mps1,L1,mps2,L2)
% Compare the difference of two MPS
% ErrL is for the entanglement spectrum
% ErrMPS is for the MPS tensor
len = length(mps1);
ErrMPS = zeros(1,len); 
ErrL = zeros(1,len);
for k = 1:length( mps1 )
    % before comparing, we should remove the signs freedom;
    v1 = abs(mps1{k}(:)); 
    v2 = abs(mps2{k}(:));
    ErrMPS(k) = norm(v1-v2);
    ErrL(k) = norm(L1{k}-L2{k});
end
end
