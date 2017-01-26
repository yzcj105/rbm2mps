function [mps] = rbm2mps(W, b, a, bpos)
%{
FUNCTION RBM2MPS v1.0 18/JAU/2017

Article
Jing Chen, Song Cheng, Haidong Xie, Lei Wang, and Tao Xiang
On the Equivalence of Restricted Boltzmann Machines and Tensor Network States
[arXiv:1701.04831](http://arxiv.org/abs/1701.04831)

A Matlab code for rbm2mps on GitHub
https://github.com/yzcj105/rbm2mps


Authors
  Song Cheng
  Jing Chen     E-mail: yzcj105@126.com
  Haidong Xie
  The Institute of Physics, Chinese Academy of Sciences

A Matlab code for rbm2mps. 
Fig.2 and Fig.3 in Sec.II of the article

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


nv = length(a);
nh = length(b);
% each piece related to one of MPS's local matrix
pieces = cell(nv,1);
% N record how many hiddens each visible connected with
N = ones(nv,1);
maxInd = 1;
% define indices of Lambda_v and Lambda_h
L_v_ind = zeros(nv,nh);
L_h_ind = L_v_ind;
for k = 1:nv*nh
    if( abs(W(k)) > eps )
        L_v_ind(k) = maxInd;
        L_h_ind(k) = maxInd +1;
        maxInd = maxInd + 2;
    end
end
% build tensor Lambda_v as Eq.5
for k = 1:nv
    index = L_v_ind(k,:);
    % to distinguish with vitual bond, physical bond indexes are negative
    pieces{k}{1} = struct('type','L','bias',a(k),'order',sum(abs(W(k,:))>eps)+1,'index',[ index(index>0),-k]);
end
% build tensor Lambda_h as Eq.6
for k = 1:nh
    % tell N the information of h2v mapping
    N(bpos(k))= N(bpos(k))+1;
    index = L_h_ind(:,k)';
    %put all hiddens connecting with a same visible into a same piece
    pieces{bpos(k)}{N(bpos(k))} = struct('type','L','bias',b(k),'order',sum(abs(W(:,k))>eps),'index',index(index>0));
end

% calculate local matrix M
% bigind is MPS's vitual bond index
bigind = cell(nv+1,1);
% this for-for circles will cut nonlocal long-term connection into 
% several local short-term connection see Fig.3
% in equation M = P*Q in Fig3, we choose P = M, and Q = I.
for i = 1:nh
    for j = 1:nv
        if (abs(W(j,i))>eps)
            % we want cutting connection always from right to left
            % if visible in the left of hidden
            if( bpos(i) > j)
                %first (MPS) piece of this connection
                k = bpos(i);
                %first index of this connection
                ind1 = L_h_ind(j,i);
                %last piece of this connection
                end_piece = j;
                %last index of this connection
                end_ind = L_v_ind(j,i);
            % visible in the right of hidden
            else
                k = j;
                ind1 = L_v_ind(j,i);
                end_piece = bpos(i);
                end_ind = L_h_ind(j,i);
            end
            % which is false only when bpos(i) == j
            % and bpos(i) == j means it's no longer a nonlocal connection
            while( k > end_piece)
                N(k)= N(k)+1;
                ind2 = maxInd;
                maxInd = maxInd + 1;
                % matrix product an identity on first piece, did nothing.
                % then direct product an identity matrix on following
                % pieces before we get to the end piece
                pieces{k}{N(k)} = struct('type','I','index',[ ind2, ind1]);
                bigind{k} = [bigind{k},ind2];
                ind1 = ind2;
                k = k - 1;
            end
            N(k)= N(k)+1;
            % put weight matrix on the end index
            pieces{k}{N(k)} = struct('type','M','Weight',W(j,i),'index',[ end_ind, ind1 ]);
        end
    end
end
% constructing mps
mps = cell(nv,1);
for k = 1:nv
    [ mps{k} ] = contract_tensor(pieces{k},bigind{k},bigind{k+1},-k);
end


function [ T ] = contract_tensor( tens,left_ind,right_ind,physical_ind )
% contract visible tensor in the same piece
% see Fig.2
temp = constructTensor(tens{1});
temp_idx = tens{1}.index;

for k = 2:length(tens)
    T2 = constructTensor(tens{k});
    [temp,temp_idx] = tensor_product(temp,temp_idx,T2,tens{k}.index);  
end
% reshape tensor into MPS pieces (3rd order tensor)
T = tensor_reshape( temp, temp_idx,left_ind,right_ind,physical_ind );

function T = constructTensor( ten )
switch ten.type
    % eq. (3-5) of paper
    case 'I' 
    % Fig.3 We choose Q = I.
        T = eye(2);
    case 'M'
        T = [ 1,1;1 exp(ten.Weight)];
    case 'L'        
        T = zeros([2*ones(1,ten.order),1,1]);
        T(1) = 1;
        T(end) = exp(ten.bias);
end



 

