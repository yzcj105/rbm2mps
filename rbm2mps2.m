function [ mps ] = rbm2mps2( W, b, a )
%{
FUNCTION RBM2MPS2 v1.0 18/JAU/2017

Article
Jing Chen, Song Cheng, Haidong Xie, Lei Wang, and Tao Xiang
On the Equivalence of Restricted Boltzmann Machines and Tensor Network States
[arXiv.17799849 ](https://arxiv.org/submit/1779849)

A Matlab code for rbm2mps on GitHub
https://github.com/yzcj105/rbm2mps


Authors
  Haidong Xie
  Jing Chen     E-mail: yzcj105@126.com
  Song Cheng
  The Institute of Physics, Chinese Academy of Sciences

A Matlab code for rbm2mps. 
Fig.4 in Sec.II of the article

Parameters:( See Eq.1 )
    a (Input) Double precision array of length Nv. 
            a(k) means the bias of v_i
    b (Input) Double precision array of length Nh.
            b(k) means the bias of h_i
    W (Input) Double precision matrix by Nv * Nh. 
            W(i,j) means the weight in the connection between  v_i and h_j
    mps (output) Cell of length Nv, Each element of cell mps{k} is a tensor
            of mps rep... with three bands: left,right,and phy.
    Nv The number of visible units.
    Nh The number of hidden uits.


Remarks:
    



%}

% Sub function rbm2data
% Get the index and position of hidden units form the given rbm.
[Lindex,Rindex,Containhi]=rbm2data(W);

% Sub function index2tensor
% Get mps (result) form the given parameters.
mps=index2tensor(a,b,W,Lindex,Rindex,Containhi);

end

function [ Lindex,Rindex,Contain_hi ]=rbm2data(W)
%{
Sub function rbm2data
Get the index and position of hidden units form the given rbm.
Input: W: a matrix by Nv * Nh
            W(i,j) means the weight in the connection between  v_i and h_j
Output: Lindex, Cell of length Nv
            Lindex{k} means the left index of tensor k.
        Rindex, Cell of length Nv
            Rindex{k} means the right index of tensor k.
        Contain_hi, Cell of length Nv
            Contain_hi{k} means the tensor k contains the hidden units list
            in Contain_hi{k}.
%}
[Nv,Nh]=size(W);

Connect_=abs(W*W')>eps;
%Connect_ is nv*nv logical matrix, Connect_(i,j) means visible units i and
%j is connected.
banddim=zeros(1,Nv);
bandindex=cell(1,Nv);
for k=2:Nv
    % Calculate each band of mps
    % k, means band between tensor:k-1 and tensor k
    index_l=1:k-1;
    % index_l: the visible units in the left of band above, which could
    % connect to the right part of band above.
    index_r=k:Nv;
    % the oppesite meanning of index_l.
    flag_find=0;

    for k_banddim=1:min(k-1,Nv+1-k)
        %k_banddim shows the band dimension of the calculate band, we wish
        %it as small as possiable. band dimension D=2^k_banddim
        all_choose=nchoosek(1:Nv,k_banddim);
        % all possiable chooses with the given k_banddim
        
        % we search the first result and just use it in this code, the
        % other result may cause degenerate sulotions which we do not
        % mention about.
        for k_=1:size(all_choose,1)
            k_choose=all_choose(k_,:);
            %k_choose, the choosen index
            index_lnew=setdiff(index_l,k_choose);
            index_rnew=setdiff(index_r,k_choose);
            Con_Mat=Connect_(index_lnew,index_rnew);
            if all(~Con_Mat(:))
                % if all element in Con_Mat are zeros, the left part and
                % the right part by this choosen index become direct
                % product. for detail see eq.8 and fig.4 in paper
                flag_find=1;
                banddim(k)=k_banddim;
                %band dimension of the band
                bandindex{k}=k_choose(:);
                %the choosen index
%                 disp(banddim(k));
%                 disp(bandindex{k});
                break;
            end
        end
        if flag_find==1
            break;
            %quit the loop, and fond sulotion
        end
        
        
    end
end
%disp
%disp(banddim);
% for k=2:Nv
    %disp(bandindex{k});
% end

% locat hi
Contain_hi=cell(1,Nv);
for k_h=1:Nh
    % find where to settle down the hidden units
    flag=0;
    for k_v=1:Nv
        % get the index of left,right,physical.
        if k_v==1
        index_l=[];
        index_r=bandindex{k_v+1};
        index_p=k_v;
        elseif k_v==Nv
            index_l=bandindex{k_v};
            index_r=[];
            index_p=k_v;
        else
            index_l=bandindex{k_v};
            index_r=bandindex{k_v+1};
            index_p=k_v;
        end
        index_tot=union(index_l,index_r);
        index_tot=union(index_tot,index_p);
        % index_tot means all the index of the tensor k
        index_hi=find(abs(W(:,k_h))>eps);
        % index_hi means the index linked to of h_i.
        if isempty(setdiff(index_hi,index_tot))
            %index_hi belong to index_tot
            flag=1;
            Contain_hi{k_v}(end+1)=k_h;
            break;
        end
    end
    if flag==0
        disp('not match, hi!');
        % show error message, there is a h_i, that no tensor could settle
        % down it.
    end
end
% for k=1:Nv
    %disp(Contain_hi{k});
% end

Lindex=cell(1,Nv);
Lindex(2:Nv)=bandindex(2:Nv);
Rindex=cell(1,Nv);
Rindex(1:Nv-1)=bandindex(2:Nv);
% for output
end

function mps=index2tensor(a,b,W,Lindex,Rindex,Containhi)
%{
Sub function index2tensor
Get mps (result) form the given parameters.
Input a,b,W, the same with the parent function
Input Lindex,Rindex,Containhi, see the output of function rbm2data
Output mps, the mps of the output in the parent function
%}
Nv=length(a);
Nh=length(b);
%{
check! 
it is must satisfied the equations below, if the function could function
good.
[Nv,Nh]==size(W);
Nv==length(Lindex)==length(Rindex)==length(Containhi);
%}

mps=cell(1,Nv);
for k=1:Nv
    % calculate each tensor one by one
    mps{k}=zeros(2^length(Lindex{k}),2^length(Rindex{k}),2);
    TOTindex=[Lindex{k};Rindex{k};k];
    %TOTindex , the index of tensor
    [uni_index,ia,ic]=unique(TOTindex(:));

    for k_nv=1:2^length(uni_index)
        % for each element of tensor, calculate the weight.
        idunique_index=myind2sub(length(uni_index),k_nv);
        % the value of each index, example:idunique_index=[0 1 1 0 0 1],Binary number
        exp_const=exp(a(k)*(idunique_index(uni_index==k)));
        % weight term of visible part.
        
        id_index=idunique_index(ic);
        allk_nv=sum((id_index(:)).*2.^(0:length(TOTindex)-1)')+1;
        % allk_nv: the tensor index
        cal_v=idunique_index;
        cal_w=W(uni_index,Containhi{k});
        cal_b=b(Containhi{k});
        
        if isempty(Containhi{k})
            mps{k}(allk_nv)=exp_const;
        else
            expindex=cal_v'*cal_w+cal_b;
            eexp=exp(expindex);
            eexp1=1+eexp;
            mps{k}(allk_nv)=prod(eexp1)*exp_const;
            %for more detail, see equation 2 in paper
        end
    end
end
end

function sub = myind2sub(n,ndx)
% change num from Decimal number to Binary number
bin_num=dec2bin(ndx-1,n);
sub=str2num(bin_num(:));
sub=sub(:);
end