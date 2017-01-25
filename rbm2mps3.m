function [ mps ] = rbm2mps3( W, b, a )
%{
FUNCTION RBM2MPS3 v1.0 18/JAU/2017

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
            of mps rep... with three bonds: left,right,and phy.
    Nv The number of visible units.
    Nh The number of hidden uits.


Remarks:
    rbm2mps3 is the most general and robust program.
    In rbm2mps1 (method1), we give a plainly method, which is very sample
    and easy to understand. But the calulation cost and bond-dimension is
    really high.
    In method2, we give a more economic method, which could give a lower
    bond-dimension.
    In this method (method3), we give the most general method, which could
    give the lowest bond-dimension: the same with given after canonicalize.



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
            Example: Contain_hi{k}(k_h) = h_index. 
                If h_index<= Nh, means this h_index is located in the
                center of tensor k.
                If h_index>Nh. means this h_index is located at the bond of
                tensor k.
%}
[Nv,Nh]=size(W);

bonddim=zeros(1,Nv);
bondindex=cell(1,Nv);
for k=2:Nv
    % Calculate each bond of mps
    % k, means bond between tensor:k-1 and tensor k
    index_l=1:k-1;
    % index_l: the visible units in the left of bond above, which could
    % connect to the right part of bond above.
    index_r=k:Nv;
    % the oppesite meanning of index_l.
    flag_find=0;

    for k_bonddim=0:min(k-1,Nv+1-k)
        %k_bonddim shows the bond dimension of the calculate bond, we wish
        %it as small as possiable. bond dimension D=2^k_bonddim
        %Here, k_bonddim==0 means the tensor of left and right has no
        %relation. In this case, the left and right could write as two
        %sperate, independence, irrelevant part.
        
        all_choose=nchoosek(1:(Nv+Nh),k_bonddim);
        % all possiable chooses with the given k_bonddim
        % here, we can choose both visible and hidden index.
        % If the choosen index <=Nv, means we choose visible index
        % If the choosen index >Nv, means we choose the hidden index =
        % choosen-index - Nv.
        
        % we search the first result and just use it in this code, the
        % other result may cause degenerate sulotions which we do not
        % mention about.
        for k_=1:size(all_choose,1)
            k_choose=all_choose(k_,:);
            k_choose_v=k_choose(k_choose<=Nv);
            k_choose_h=k_choose(k_choose>Nv)-Nv;
            %k_choose, the choosen index
            %k_choose_v, the visible index of k_choose
            %k_choose_h, the hidden index of k_choose
            index_lnew=setdiff(index_l,k_choose_v);
            index_rnew=setdiff(index_r,k_choose_v);
            index_hnew=setdiff(1:Nh,k_choose_h);
            %index_*new means the index that git rid of the choosen index.
            W_mat_1=W(index_lnew,index_hnew);
            W_mat_2=W(index_rnew,index_hnew);
            Con_Mat=W_mat_1*W_mat_2';
            
            if all(~Con_Mat(:))
                % if all element in Con_Mat are zeros, the left part (index_l) and
                % the right part (index_r) by fixed this choosen index become direct
                % product. for detail see eq.8 and fig.4 in paper
                flag_find=1;
                bonddim(k)=k_bonddim;
                %bond dimension of the bond
                bondindex{k}=k_choose(:);
                if k_bonddim==0;
                    bondindex{k}=[];
                end
                %the choosen index
%                 disp(bonddim(k));
%                 disp(bondindex{k});
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
%disp(bonddim);
% for k=2:Nv
    %disp(bondindex{k});
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
        index_r=bondindex{k_v+1};
        index_p=k_v;
        elseif k_v==Nv
            index_l=bondindex{k_v};
            index_r=[];
            index_p=k_v;
        else
            index_l=bondindex{k_v};
            index_r=bondindex{k_v+1};
            index_p=k_v;
        end
        index_tot=union(index_l,index_r);
        index_tot=union(index_tot,index_p);
        % index_tot means all the index of the tensor k, sum of left,
        % right, physical.
        
        % situation 1, hi in bond of tensor
        if ismember(k_h+Nv,index_tot)
            %index_hi belong to index_tot
            flag=1;
            Contain_hi{k_v}(end+1)=k_h+Nh;
            % with +Nh, means k_h index is a bond index
%             disp('use method 3')
%             disp(Contain_hi{k_v})
            break;
        end
        % situation 2, hi in center of tensor
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
Lindex(2:Nv)=bondindex(2:Nv);
Rindex=cell(1,Nv);
Rindex(1:Nv-1)=bondindex(2:Nv);
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
cal_need_wjk=abs(W)>eps;
%If cal_need_wjk(i,j)~=0 means the connection of visible v_i and hidden h_j
%is exist and need to be calculated.
for k=1:Nv
    % calculate each tensor one by one
    mps{k}=zeros(2^length(Lindex{k}),2^length(Rindex{k}),2);
    TOTindex=[Lindex{k};Rindex{k};k];
    %TOTindex , the index of tensor
    [uni_index,ia,ic]=unique(TOTindex(:));
    %
    cal_need_v=uni_index(uni_index<=Nv);
    %cal_need_v means visible index in the tensor index
    cal_need_h=uni_index(uni_index>Nv)-Nv;
    %cal_need_h means hidden index in the tensor index
    cal_need_k=zeros(Nv,Nh);
    cal_need_k(cal_need_v(:)',cal_need_h(:)')=cal_need_wjk(cal_need_v(:)',cal_need_h(:)');
    %If cal_need_k(i,j)~=0 means the connection of visible v_i and hidden h_j
    %has not calculated and need to be calculated in this tensor k
    cal_need_wjk(cal_need_v(:)',cal_need_h(:)')=0;
    %set the connection that has been calculated to 0,
    
    %disp(cal_need_k)
    for k_nv=1:2^length(uni_index)
        % for each element of tensor, calculate the weight.
        idunique_index=myind2sub(length(uni_index),k_nv);
        % the value of each index, example:idunique_index=[0 1 1 0 0 1],Binary number
        
        vk=idunique_index(uni_index==k);
        exptop=a(k)*vk;
        % the bond weight of visible unit
        
        calout_hi=Containhi{k}(Containhi{k}>Nh)-Nh;
        %calout_hi means the index which weight of hidden index need to be calculated
        if ~isempty(calout_hi)
            for k_hi=1:length(calout_hi)
                hk(k_hi)=idunique_index(uni_index==(calout_hi(k_hi)+Nv));
                exptop=exptop+b(calout_hi(k_hi))*hk(k_hi);
            end
        end
        
        for k_v=cal_need_v(:)'
            for k_h=cal_need_h(:)'
                if cal_need_k(k_v,k_h)
                    v_need=idunique_index(uni_index==k_v);
                    h_need=idunique_index(uni_index==((k_h)+Nv));
                    exptop=exptop+v_need*W(k_v,k_h)*h_need;
                    % calculate each connection of cal_need_k, to the
                    % weight term.
                end
            end
        end
        exp_const=exp(exptop);
        % weight term of index part.
        
        id_index=idunique_index(ic);
        allk_nv=sum((id_index(:)).*2.^(0:length(TOTindex)-1)')+1;
        % allk_nv: the tensor index
        
        calin_hi=Containhi{k}(Containhi{k}<=Nh);
        %calculate weight for the hi located at the center of tensor k
        if isempty(calin_hi)
            mps{k}(allk_nv)=exp_const;
        else
            cal_v=idunique_index(uni_index<=Nv);
            cal_w=W(uni_index(uni_index<=Nv),calin_hi);
            cal_b=b(calin_hi);
            % v,w,b same with equation 2 in paper
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