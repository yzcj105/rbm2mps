# Transforming an RBM into an MPS

These are demonstration codes for the Sec.II of [arXiv.17799849 ](https://arxiv.org/submit/1779849). You are free to use them. Please kindly cite the paper 
- Jing Chen, Song Cheng, Haidong Xie, Lei Wang, and Tao Xiang, "On the Equivalence of Restricted Boltzmann Machines and Tensor Network States", [arXiv.17799849 ](https://arxiv.org/submit/1779849)


We offer two versions of the code. Both are written in MATLAB. 
* rbm2mps.m
Constructs the MPS with the bond dimension D equals the number of cut connections at the bipartition. See Fig.2.
    * Input:	
        * W:   n_v by n_h weight matrix W
    	* a:   bias vector of n_v for visible units
	    * b:   bias vector of n_h for hidden units
	    * bpos: the vector telling which piece the hidden units belongs to. See Fig.2 
    * Output: 
        * mps: a cell of of mps tensors

* rbm2mps2.m
Constructs the MPS with bond dimension D equals to the size of min(A_{1},B_{1}). The internal bond are just copies of the physical bonds. See Fig.4. 
    * Input:      
      * W:  n\_v by n\_h weight matrix W
      * a:  bias vector of n_v for visible units
      * b:  bias vector of n_h for hidden units
    * Output: 
      * mps: a cell of of mps tensors

The bond dimension D of rbm2mps2.m is smaller than that of rbm2mps.m, but it is not the smallest. In rbm2mps2.m, the space of the internal bond is only visible units. Actully we are not limited to this case, we can set the internal bond the copies of visible or hidden units. We will add the version latter. 

## Auxillary tensor programs ##
* MPS\_Canonicalize.m
Canonicalize a finite MPS and give us the entanglement spectrum in each MPS bond 

* tensor\_product.m 
Calculate the contraction of two tensors and permute the index with the given order.
For example 
C(a,c,m) = sum_{c}A_{a,b,m)X_{c,b}
MATLAB code:
`C = tensor_product('acm',A,'abm',X,'cb');`

* tensor\_reshape.m
For example a 4 index tensor A can be transformed into a matrix by "permute" and "reshape" 
B( ac,bd ) =  A(a,b,c,d);
MATLAB code:
`B = tensor\_reshape( A,'abcd','ac','bd')`


## Test program ##
* Example.m
Take the RBM in Fig.1(a) as an example and construct MPS by two different ways as Fig.2 and Fig.4 in Sec.II . The bond diemensions are shown, consistant with Fig.2(c) and Fig.4(c). The MPS are exactly the same in their canonical form. It helps to simplify an RBM as in Sec.V, Fig.8. 
