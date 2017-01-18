# Transforming an RBM into an MPS

Demo code for the paper [arXiv.17799849 ](https://arxiv.org/submit/1779849). You are free to use these codes. Please kindly cite the paper 
- Jing Chen, Song Cheng, Haidong Xie, Lei Wang, and Tao Xiang, *On the Equivalence of Restricted Boltzmann Machines and Tensor Network States*, [arXiv.17799849 ](https://arxiv.org/submit/1779849)

We implement two approaches using Matlab 

* rbm2mps.m: The algorithm of Fig. 2. The MPS bond dimension is determined by the number of cut RBM connections. 
    * Input:	
        * W:   n_v by n_h weight matrix W
    	* a:   bias vector of n_v for visible units
        * b:   bias vector of n_h for hidden units
   	    * bpos: the vector telling which piece the hidden units belongs to. See Fig. 2.
    * Output: 
        * mps: a cell of MPS tensors

* rbm2mps2.m: The algorithm of Fig. 4. The MPS bond dimension is determined by min(|A_{1}|, |B_{1}|). 
    * Input:      
      * W:  n\_v by n\_h weight matrix W
      * a:  bias vector of n_v for visible units
      * b:  bias vector of n_h for hidden units
    * Output: 
      * mps: a cell of MPS tensors

The MPS bond dimensions of rbm2mps2.m will be smaller or equal to that of rbm2mps.m.

## Auxillary tensor programs ##
* MPS\_Canonicalize.m: Canonicalize a finite MPS and return the entanglement spectrum of each bond.

* tensor\_product.m: Contract two tensors and permute the indices to the given order. 
`C = tensor_product('acm',A,'abm',X,'cb')` does C(a,c,m) = sum_{c}A_{a,b,m)X_{c,b}

* tensor\_reshape.m: 
 `B = tensor\_reshape( A,'abcd','ac','bd')` transforms a tensor A into a matrix B( ac,bd ) =  A(a,b,c,d).


## Test program ##
* Test.m: Using the RBM architecure in Fig.1(a) as an example, we construct the MPS with two approaches. The bond dimensions are consistent with Fig. 2(c) and Fig. 4(c). The two MPS are identical in their canonical form, which also demostrates simplifying RBM using tensor techniques discussed in Sec.V of the paper. 
