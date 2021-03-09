# LinearCryptanalysis

### AES Bias Visualization
Write a program to compute the bias of each combination of inputs and outputs to/from
the AES S-Box (note that each combination should include at least one of the 8 inputs and at
least one of the 8 outputs). Plot a histogram of the biases and also include a table showing
the biases and the number of combinations with that bias. (The biases should be positive
(even – why?) integers).

### Maximum Bias Subgraph Search
Given an n X n S-Box, the number of stages in the SPN (Substitution Permutation Network)
and the plaintext/ciphertext block size, write a program to obtain a linear expression with
the maximum possible bias. The expression should be a XOR of bits of the plaintext, bits of
the ciphertext and bits of the round keys (all but the last round).
“A greedy strategy will always work.” Examine the validity of this claim.
Extra Credit: Implement Step 2 of the linear cryptanalysis attack to deduce bits of the last
round key


### Team Members  -

[Neel Aryan Gupta](https://www.cse.iitb.ac.in/~neelaryan) - 180050067  
Pulkit Agrawal - 180050081\
Tathagat Verma - 180050111\
