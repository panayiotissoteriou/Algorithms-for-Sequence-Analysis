#!/usr/bin/python3
# python hmm.py viterbi input/sequence_X1.fasta input/A1_prior.tsv input/E1_prior.tsv
#python hmm.py -vv viterbi input/sequence_X1.fasta input/A1_prior.tsv input/E1_prior.tsv
#python hmm.py -vv forward input/sequence_X1.fasta input/A1_prior.tsv input/E1_prior.tsv
# python hmm.py -vv backward input/sequence_X1.fasta input/A1_prior.tsv input/E1_prior.tsv
#python hmm.py -h
# for 1 iteration
# python hmm.py -vv baumwelch -i 1 input/sequence_X1.fasta input/A1_prior.tsv input/E1_prior.tsv
# python hmm.py -vv baumwelch input/sequence_X1.fasta input/A2_prior.tsv input/E2_prior.tsv

"""
DESCRIPTION:
    Template code for the Hidden Markov Models assignment in the Algorithms in Sequence Analysis course at the VU.

INSTRUCTIONS:
    Complete the code (compatible with Python 3!) upload to CodeGrade via corresponding Canvas assignment.

AUTHOR:
    <Panayiotis Soteriou 2738068>
"""

import os.path as op

from os import makedirs
from math import log10
from hmm_utility import parse_args, load_fasta, load_tsv, print_trellis, print_params, serialize



def viterbi(X,A,E):
    """Given a single sequence, with Transition and Emission probabilities,
    return the most probable state path, the corresponding P(X), and trellis."""

    allStates = A.keys() # states (D = domain, L = linked, B = beginnign, E = end)
    emittingStates = E.keys()   # emitting states (D and L)
    L = len(X) + 2

    # Initialize
    V = {k:[0] * L for k in allStates} # The Viterbi trellis
    V['B'][0] = 1.

    # Middle columns
    for i,s in enumerate(X):        # i = index, s = value/AA; loop goes over the row
        for l in emittingStates:        # goes down the column / 1 value per row (emitting states ONLY); l = L or D
            terms = [V[k][i] * A[k][l] for k in allStates]      # V[k][i] (= cell in previous column) * A[k][i] (= transition prob), all values inserted in a list
            V[l][i+1] = max(terms) * E[l][s]        # k = B, L, D, E; selects the highest value and adds it in the cell of interest
                                                    # E[l][s] = emission probability of L or D (l) at AA (s)

    # Last column
    for k in allStates:
        term = V[k][i+1] * A[k]['E'] # calculates the probabilities using the termination formula
        if term > V['E'][-1]:   # it chooses the largest value of the 2nd to last column
            V['E'][-1] = term
            pi = k # Last state of the State Path   # path (pi = E)

    # FOR VITERBI ONLY: Trace back the State Path
    l = pi      # 'E' is the first letter, and the
    i = L-2
    while i:
        i -= 1
        for k in emittingStates:
            if V[k][i] * A[k][l] * E[l][X[i]] == V[l][i+1]:
                pi = k + pi         # you add string on the pi string (the path)
                l = k
                break

    P = V['E'][-1] # The Viterbi probability: P(X,pi|A,E)
    return(pi,P,V) # Return the state path, Viterbi probability, and Viterbi trellis



def forward(X,A,E):
    """Given a single sequence, with Transition and Emission probabilities,
    return the Forward probability and corresponding trellis."""

    allStates = A.keys()
    emittingStates = E.keys()
    L = len(X) + 2

    # Initialize
    F = {k:[0] * L for k in allStates}
    F['B'][0] = 1

    #####################
    # START CODING HERE #
    #####################
    # HINT: The Viterbi and Forward algorithm are very similar! 
    # Adapt the viterbi() function to account for the differences.

    # Middle columns
    # for ...
    for i,s in enumerate(X):        # i = index, s = value/AA; loop goes over the row
        for l in emittingStates:        # goes down the column / 1 value per row (emitting states ONLY); l = L or D
            terms = [F[k][i] * A[k][l] for k in allStates]      # F[k][i] (= cell in previous column) * A[k][i] (= transition prob), all values inserted in a list
            F[l][i+1] = sum(terms) * E[l][s]


    # Last columns
    for k in allStates:
        F['E'][-1] += F[k][i+1] * A[k]['E'] # calculates/sums the probabilities using the termination formula

    #####################
    #  END CODING HERE  #
    #####################

    P = F['E'][-1] # The Forward probability: P(X|A,E)
    return(P,F)



def backward(X,A,E):
    """Given a single sequence, with Transition and Emission probabilities,
    return the Backward probability and corresponding trellis."""

    allStates = A.keys()
    emittingStates = E.keys()
    L = len(X) + 2

    # Initialize
    B = {k:[0] * L for k in allStates} # The Backward trellis
    for k in allStates:
        B[k][-2] = A[k]['E']    #

    #####################
    # START CODING HERE #
    #####################s
    # Remaining columns
    #print(len(X))
    for i in range(L-3,-1,-1): # up to and not including -1, which means that the loop stops at 0
        s = X[i] # value
        for k in allStates:
            terms = [B[l][i+1] * E[l][s] * A[k][l] for l in emittingStates]
            B[k][i] = sum(terms)

    #     ...

    #####################
    #  END CODING HERE  #
    #####################

    P = B['B'][0] # The Backward probability -- should be identical to Forward!
    return(P,B)

# python hmm.py -vv baumwelch input/sequence_X1.fasta input/A1_prior.tsv input/E1_prior.tsv
# python hmm.py -vv baumwelch input/training_sequences.fasta input/A1_prior.tsv input/E1_prior.tsv
# python hmm.py -vv baumwelch input/training_sequences.fasta input/A2_prior.tsv input/E2_prior.tsv
# python hmm.py -vv baumwelch input/training_sequences.fasta input/A2_prior.tsv input/E3_prior.tsv
# python hmm.py -vv -i 150 baumwelch input/training_sequences.fasta input/A1_prior.tsv input/E2_prior.tsv
# python hmm.py -vv baumwelch input/training_sequences.fasta input/A2_prior_2.tsv input/E2_prior.tsv
# python hmm.py -vv baumwelch input/training_sequences.fasta input/A2_prior_2.tsv input/E3_prior.tsv
# python hmm.py -vv -i 150 baumwelch input/training_sequences.fasta input/A2_prior_2.tsv input/E4_prior.tsv


def baumwelch(set_X,A,E):
    """Given a set of sequences X and priors A and E,
    return the Sum Log Likelihood of X given the priors,
    along with the calculated posteriors for A and E."""

    allStates = A.keys()
    emittingStates = E.keys()
    
    # Initialize a new (posterior) Transition and Emission matrix
    new_A = {}
    for k in A:
        new_A[k] = {l:0 for l in A[k]}
    
    new_E = {}
    for k in E:
        new_E[k] = {s:0 for s in E[k]}

    # Iterate through all sequences in X
    SLL = 0 # Sum Log-Likelihood
    for X in set_X:
        P,F = forward(X,A,E)  # Save both the forward probability and the forward trellis
        _,B = backward(X,A,E) # Forward P == Backward P, so only save the backward trellis
        SLL += log10(P)

        #####################
        # START CODING HERE #
        #####################
        ## ATTEMPT FOR MATRIX [A]
        #print(len(X))
        for k in allStates:  # up to and not including -1, which means that the loop stops at 0
            #print(k)
            for l in allStates:
                    sum_akl = 0
                    if l == 'B':
                        sum_akl += 0

                    elif l in emittingStates:
                        sum_akl += sum([F[k][i] * A[k][l] * E[l][X[i]] * B[l][i + 1] for i in range(len(X))]) / P  # for each k

                    else:# l == 'E':
                        # sum_akl += F[k][len(X)] * A[k][l] /P
                        sum_akl += F[k][len(X)] * A[k][l] / P

                    new_A[k][l] += sum_akl





        #"""E MATRIX"""
        ### E MATRIX ###

        for k in emittingStates:
            for i in range(len(X)):
                new_E[k][X[i]] += (F[k][i+1] * B[k][i+1]) / P


    for k in allStates:
        sum_overall = 0
        for l in allStates:
            sum_overall += new_A[k][l]
        for i in allStates:
            if sum_overall != 0:
                new_A[k][i] = new_A[k][i] / sum_overall

    for k in emittingStates:
        sum_E = 0
        for l in new_E[k]:
            sum_E += new_E[k][l]
        for i in new_E[k]:
            new_E[k][i] = new_E[k][i] / sum_E


        # Inside the for loop: Expectation
        # Calculate the expected transitions and emissions for the sequence.
        # Add the contributions to your posterior matrices.
        # Remember to normalize to the sequence's probability P!
        
    # Outside the for loop: Maximization
    # Normalize row sums to 1 (except for one row in the Transition matrix!)
    # new_A = ...
    # new_E = ...

    #####################
    #  END CODING HERE  #
    #####################

    return(SLL,new_A,new_E)



def main(args = False):
    "Perform the specified algorithm, for a given set of sequences and parameters."
    
    # Process arguments and load specified files
    if not args: args = parse_args()

    cmd = args.command            # viterbi, forward, backward or baumwelch
    verbosity = args.verbosity
    set_X, labels = load_fasta(args.fasta)  # List of sequences, list of labels
    A = load_tsv(args.transition) # Nested Q -> Q dictionary
    E = load_tsv(args.emission)   # Nested Q -> S dictionary
    
    def save(filename, contents):
        if args.out_dir:
            makedirs(args.out_dir, exist_ok=True) # Make sure the output directory exists.
            path = op.join(args.out_dir,filename)
            with open(path,'w') as f: f.write(contents)
        # Note this function does nothing if no out_dir is specified!



    # VITERBI
    if cmd == 'viterbi':
        for j,X in enumerate(set_X): # For every sequence:
            # Calculate the most probable state path, with the corresponding probability and matrix
            Q, P, T = viterbi(X,A,E)

            # Save and/or print relevant output
            label = labels[j]
            save('%s.path' % label, Q)
            save('%s.matrix' % label, serialize(T,X))
            save('%s.p' % label, '%1.2e' % P)
            print('>%s\n Path = %s' % (label,Q))
            if verbosity: print(' Seq  = %s\n P    = %1.2e\n' % (X,P))
            if verbosity >= 2: print_trellis(T, X)
            


    # FORWARD or BACKWARD
    elif cmd in ['forward','backward']:
        if cmd == 'forward':
            algorithm = forward
        elif cmd == 'backward':
            algorithm = backward

        for j,X in enumerate(set_X): # For every sequence:
            # Calculate the Forward/Backward probability and corresponding matrix
            P, T = algorithm(X,A,E)

            # Save and/or print relevant output
            label = labels[j]
            save('%s.matrix' % label, serialize(T,X))
            save('%s.p' % label, '%1.2e' % P)
            if verbosity >= 2:
                print('\n>%s\n P = %1.2e\n' % (label,P))
                print_trellis(T, X)
            elif verbosity: print('>%-10s\tP = %1.2e' % (label,P))



    # BAUM-WELCH TRAINING
    elif cmd == 'baumwelch':
        # Initialize
        i = 1
        i_max = args.max_iter
        threshold = args.conv_thresh

        current_SLL, A, E = baumwelch(set_X,A,E)
        if verbosity: print('Iteration %i, prior SLL = %1.2e' % (i,current_SLL))
        if verbosity >= 2: print_params(A,E)
        
        last_SLL = current_SLL - threshold - 1 # Iterate at least once

        # Iterate until convergence or limit
        while i < i_max and current_SLL - last_SLL > threshold:
            i += 1
            last_SLL = current_SLL

            # Calculate the Sum Log-Likelihood of X given A and E,
            # and update the estimates (posteriors) for A and E.
            current_SLL, A, E = baumwelch(set_X,A,E)

            if verbosity: print('Iteration %i, prior SLL = %1.2e' % (i,current_SLL))
            if verbosity >= 2: print_params(A,E)

        converged = current_SLL - last_SLL <= threshold
        try:
            final_SLL = sum([log10(forward(X,A,E)[0]) for X in set_X])
        except ValueError:
            final_SLL = 0

        # Save and/or print relevant output
        save('SLL','%1.2e\t%i\t%s' % (final_SLL, i, converged))
        save('posterior_A',serialize(A))
        save('posterior_E',serialize(E))
        if verbosity: print('========================================\n')

        if converged:
            print('Converged after %i iterations.' % i)
        else:
            print('Failed to converge after %i iterations.' % i_max)

        if verbosity:
            print('Final SLL: %1.2e' % final_SLL)
            print('Final parameters:')
            print_params(A,E)



if __name__ == '__main__':
	main()