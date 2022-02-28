#!/usr/bin/python3
# python sequence_generator.py ../input/A1_prior.tsv ../input/E1_prior.tsv -n 10 -o random_seq.fasta

"""
DESCRIPTION:
    Template code for the FIRST Advanced Question of the Hidden Markov Models
    assignment in the Algorithms in Sequence Analysis course at the VU.

INSTRUCTIONS:
    Complete the code (compatible with Python 3!) upload to CodeGrade via
    corresponding Canvas assignment. Note this script will be graded manually,
    if and only if your "hmm.py" script succesfully implements Baum-Welch
    training! Continuous Feedback will not be available for this script.

AUTHOR:
    <Panayiotis Soteriou 2738068>
"""

from argparse import ArgumentParser, RawTextHelpFormatter
from hmm_utility import load_tsv
from numpy.random import choice



def parse_args():
    #####################
    # START CODING HERE #
    #####################
    # Implement a simple argument parser (WITH help documentation!) that parses
    # the information needed by main() from commandline. Take a look at the
    # argparse documentation, the parser in hmm_utility.py or align.py
    # (from the Dynamic Programming exercise) for hints on how to do this.

    parser = ArgumentParser()
    parser.add_argument('transition', help='path to a TSV formatted transition matrix')
    parser.add_argument('emission', help='path to a TSV formatted emission matrix')
    parser.add_argument('-n', dest='number', help="number of sequences to generate", type=int, default=1)
    parser.add_argument('-o', dest='out_file', help='filename to save the output', default='random_sequences.fasta')

    args = parser.parse_args()
    return args
    
    #####################
    #  END CODING HERE  #
    #####################


def generate_sequence(A,E):
    #####################
    # START CODING HERE #
    #####################
    # Implement a function that generates a random sequence using the choice()
    # function, given a Transition and Emission matrix.
    
    # Look up its documentation online:
    # https://docs.scipy.org/doc/numpy-1.15.0/reference/generated/numpy.random.choice.html
            
    letters = []
    for letter in list(E.values())[0].keys():
        letters.append(letter)

    states = []
    for state in A.keys():
        states.append(state)

    state = 'B'
    probability_a = []
    for prob in A[state].values():
        probability_a.append(prob)

    sequence = []


    state = choice(states, p=probability_a)
    while True:
        if state == 'E':
            break
        else:
            probability_e = list(E[state].values())
            sequence.append(choice(letters, p=probability_e))
            probability_a = list(A[state].values())
            state = choice(states, p=probability_a)
    
    #####################
    #  END CODING HERE  #
    #####################
    
    return sequence



def main():
    args = parse_args()
    #####################
    # START CODING HERE #
    #####################
    # Uncomment and complete (i.e. replace '?' in) the lines below:

    N = args.number  # The number of sequences to generate
    A = load_tsv(args.transition)  # Transition matrix
    E = load_tsv(args.emission)  # Emission matrix
    out_file = args.out_file
    sequences = []
    for i in range(N):
        sequences.append(generate_sequence(A, E))

    if args.out_file:
        with open(out_file, 'w') as f:
            for i in range(N):
                seq = sequences[i]
                f.write('>random_sequence_%i\n%s\n' % (i, ''.join(seq)))
        
    #####################
    #  END CODING HERE  #
    #####################
    


if __name__ == "__main__":
    main()
