#!/usr/bin/python3

"""
DESCRIPTION:
    Template code for the Dynamic Programming assignment in the Algorithms in Sequence Analysis course at the VU.
    
INSTRUCTIONS:
    Complete the code (compatible with Python 3!) upload to CodeGrade via corresponding Canvas assignment.

AUTHOR:
    <Panayiotis Soteriou 2738068>
"""



import argparse
import pickle



def parse_args():
    "Parses inputs from commandline and returns them as a Namespace object."

    parser = argparse.ArgumentParser(prog = 'python3 align.py',
        formatter_class = argparse.RawTextHelpFormatter, description =
        '  Aligns the first two sequences in a specified FASTA\n'
        '  file with a chosen strategy and parameters.\n'
        '\n'
        'defaults:\n'
        '  strategy = global\n'
        '  substitution matrix = pam250\n'
        '  gap penalty = 2')
        
    parser.add_argument('fasta', help='path to a FASTA formatted input file')
    parser.add_argument('output', nargs='*', 
        help='path to an output file where the alignment is saved\n'
             '  (if a second output file is given,\n'
             '   save the score matrix in there)')
    parser.add_argument('-v', '--verbose', dest='verbose', action='store_true',
        help='print the score matrix and alignment on screen', default=False)
    parser.add_argument('-s', '--strategy', dest='strategy',
        choices=['global','semiglobal','local'], default="global")
    parser.add_argument('-m', '--matrix', dest='substitution_matrix',
        choices=['pam250','blosum62','identity'], default='pam250')
    parser.add_argument('-g', '--gap_penalty', dest='gap_penalty', type=int,
        help='must be a positive integer', default=2)

    args = parser.parse_args()

    args.align_out = args.output[0] if args.output else False
    args.matrix_out = args.output[1] if len(args.output) >= 2 else False
                      # Fancy inline if-else statements. Use cautiously!
                      
    if args.gap_penalty <= 0:
        parser.error('gap penalty must be a positive integer')

    return args



def load_substitution_matrix(name):
    "Loads and returns the specified substitution matrix from a pickle (.pkl) file."
    # Substitution matrices have been prepared as nested dictionaries:
    # the score of substituting A for Z can be found with subst['A']['Z']
    # NOTE: Only works if working directory contains the correct folder and file!
    
    with open('substitution_matrices/%s.pkl' % name, 'rb') as f:
        subst = pickle.load(f)
    return subst
    

def load_sequences(filepath):
    "Reads a FASTA file and returns the first two sequences it contains."
    
    seq1 = []
    seq2 = []
    with open(filepath,'r') as f:
        for line in f:
            if line.startswith('>'):
                if not seq1:
                    current_seq = seq1
                elif not seq2:
                    current_seq = seq2
                else:
                    break # Stop if a 3rd sequence is encountered
            else:
                current_seq.append(line.strip())
    
    if not seq2:
        raise Exception('Error: Not enough sequences in specified FASTA file.')
    
    seq1 = ''.join(seq1)
    seq2 = ''.join(seq2)
    return seq1, seq2



def align(seq1, seq2, strategy, substitution_matrix, gap_penalty):
    "Do pairwise alignment using the specified strategy and parameters."
    # This function consists of 3 parts:
    #
    #   1) Initialize a score matrix as a "list of lists" of the appropriate length.
    #      Fill in the correct values for the first row and column given the strategy.
    #        (local / semiglobal = 0  --  global = stacking gap penalties)
    #   2) Fill in the rest of the score matrix using Dynamic Programming, accounting
    #      for the selected alignment strategy, substitution matrix and gap penalty.
    #   3) Perform the correct traceback routine on your filled in score matrix.
    #
    # Both the resulting alignment (sequences with gaps and the corresponding score)
    # and the filled in score matrix are returned as outputs.
    #
    # NOTE: You are strongly encouraged to think about how you can reuse (parts of)
    #       your code between steps 2 and 3 for the different strategies!
    
    
    ### 1: Initialize
    M = len(seq1)+1         #seq1 = x-axis/vertical axis
    N = len(seq2)+1         #seq2 = y-axis/horizontal axis
    score_matrix = []
    for i in range(M):
        row = []
        score_matrix.append(row)
        for j in range(N):
            row.append(0)

    if strategy == 'global':
        #####################
        # START CODING HERE #
        for i in range(N):
            score_matrix[0][i] = -gap_penalty*i #the gap penalty is added in the top row

        for i in range(M):
            score_matrix[i][0] = -gap_penalty*i #gap penalty added in leftmost column

        #####################
        #####################
        #  END CODING HERE  #
        #####################


    
    ### 2: Fill in Score Matrix
 
    #####################
    # START CODING HERE #
    #####################
    def dp_function(letter_seq1, letter_seq2, diagonal, horizontal, vertical):
        values_list = []

        gap_horizontal_value = horizontal - gap_penalty
        values_list.append(int(gap_horizontal_value))

        diagonal_value = diagonal + substitution_matrix[letter_seq1][letter_seq2] #match/mismatch
        values_list.append(int(diagonal_value))

        gap_vertical_value = vertical - gap_penalty
        values_list.append(int(gap_vertical_value))

        # LOCAL ALIGNMENT, NO NEGATIVE NUMBERS ALLOWED
        if strategy == 'local':
            for i in range(len(values_list)):
                if values_list[i] < 0:
                    values_list[i] = 0

        return max(values_list)



    for i in range(1,M):            #seq1 = x-axis
        for j in range(1,N):       #seq2 = y-axis
            score_matrix[i][j] = dp_function(seq1[i-1], seq2[j-1], score_matrix[i-1][j-1], score_matrix[i][j-1], score_matrix[i-1][j]) #adds the max value to the cell


    #####################
    #  END CODING HERE  #
    #####################   
    

    ### 3: Traceback
    
    #####################
    # START CODING HERE #
    #####################

    # Score of alignment
    align_score = 0

    if strategy == 'global':
        align_score = score_matrix[-1][-1]        # with gaps inserted at the appropriate positions.

    #getting alignment score for LOCAL
    elif strategy == 'local':
        max_score = 0
        count = -1
        x_axis_index = 0
        y_axis_index = 0

        for scores in score_matrix:
            count += 1
            y_axis_index_current = 0
            for number in scores:
                if number > max_score:          #getting a high score and the coordinates of that score
                    max_score = number
                    x_axis_index = count
                    y_axis_index = y_axis_index_current
                    align_score = max_score

                elif max_score == number:
                    if y_axis_index < y_axis_index_current:     #prioritising the high road (higher y axis index)
                        max_score = number
                        x_axis_index = count
                        y_axis_index = y_axis_index_current
                        align_score = max_score

                    elif y_axis_index == y_axis_index_current:      #prioritising the high road (lower x axis index)
                        if x_axis_index > count:
                            x_axis_index = count
                y_axis_index_current += 1


    #SEMIGLOBAL
    elif strategy == 'semiglobal':
        max_score_column = 0
        max_score_row = 0
        max_score = 0
        count = -1
        count_2 = -1
        highest_value_x_axis_index_row = 0
        highest_value_y_axis_index_row = 0
        highest_value_x_axis_index_column = 0
        highest_value_y_axis_index_column = 0
        highest_value_x_axis_index = 0
        highest_value_y_axis_index = 0

        for row in score_matrix: # gets the high value for rightmost column
            count_2 += 1
            if row[-1] > max_score_column:
                highest_value_x_axis_index_column = count_2
                highest_value_y_axis_index_column = len(row)-1
                max_score_column = row[-1]

        for score in score_matrix[-1]: # gets the high value for bottom row
            count += 1
            if score >= max_score_row:
                highest_value_x_axis_index_row = len(score_matrix)-1
                highest_value_y_axis_index_row = count
                max_score_row = score

        if max_score_row > max_score_column:        #prioritising the high road
            max_score = max_score_row
            highest_value_x_axis_index = highest_value_x_axis_index_row
            highest_value_y_axis_index = highest_value_y_axis_index_row
            align_score = max_score
        else:
            max_score = max_score_column
            highest_value_x_axis_index = highest_value_x_axis_index_column
            highest_value_y_axis_index = highest_value_y_axis_index_column
            align_score = max_score


    aligned_seq1 = ''  # creating the aligned_seq1 variable
    aligned_seq2 = ''  # creating the aligned_seq2 variable

    #Global
    if strategy == 'global':
        i = M
        j = N
        while i > 1 and j > 1:  # goes through the matrix starting from the lowest rightmost cell
            score_current = score_matrix[i-1][j-1]
            score_diagonal = score_matrix[i-2][j-2]
            score_up = score_matrix[i-2][j-1]
            score_left = score_matrix[i-1][j-2]

            # checking from which direction the score on the current cell was obtained and creating alignment
            if score_current == score_up - gap_penalty:
                aligned_seq1 += seq1[i-2]
                aligned_seq2 += '-'
                i = i-1

            elif score_current == score_diagonal + substitution_matrix[seq1[i-2]][seq2[j-2]]:
                aligned_seq1 += seq1[i-2]
                aligned_seq2 += seq2[j-2]
                i = i-1
                j = j-1

            elif score_current == score_left - gap_penalty:
                aligned_seq1 += '-'
                aligned_seq2 += seq2[j-2]
                j = j-1

            #adding gaps if the top row or leftmost column is reached until the top leftmost cell is reached
            if i == 1 and j != 1:
                while j > i:
                    aligned_seq1 += '-'
                    aligned_seq2 += seq2[j - 2]
                    j = j - 1

            if j ==1 and i != 1:
                while i > j:
                    aligned_seq1 += seq1[i - 2]
                    aligned_seq2 += '-'
                    i = i - 1

    #Local
    if strategy == 'local':
        i = x_axis_index
        j = y_axis_index

        #
        while i > 0 and j > 0:  # goes through the matrix starting from the cell with highest score
            score_current = score_matrix[i][j]
            score_diagonal = score_matrix[i-1][j-1]
            score_up = score_matrix[i-1][j]
            score_left = score_matrix[i][j-1]

            # checking from which direction the score on the current cell was obtained and creating alignment sequences
            if score_current == score_up - gap_penalty:
                aligned_seq1 += seq1[i-1]
                aligned_seq2 += '-'
                i = i-1

            elif score_current == score_diagonal + substitution_matrix[seq1[i-1]][seq2[j-1]]:
                aligned_seq1 += seq1[i-1]
                aligned_seq2 += seq2[j-1]
                i = i-1
                j = j-1

            elif score_current == score_left - gap_penalty:
                aligned_seq1 += '-'
                aligned_seq2 += seq2[j-1]
                j = j-1

            else:   #reaching a hard 0 stops the traceback
                i = 0
                j = 0



    #Semiglobal
    if strategy == 'semiglobal':
        i = M-1
        j = N-1


        if highest_value_x_axis_index == 0 and highest_value_y_axis_index == 0: #prioritising the high road if only gaps in the alignment
            highest_value_y_axis_index = j

        while i > 0 and j > 0: # goes through the matrix starting from the rightmost, bottom cell
            score_current = score_matrix[i][j]
            score_diagonal = score_matrix[i - 1][j - 1]
            score_up = score_matrix[i - 1][j]
            score_left = score_matrix[i][j - 1]

            if i > highest_value_x_axis_index and j == highest_value_y_axis_index:  # going up in the rightmost column until max in the rightmost column and bottom row is reached
                aligned_seq1 += seq1[i - 1]
                aligned_seq2 += '-'
                i = i - 1

            elif i == highest_value_x_axis_index and j > highest_value_y_axis_index: #going left in the bottom row until max in the rightmost column and bottom row is reached
                aligned_seq1 += '-'
                aligned_seq2 += seq2[j - 1]
                j = j - 1

            elif i <= highest_value_x_axis_index and j <= highest_value_y_axis_index: #when coordinates of max reached, all 3 directions are options
                if score_current == score_up - gap_penalty:
                    aligned_seq1 += seq1[i - 1]
                    aligned_seq2 += '-'
                    i = i - 1

                elif score_current == score_diagonal + substitution_matrix[seq1[i - 1]][seq2[j - 1]]:
                    aligned_seq1 += seq1[i - 1]
                    aligned_seq2 += seq2[j - 1]
                    i = i - 1
                    j = j - 1

                elif score_current == score_left - gap_penalty:
                    aligned_seq1 += '-'
                    aligned_seq2 += seq2[j - 1]
                    j = j - 1

            # adding gaps if the top row or leftmost column is reached until the top leftmost cell is reached
            if i == 0 and j != 0:
                while j > i:
                    aligned_seq1 += '-'
                    aligned_seq2 += seq2[j - 1]
                    j = j - 1

            if j ==0 and i != 0:
                while i > j:
                    aligned_seq1 += seq1[i - 1]
                    aligned_seq2 += '-'
                    i = i - 1

    aligned_seq1 = aligned_seq1[::-1]
    aligned_seq2 = aligned_seq2[::-1]

    alignment = (aligned_seq1, aligned_seq2, align_score)
    return (alignment, score_matrix)
    #####################
    #  END CODING HERE  #
    #####################   



def print_score_matrix(s1,s2,mat):
    "Pretty print function for a score matrix."
    
    # Prepend filler characters to seq1 and seq2
    s1 = '-' + s1
    s2 = ' -' + s2
    
    # Print them around the score matrix, in columns of 5 characters
    print(''.join(['%5s' % aa for aa in s2])) # Convert s2 to a list of length 5 strings, then join it back into a string
    for i,row in enumerate(mat):               # Iterate through the rows of your score matrix (and keep count with 'i').
        vals = ['%5i' % val for val in row]    # Convert this row's scores to a list of strings.
        vals.insert(0,'%5s' % s1[i])           # Add this row's character from s2 to the front of the list
        print(''.join(vals))                   # Join the list elements into a single string, and print the line.



def print_alignment(a):
    "Pretty print function for an alignment (and alignment score)."
    
    # Unpack the alignment tuple
    seq1 = a[0]
    seq2 = a[1]
    score = a[2]
    
    # Check which positions are identical
    match = ''
    for i in range(len(seq1)): # Remember: Aligned sequences have the same length!
        match += '|' if seq1[i] == seq2[i] else ' ' # Fancy inline if-else statement. Use cautiously!
            
    # Concatenate lines into a list, and join them together with newline characters.
    print('\n'.join([seq1,match,seq2,'','Score = %i' % score]))



def save_alignment(a,f):
    "Saves two aligned sequences and their alignment score to a file."
    with open(f,'w') as out:
        out.write(a[0] + '\n') # Aligned sequence 1
        out.write(a[1] + '\n') # Aligned sequence 2
        out.write('Score: %i' % a[2]) # Alignment score


    
def save_score_matrix(m,f):
    "Saves a score matrix to a file in tab-separated format."
    with open(f,'w') as out:
        for row in m:
            vals = [str(val) for val in row]
            out.write('\t'.join(vals)+'\n')
    


def main(args = False):
    # Process arguments and load required data
    if not args: args = parse_args()
    
    sub_mat = load_substitution_matrix(args.substitution_matrix)
    seq1, seq2 = load_sequences(args.fasta)

    # Perform specified alignment
    strat = args.strategy
    gp = args.gap_penalty
    alignment, score_matrix = align(seq1, seq2, strat, sub_mat, gp)

    # If running in "verbose" mode, print additional output
    if args.verbose:
        print_score_matrix(seq1,seq2,score_matrix)
        print('') # Insert a blank line in between
        print_alignment(alignment)
    
    # Save results
    if args.align_out: save_alignment(alignment, args.align_out)
    if args.matrix_out: save_score_matrix(score_matrix, args.matrix_out)



if __name__ == '__main__':
    main()

