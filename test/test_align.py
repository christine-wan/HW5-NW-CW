# Importing Dependencies
import pytest
from align import NeedlemanWunsch, read_fasta
import numpy as np

def test_nw_alignment():
    """
    Write your unit test for NW alignment
    using test_seq1.fa and test_seq2.fa by
    asserting that you have correctly filled out
    the your 3 alignment matrices.
    Use the BLOSUM62 matrix and a gap open penalty
    of -10 and a gap extension penalty of -1.
    """
    seq1, _ = read_fasta("./data/test_seq1.fa")
    seq2, _ = read_fasta("./data/test_seq2.fa")

    # Define the expected matrices
    true_align = np.array([[0., -11., -12., -13.],
                           [-11., 5., -6., -7.],
                           [-12., -6., 4., -7.],
                           [-13., -7., -1., 5.],
                           [-14., -8., -6., 4.]])

    true_gaps = np.array([[0, 1, 1, 1],
                            [1, 0, 1, 1],
                            [1, 1, 0, 1],
                            [1, 1, 0, 0],
                            [1, 1, 0, 0]])

    true_back = np.array([[None, None, None, None],
                                [None, (0, 0), (1, 1), (1, 2)],
                                [None, (1, 1), (1, 1), (2, 2)],
                                [None, (2, 1), (2, 1), (2, 2)],
                                [None, (3, 1), (3, 1), (3, 2)]], dtype=object)

    # Initialize the Needleman-Wunsch class
    nw = NeedlemanWunsch('./substitution_matrices/BLOSUM62.mat', -10, -1)
    nw.align(seq1, seq2)

    # Assert the alignment matrices are as expected
    assert np.array_equal(nw._align_matrix, true_align), "Alignment matrix does not match"
    assert np.array_equal(nw._gaps, true_gaps), "Gap matrix does not match"
    assert np.array_equal(nw._back, true_back), "Backtrace matrix does not match"
    

def test_nw_backtrace():
    """
    Write your unit test for NW backtracing
    using test_seq3.fa and test_seq4.fa by
    asserting that the backtrace is correct.
    Use the BLOSUM62 matrix. Use a gap open
    penalty of -10 and a gap extension penalty of -1.
    """
    # Read sequences
    seq3, _ = read_fasta("./data/test_seq3.fa")
    seq4, _ = read_fasta("./data/test_seq4.fa")

    # Initialize NeedlemanWunsch with BLOSUM62 matrix and gap penalties
    nw = NeedlemanWunsch("./substitution_matrices/BLOSUM62.mat", gap_open=-10, gap_extend=-1)

    # Perform alignment and backtrace
    score, aligned3, aligned4 = nw.align(seq3, seq4)

    # Expected alignment
    expected_aligned_seq3 = "MAVHQLIRRP"
    expected_aligned_seq4 = "M---QLIRHP"

    # Expected alignment score
    expected_score = 17

    # Assert that the aligned sequences match the expected result
    assert aligned3 == expected_aligned_seq3, "Aligned seq3 is incorrect."
    assert aligned4 == expected_aligned_seq4, "Aligned seq4 is incorrect."

    # Assert that the alignment score matches the expected result
    assert score == expected_score, f"Alignment score is incorrect. Expected {expected_score}, got {score}."


def test_empty_sequence():
    """Test edge case for empty sequences."""
    nw = NeedlemanWunsch("./substitution_matrices/BLOSUM62.mat", gap_open=-10, gap_extend=-1)

    # Both sequences empty should raise an error
    with pytest.raises(ValueError, match="Both input sequences are empty"):
        nw.align("", "")


