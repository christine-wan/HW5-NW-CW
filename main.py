# Import NeedlemanWunsch class and read_fasta function
from align import read_fasta, NeedlemanWunsch

def main():
    """
    This function should
    (1) Align all species to humans and print species in order of most similar to human BRD
    (2) Print all alignment scores between each species BRD2 and human BRD2
    """
    hs_seq, hs_header = read_fasta("./data/Homo_sapiens_BRD2.fa")
    gg_seq, gg_header = read_fasta("./data/Gallus_gallus_BRD2.fa")
    mm_seq, mm_header = read_fasta("./data/Mus_musculus_BRD2.fa")
    br_seq, br_header = read_fasta("./data/Balaeniceps_rex_BRD2.fa")
    tt_seq, tt_header = read_fasta("./data/tursiops_truncatus_BRD2.fa")

    # Align all species to humans and print species in order of most similar to human BRD
    # using gap opening penalty of -10 and a gap extension penalty of -1 and BLOSUM62 matrix

    # Create an instance of NeedlemanWunsch with BLOSUM62 and the given penalties
    nw = NeedlemanWunsch('./substitution_matrices/BLOSUM62.mat', gap_open=-10, gap_extend=-1)

    # Align each species to human and store scores with corresponding sequences
    hs_gg_score, hs_gg_a, hs_gg_b = nw.align(hs_seq, gg_seq)
    hs_mm_score, hs_mm_a, hs_mm_b = nw.align(hs_seq, mm_seq)
    hs_br_score, hs_br_a, hs_br_b = nw.align(hs_seq, br_seq)
    hs_tt_score, hs_tt_a, hs_tt_b = nw.align(hs_seq, tt_seq)

    # Store all scores in a list (or dictionary) for sorting
    species_scores = [
        ('Gallus gallus', hs_gg_score),
        ('Mus musculus', hs_mm_score),
        ('Balaeniceps rex', hs_br_score),
        ('Tursiops truncatus', hs_tt_score)
    ]

    # Sort species based on alignment score (highest score first)
    species_scores.sort(key=lambda x: x[1], reverse=True)

    # Print species sorted by similarity to human BRD2
    print("Species in order of most similar to human BRD2:")
    for species, score in species_scores:
        print(f"{species}: Score = {score}")

    print("\nAlignment details between each species BRD2 and human BRD2:")
    # Print all alignment details
    print_alignment('Gallus gallus', hs_gg_a, hs_gg_b, hs_gg_score)
    print_alignment('Mus musculus', hs_mm_a, hs_mm_b, hs_mm_score)
    print_alignment('Balaeniceps rex', hs_br_a, hs_br_b, hs_br_score)
    print_alignment('Tursiops truncatus', hs_tt_a, hs_tt_b, hs_tt_score)


def print_alignment(species_name, seq_a, seq_b, score):
    """
    Helper function to print alignment results.
    """
    print(f"\nAlignment with {species_name}:")
    print(f"Sequence A: {seq_a}")
    print(f"Sequence B: {seq_b}")
    print(f"Alignment Score: {score}\n")
    

if __name__ == "__main__":
    main()
