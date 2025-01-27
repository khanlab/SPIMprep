import numpy as np

def compute_pairwise_offsets_from_truth(groundtruth_file, pairs_file, output_file):
    """
    Calculate ground truth pairwise offsets based on overlapping pairs and ground truth translations.

    Parameters:
    - groundtruth_file (str): Path to the file containing ground truth translations.
    - pairs_file (str): Path to the file containing overlapping tile pairs.
    - output_file (str): Path to save the calculated pairwise offsets.
    """
    # Load ground truth translations
    groundtruth = np.loadtxt(groundtruth_file)

    # Load overlapping pairs
    pairs = np.loadtxt(pairs_file, dtype=int)

    # Compute pairwise offsets
    pairwise_offsets = []
    for pair in pairs:
        t1 = groundtruth[pair[0]]
        t2 = groundtruth[pair[1]]
        offset = t2 - t1
        pairwise_offsets.append(offset)

    # Save pairwise offsets to file
   
    np.savetxt(output_file, pairwise_offsets, fmt="%.6f")
    print(f"Pairwise offsets saved to {output_file}")

compute_pairwise_offsets_from_truth(snakemake.input.true_translations,snakemake.input.pairs, snakemake.output.pairwise_offsets)

