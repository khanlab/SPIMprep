import numpy as np
from lib.utils import *

overlapping_pairs = find_overlapping_pairs(snakemake.input)
np.savetxt(snakemake.output.txt,np.array(overlapping_pairs),fmt='%d')

