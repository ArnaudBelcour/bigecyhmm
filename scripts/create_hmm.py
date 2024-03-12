import argparse
import os

try:
    import pyhmmer
except:
    print('This script requires pyhmmer, try to install it with "pip install pyhmmer"')

try:
    from pyfamsa import Aligner, Sequence
except:
    print('This script requires pyfamsa, try to install it with "pip install pyfamsa"')

try:
    from pytrimal import Alignment, AutomaticTrimmer
except:
    print('This script requires pytrimal, try to install it with "pip install pytrimal"')

try:
    from Bio import SeqIO
except:
    print('This script requires biopython, try to install it with "pip install biopython"')

parser = argparse.ArgumentParser()
parser.add_argument(
    '-i',
    '--input',
    dest='input',
    required=True,
    help='Input protein fasta files to use to create a HMM profile.',
    metavar='INPUT_FILE')

args = parser.parse_args()
input_fasta_file = args.input
input_basename = os.path.splitext(os.path.basename(input_fasta_file))[0]
output_file = input_basename + '.hmm'

sequences = [Sequence(r.id.encode(), str(r.seq).encode()) for r in SeqIO.parse(input_fasta_file, "fasta")]

aligner = Aligner(guide_tree="upgma")
msa = aligner.align(sequences)
sequence_names = [sequences.id for sequences in msa]
sequences = [sequences.sequence for sequences in msa]

alignment = Alignment(names=sequence_names, sequences=sequences)
trimmer = AutomaticTrimmer(method="strictplus")

trimmed = trimmer.trim(alignment)

trimmed_msa = trimmed.to_pyhmmer()
trimmed_msa.name = input_basename.encode()

alphabet = pyhmmer.easel.Alphabet.amino()
builder = pyhmmer.plan7.Builder(alphabet)
background = pyhmmer.plan7.Background(alphabet)

hmm_trimmed, _, _ = builder.build_msa(trimmed_msa.digitize(alphabet), background)

with open(output_file, "wb") as output_file:
    hmm_trimmed.write(output_file)
