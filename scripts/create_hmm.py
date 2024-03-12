import argparse
import os
import json
import sys

try:
    import pyhmmer
except:
    print('This script requires pyhmmer, try to install it with "pip install pyhmmer"')

try:
    from pyfamsa import Aligner, Sequence
    from pyfamsa import __version__ as pyfamsa_version

except:
    print('This script requires pyfamsa, try to install it with "pip install pyfamsa"')

try:
    from pytrimal import Alignment, AutomaticTrimmer
    from pytrimal import __version__ as pytrimal_version

except:
    print('This script requires pytrimal, try to install it with "pip install pytrimal"')

try:
    from Bio import SeqIO
    from Bio import __version__ as biopython_version
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
parser.add_argument(
    '-o',
    '--output',
    dest='output',
    required=True,
    help='Output folder to store result file.',
    metavar='OUTPUT_FOLDER')

args = parser.parse_args()
input_fasta_file = args.input
output_folder = args.output

input_basename = os.path.splitext(os.path.basename(input_fasta_file))[0]
output_file = os.path.join(output_folder, input_basename + '.hmm')
metadata_json_file = os.path.join(output_folder, input_basename + '_metadata.json')

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

with pyhmmer.easel.SequenceFile(input_fasta_file, digital=True) as seq_file:
    sequences = list(seq_file)

# Compute the trusted cutoffs as the lowest cutoff for the trusted proteins of the dataset.
hmmer_hits = []
for hits in pyhmmer.hmmsearch(hmm_trimmed, sequences, cpus=1):
    for hit in hits:
        hmmer_hits.append(hit.score)
trusted_cutoff = min(hmmer_hits)

metadata = {}
metadata['Python_version'] = sys.version
metadata['dependencies'] = {'pyhmmer': pyhmmer.__version__, 'pyfamsa': pyfamsa_version, 'pytrimal': pytrimal_version, 'biopython': biopython_version}
metadata['input_parameters'] = {'input': input_fasta_file, 'output': output_folder}
metadata['trusted_cutoff'] = trusted_cutoff

with open(metadata_json_file, 'w') as dumpfile:
    json.dump(metadata, dumpfile, indent=4)