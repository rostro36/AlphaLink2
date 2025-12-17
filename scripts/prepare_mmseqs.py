from unifold.msa.pipeline import make_msa_features, make_sequence_features
from unifold.msa import parsers
from unifold.colab.data import load_fasta, validate_input
import pickle
import os
import gzip
import re
import hashlib
from pathlib import Path
from unifold.msa.utils import divide_multi_chains
from unifold.data.utils import compress_features
from unifold.data.protein import PDB_CHAIN_IDS
from unifold.colab.mmseqs import get_msa_and_templates
import shutil
import argparse

MIN_SINGLE_SEQUENCE_LENGTH = 6
MAX_SINGLE_SEQUENCE_LENGTH = 3000
MAX_MULTIMER_LENGTH = 3000

def make_mmseqs(fasta_path, output_dir_base, use_templates):
    jobname = "curr_job"
    input_sequences, descriptions = load_fasta(fasta_path,jobname)
    def add_hash(x,y):
        return x+"_"+hashlib.sha1(y.encode()).hexdigest()[:5]

    basejobname = "".join(input_sequences)
    basejobname = re.sub(r'\W+', '', basejobname)
    target_id = add_hash(jobname, basejobname)

    # Validate the input.
    result_dir = Path(output_dir_base)
    output_dir = os.path.join(output_dir_base, target_id)
    output_dir = output_dir_base
    sequences, is_multimer, _ = validate_input(
        input_sequences=input_sequences,
        symmetry_group='C1',
        min_length=MIN_SINGLE_SEQUENCE_LENGTH,
        max_length=MAX_SINGLE_SEQUENCE_LENGTH,
        max_multimer_length=MAX_MULTIMER_LENGTH)

    if is_multimer:
        _, _, chain_mapping, _ = divide_multi_chains(target_id, output_dir, sequences, descriptions)
    else:
        chain_mapping = { 'A' : [] }

    s = []
    for des, seq in zip(descriptions, sequences):
        s += [des, seq]

    unique_sequences = []
    [unique_sequences.append(x) for x in sequences if x not in unique_sequences]

    if len(unique_sequences)==1:
        homooligomers_num = len(sequences)
    else:
        homooligomers_num = 1

    with open(f"{jobname}.fasta", "w") as f:
        f.write("\n".join(s))
    (
        unpaired_msa,
        paired_msa,
        template_results,
    ) = get_msa_and_templates(
        target_id,
        sequences,
        result_dir=result_dir,
        msa_mode="MMseqs2",
        use_templates=use_templates,
        homooligomers_num = homooligomers_num
    )
    #_, homomers = chain_mapping.items()
    for idx, seq in enumerate(sequences):
        chain_id = PDB_CHAIN_IDS[idx]
        head = chain_id
        sequence_features = make_sequence_features(
                sequence=seq, description=f'> {jobname} seq {chain_id}', num_res=len(seq)
            )
        monomer_msa = parsers.parse_a3m(unpaired_msa[idx])
        msa_features = make_msa_features([monomer_msa])
        template_features = template_results[idx]
        feature_dict = {**sequence_features, **msa_features, **template_features}
        feature_dict = compress_features(feature_dict)
        os.makedirs(output_dir, exist_ok=True)
        output_dir = Path(output_dir)
        features_output_path = output_dir/f"{chain_id}.feature.pkl.gz"
        features_h_path = output_dir/f"{head}.feature.pkl.gz"
        pickle.dump(feature_dict,gzip.GzipFile(features_output_path, "wb"),protocol=4)
        #shutil.copy(features_output_path, features_h_path)
        
        pickle.dump(
            feature_dict,
            gzip.GzipFile(features_output_path, "wb"),
            protocol=4
            )
        if is_multimer:
            multimer_msa = parsers.parse_a3m(paired_msa[idx])
            pair_features = make_msa_features([multimer_msa])
            pair_feature_dict = compress_features(pair_features)
            uniprot_output_path = output_dir/f"{chain_id}.uniprot.pkl.gz"
            uniprot_h_path = output_dir/f"{head}.uniprot.pkl.gz"
            pickle.dump(
                pair_feature_dict,
                gzip.GzipFile(uniprot_output_path, "wb"),
                protocol=4,
            )
            #shutil.copy(uniprot_output_path, uniprot_h_path)

def parse_args():
    parser = argparse.ArgumentParser(
        description="Process a FASTA file with a max template date"
    )
    parser.add_argument(
        "--fasta_path",
        type=str,
        required=True,
        help="Path to the input FASTA file"
    )

    parser.add_argument(
        "--output_dir",
        type=str,
        required=True,
        help="Directory to write output files"
    )

    parser.add_argument(
        "--max_template_date",
        type=str,
        required=True,
        help="Maximum template date (YYYY-MM-DD)"
    )

    return parser.parse_args()

if __name__ == "__main__":
    args = parse_args()
    fasta_path = Path(args.fasta_path)
    output_dir = Path(args.output_dir)
    use_templates = str(args.max_template_date)>"1901"
    make_mmseqs(fasta_path, output_dir, use_templates)