from unifold.msa.templates import HhsearchHitFeaturizer
from unifold.msa.pipeline import run_msa_tool, make_msa_features, make_sequence_features
from unifold.msa import parsers
from unifold.colab.data import read_a3m_file
import pickle
import os
import gzip
from pathlib import Path
from unifold.msa.utils import get_chain_id_map
from unifold.data.utils import compress_features
import shutil
import argparse

MIN_SINGLE_SEQUENCE_LENGTH = 6
MAX_SINGLE_SEQUENCE_LENGTH = 3000
MAX_MULTIMER_LENGTH = 3000

def create_dummy_file(file_path):
    """Using touch method style (like Unix touch command)"""
    dummy_path = file_path/"dummydummydummy.cif"
    try:
        # If file doesn't exist, create it
        if not os.path.exists(dummy_path):
            # Create directories first
            os.makedirs(os.path.dirname(dummy_path), exist_ok=True)
            # Create empty file
            open(dummy_path, 'w').close()
        # If file exists, just update its modification time
        else:
            os.utime(dummy_path, None)
        return True
    except Exception as e:
        print(f"Error: {e}")
        return False

def get_empty_templates(query_sequence, parent_path):
    create_dummy_file(parent_path)
    pdb_template_hits = HhsearchHitFeaturizer(parent_path, "1900-01-01",0,None,None,None).get_templates(query_sequence, [])
    return pdb_template_hits

#like pipeline
def prepare_msa(msa_path):
    result = run_msa_tool(None, None, msa_path, "a3m", True)
    parsed = parsers.parse_a3m(result["a3m"])
    msa_features = make_msa_features([parsed])
    return msa_features


def write_separate_a3m(msa, msa_path:Path):
    """
    Write a string to a file at the specified path.
    
    Args:
        msa (str): The string content to write to the file
        msa_path (str): The file path where the string should be written
    """
    with open(msa_path, 'w', encoding='utf-8') as file:
        file.write(msa)


def process_precomputed(a3m_path:Path, parent_path:Path):
    first_sequences, msas = read_a3m_file(a3m_path)
    #TODO prepare a3m with deletion of inserts? and adding of gaps if needed?
    make_chains_txt(parent_path, first_sequences)
    #msas = msas as description\nseq\n
    chain_ids = "ABCDEFGHIJKLMNOPQRSTUVXYZ"[:len(first_sequences)]
    descriptions = [f"seq_{it}" for it, _ in enumerate(first_sequences)]
    #_, _, mapping = get_chain_id_map(first_sequences, descriptions)
    for sequence, msa, chain_id in zip(first_sequences, msas, chain_ids):
        msa_path = parent_path/f"{sequence[:10]}_single.a3m"
        #head = homomers[0] if len(homomers)>0 else chain_id
        msa = "".join(msa)
        write_separate_a3m(msa, msa_path)
        # how to concat msa_features?
        msa_features = prepare_msa(msa_path)
        sequence_features = make_sequence_features(sequence=sequence, description="first", num_res=len(sequence))
        templates_results = get_empty_templates(sequence, parent_path)
        feature_dict = {**sequence_features, **msa_features, **templates_results.features}
        feature_dict = compress_features(feature_dict)
        os.makedirs(parent_path, exist_ok=True)

        features_output_path = parent_path/f"{chain_id}.feature.pkl.gz"
        #features_h_path = parent_path/f"{head}.feature.pkl.gz"

        pickle.dump(
            feature_dict,
            gzip.GzipFile(features_output_path, "wb"),
            protocol=4
            )
        #shutil.copy(features_output_path, features_h_path)

        if len(first_sequences)>1:
            multimer_msa = parsers.parse_a3m(msa)
            pair_features = make_msa_features([multimer_msa])
            pair_feature_dict = compress_features(pair_features)
            uniprot_output_path = parent_path/f"{chain_id}.uniprot.pkl.gz"
            #uniprot_h_path = parent_path/f"{head}.uniprot.pkl.gz"
            pickle.dump(
                pair_feature_dict,
                gzip.GzipFile(uniprot_output_path, "wb"),
                protocol=4,
            )
            #shutil.copy(uniprot_output_path, uniprot_h_path)
    #return {**sequence_features, **msa_features, **templates_results.features}


def make_chains_txt(output_dir, input_sequences):
    chain_order_path = os.path.join(output_dir, "chains.txt")
    ALPHABET = "ABCDEFGHIJKLMNOPQRSTUVXYZ"
    with open(chain_order_path, "w") as f:
        f.write(" ".join([char for char in ALPHABET[:len(input_sequences)]] ))

def parse_args():
    parser = argparse.ArgumentParser(
        description="Process a precomputed a3m file"
    )
    parser.add_argument(
        "--input_path",
        type=str,
        required=True,
        help="Path to the input a3m file"
    )

    parser.add_argument(
        "--output_dir",
        type=str,
        required=True,
        help="Directory to write output files"
    )
    return parser.parse_args()  

if __name__ == "__main__":
    args = parse_args()
    input_path = Path(args.input_path)
    output_dir = Path(args.output_dir)
    process_precomputed(input_path, output_dir)