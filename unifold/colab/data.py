import hashlib
import os
from typing import Dict, List, Sequence, Tuple, Union, Any, Optional


from unifold.data import residue_constants, protein
from unifold.msa.utils import divide_multi_chains
from unifold.dataset import load_and_process, UnifoldDataset
from unifold.symmetry import load_and_process_symmetry

import numpy as np
from Bio import SeqIO
import re
import random
import hashlib
from pathlib import Path

def read_a3m_file(file_path: Union[str, Path, Dict[str, bytes]]) -> Tuple[List[str], List[List[str]]]:
    """
    Read an A3M file and parse its contents.
    
    The A3M format is used for storing multiple sequence alignments with optional
    header information about sequence lengths and counts.
    
    Args:
        file_path: Can be either:
            - A string or Path object pointing to an A3M file
            - A dictionary with a single key-value pair where the value is bytes
              (useful for reading from compressed or in-memory sources)
    
    Returns:
        A tuple containing:
            - first_sequences: List of sequences from the first MSA line
            - msas: List of lists, where each inner list contains formatted
              sequences with their descriptions in the format "description\\nsequence\\n"
    
    Raises:
        FileNotFoundError: If the file doesn't exist (when file_path is a string/Path)
        ValueError: If the file format is invalid or malformed
    """
    # Read file contents
    lines = _read_file_lines(file_path)
    
    # Parse header if present
    lengths, counts = _parse_header(lines)
    
    # Extract descriptions and sequences
    descriptions, sequences = _extract_descriptions_and_sequences(lines)
    
    # Parse MSAs
    first_sequences, msas = _parse_msas(descriptions, sequences, lengths, counts)
    
    return first_sequences, msas


def _read_file_lines(file_path: Union[str, Path, Dict[str, bytes]]) -> List[str]:
    """Read lines from file or bytes dictionary."""
    if isinstance(file_path, dict):
        if not file_path:
            raise ValueError("Dictionary file_path must contain at least one key-value pair")
        
        bytes_content = list(file_path.values())[0]
        text_content = bytes_content.decode("utf-8")
        lines = [
            line.strip() 
            for line in text_content.split("\n") 
            if line.strip()
        ]
    else:
        file_path = Path(file_path) if isinstance(file_path, str) else file_path
        
        try:
            with open(file_path, "r", encoding="utf-8") as file:
                lines = [line.strip() for line in file if line.strip()]
        except FileNotFoundError:
            raise FileNotFoundError(f"A3M file not found: {file_path}")
        except UnicodeDecodeError:
            raise ValueError(f"File is not valid UTF-8: {file_path}")
    
    if not lines:
        raise ValueError("A3M file is empty")
    
    return lines


def _parse_header(lines: List[str]) -> Tuple[List[int], List[int]]:
    """Parse the A3M file header if present."""
    if lines[0].startswith("#"):
        header_line = lines[0][1:].strip()  # Remove '#' and strip whitespace
        header_parts = header_line.split()
        
        if len(header_parts) == 2:
            lengths_str, counts_str = header_parts
            lengths = _parse_int_list(lengths_str)
            counts = _parse_int_list(counts_str)
            
            if len(lengths) != len(counts):
                raise ValueError(
                    f"Mismatch between lengths ({len(lengths)}) and counts ({len(counts)}) "
                    f"in header: {lines[0]}"
                )
                
        elif len(header_parts) == 1:
            lengths = _parse_int_list(header_parts[0])
            counts = [1] * len(lengths)
        else:
            raise ValueError(f"Malformed header line: {lines[0]}")
        
        # Remove header from lines
        lines.pop(0)
    else:
        # No header, assume single sequence
        if len(lines) < 2:
            raise ValueError("A3M file must contain at least description and sequence")
        
        sequence_length = len(lines[1])
        lengths = [sequence_length]
        counts = [1]
    
    return lengths, counts


def _parse_int_list(int_list_str: str) -> List[int]:
    """Parse comma-separated string of integers into list."""
    if not int_list_str:
        return []
    
    result = []
    for item in int_list_str.split(","):
        item = item.strip()
        if item:  # Skip empty strings
            try:
                result.append(int(item))
            except ValueError:
                raise ValueError(f"Invalid integer in list: '{item}'")
    
    return result


def _extract_descriptions_and_sequences(lines: List[str]) -> Tuple[List[str], List[str]]:
    """Extract descriptions and sequences from alternating lines."""
    if len(lines) % 2 != 0:
        raise ValueError(
            f"A3M file must have alternating description/sequence lines. "
            f"Found {len(lines)} lines (expected even number)."
        )
    
    descriptions = [lines[i] for i in range(0, len(lines), 2)]
    sequences = [lines[i] for i in range(1, len(lines), 2)]
    
    # Validate that we have pairs
    if len(descriptions) != len(sequences):
        raise ValueError(
            f"Mismatch between descriptions ({len(descriptions)}) "
            f"and sequences ({len(sequences)})"
        )
    
    return descriptions, sequences


def _parse_msas(
    descriptions: List[str],
    sequences: List[str],
    lengths: List[int],
    counts: List[int]
) -> Tuple[List[str], List[List[str]]]:
    """
    Parse multiple sequence alignments from descriptions and sequences.
    
    Returns:
        Tuple of (first_sequences, msas) where:
            first_sequences: List of sequences from the first MSA line
            msas: List of MSA groups, each containing formatted sequences
    """
    msas_dict: Dict[str, List[str]] = {}
    first_sequences: List[str] = []
    is_first_iteration = True
    
    for desc, seq in zip(descriptions, sequences):
        accumulated_length = 0
        
        for cluster_idx, (segment_length, segment_count) in enumerate(zip(lengths, counts)):
            # Extract sequence segment
            segment_start = accumulated_length
            segment_end = segment_start + segment_length
            segment_seq = seq[segment_start:segment_end]
            
            # Store first sequences
            if is_first_iteration:
                first_sequences.extend([segment_seq] * segment_count)
            
            # Create description for this segment
            segment_desc = f"{desc}_{cluster_idx}"
            
            # Skip empty segments
            if segment_seq != "-" * len(segment_seq):
                # Add to each count in the segment
                for count_idx in range(segment_count):
                    dict_key = f"{cluster_idx}_{count_idx}"
                    
                    if dict_key not in msas_dict:
                        msas_dict[dict_key] = []
                    
                    # Format: "description\nsequence\n"
                    formatted_entry = f"{segment_desc}\n{segment_seq}\n"
                    msas_dict[dict_key].append(formatted_entry)
            
            accumulated_length += segment_length
        
        # Only first description-sequence pair contributes to first_sequences
        is_first_iteration = False
    
    # Convert dictionary to list of MSA groups
    msas = list(msas_dict.values())
    
    return first_sequences, msas

def load_crosslinks(fname):
  links = np.loadtxt(fname,dtype=str)

  if len(links.shape) == 1:
      links = np.array([links])

  crosslinks = {}

  for i,chain1,j,chain2,fdr in links:
      i = int(i)
      j = int(j)
      fdr = float(fdr)
      if not chain1 in crosslinks:
          crosslinks[chain1] = {}
      if not chain2 in crosslinks[chain1]:
          crosslinks[chain1][chain2] = []

      crosslinks[chain1][chain2].append((i-1,j-1,fdr))

  return crosslinks


def load_fasta(fname,jobname):
  def add_hash(x,y):
      return x+"_"+hashlib.sha1(y.encode()).hexdigest()[:5]

  input_sequences = []
  descriptions = []
  for seq in SeqIO.parse(fname, 'fasta'):
      input_sequences.append(str(seq.seq))
      descriptions.append(seq.description)

  basejobname = "".join(input_sequences)
  basejobname = re.sub(r'\W+', '', basejobname)
  target_id = add_hash(jobname, basejobname)

  descriptions = ['>'+target_id+' chain '+d for d in descriptions]

  return input_sequences, descriptions



def clean_and_validate_sequence(
    input_sequence: str, min_length: int, max_length: int) -> str:
  """Checks that the input sequence is ok and returns a clean version of it."""
  # Remove all whitespaces, tabs and end lines; upper-case.
  clean_sequence = input_sequence.translate(
      str.maketrans('', '', ' \n\t')).upper()
  aatypes = set(residue_constants.restypes)  # 20 standard aatypes.
  if not set(clean_sequence).issubset(aatypes):
    raise ValueError(
        f'Input sequence contains non-amino acid letters: '
        f'{set(clean_sequence) - aatypes}. AlphaFold only supports 20 standard '
        'amino acids as inputs.')
  if len(clean_sequence) < min_length:
    raise ValueError(
        f'Input sequence is too short: {len(clean_sequence)} amino acids, '
        f'while the minimum is {min_length}')
  if len(clean_sequence) > max_length:
    raise ValueError(
        f'Input sequence is too long: {len(clean_sequence)} amino acids, while '
        f'the maximum is {max_length}. You may be able to run it with the full '
        f'Uni-Fold system depending on your resources (system memory, '
        f'GPU memory).')
  return clean_sequence


def validate_input(
    input_sequences: Sequence[str],
    symmetry_group: str,
    min_length: int,
    max_length: int,
    max_multimer_length: int) -> Tuple[Sequence[str], bool]:
  """Validates and cleans input sequences and determines which model to use."""
  sequences = []

  for input_sequence in input_sequences:
    if input_sequence.strip():
      input_sequence = clean_and_validate_sequence(
          input_sequence=input_sequence,
          min_length=min_length,
          max_length=max_length)
      sequences.append(input_sequence)
  
  if symmetry_group != 'C1':
    if symmetry_group.startswith('C') and symmetry_group[1:].isnumeric():
      print(f'Using UF-Symmetry with group {symmetry_group}. If you do not '
            f'want to use UF-Symmetry, please use `C1` and copy the AU '
            f'sequences to the count in the assembly.')
      is_multimer = (len(sequences) > 1)
      return sequences, is_multimer, symmetry_group
    else:
      raise ValueError(f"UF-Symmetry does not support symmetry group "
                       f"{symmetry_group} currently. Cyclic groups (Cx) are "
                       f"supported only.")

  elif len(sequences) == 1:
    print('Using the single-chain model.')
    return sequences, False, None

  elif len(sequences) > 1:
    total_multimer_length = sum([len(seq) for seq in sequences])
    if total_multimer_length > max_multimer_length:
      raise ValueError(f'The total length of multimer sequences is too long: '
                       f'{total_multimer_length}, while the maximum is '
                       f'{max_multimer_length}. Please use the full AlphaFold '
                       f'system for long multimers.')
    print(f'Using the multimer model with {len(sequences)} sequences.')
    return sequences, True, None

  else:
    raise ValueError('No input amino acid sequence provided, please provide at '
                     'least one sequence.')


def load_feature_for_one_target(
    config, data_folder, seed=0, is_multimer=False, use_uniprot=False, symmetry_group=None,
):
    if not is_multimer:
        uniprot_msa_dir = None
        sequence_ids = ["A"]
        if use_uniprot:
            uniprot_msa_dir = data_folder

    else:
        uniprot_msa_dir = data_folder
        sequence_ids = open(os.path.join(data_folder, "chains.txt")).readline().split()
    
    batch, _ = load_and_process(
        config=config.data,
        mode="predict",
        seed=seed,
        batch_idx=None,
        data_idx=0,
        is_distillation=False,
        crosslinks='crosslinks.pkl.gz',
        sequence_ids=sequence_ids,
        monomer_feature_dir=data_folder,
        uniprot_msa_dir=uniprot_msa_dir,
    )

    batch = UnifoldDataset.collater([batch])
    return batch
