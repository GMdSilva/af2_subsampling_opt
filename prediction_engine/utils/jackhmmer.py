"""
Python wrapper for jackhmmer,
taken from the original AlphaFold distribution
"""

import os
import pickle
from concurrent import futures
from sys import version_info
from urllib import request

import tqdm.notebook
from alphafold.data import parsers
from alphafold.data.tools import jackhmmer

python_version = f"{version_info.major}.{version_info.minor}"
os.environ['TF_FORCE_UNIFIED_MEMORY'] = '1'
os.environ['XLA_PYTHON_CLIENT_MEM_FRACTION'] = '2.0'

GIT_REPO = 'https://github.com/deepmind/alphafold'
SOURCE_URL = 'https://storage.googleapis.com/alphafold/alphafold_params_2021-07-14.tar'
TQDM_BAR_FORMAT = '{l_bar}{bar}| {n_fmt}/{total_fmt} [elapsed: {elapsed} remaining: {remaining}]'

os.system("mkdir -m 777 --parents /tmp/ramdisk")
os.system("mount -t tmpfs -o size=9G ramdisk /tmp/ramdisk")

def run_hmmer(prefix, sequence):
    pickled_msa_path = f"{prefix}.jackhmmer.pickle"
    if os.path.isfile(pickled_msa_path):
        msas_dict = pickle.load(open(pickled_msa_path, "rb"))
        msas, deletion_matrices = (msas_dict[k] for k in ['msas', 'deletion_matrices'])
        full_msa = []
        for msa in msas:
            full_msa += msa
    else:
        # --- Find the closest source ---
        test_url_pattern = 'https://storage.googleapis.com/alphafold-colab{:s}/latest/uniref90_2021_03.fasta.1'
        ex = futures.ThreadPoolExecutor(3)

        def fetch(source):
            request.urlretrieve(test_url_pattern.format(source))
            print(source)
            return source

        fs = [ex.submit(fetch, source) for source in ['', '-europe', '-asia']]
        source = None
        for f in futures.as_completed(fs):
            source = f.result()
            print(source)
            ex.shutdown()
            break

        jackhmmer_binary_path = 'jackhmmer'
        dbs = []

        num_jackhmmer_chunks = {'uniref90': 59, 'smallbfd': 17, 'mgnify': 71}
        total_jackhmmer_chunks = sum(num_jackhmmer_chunks.values())

        with tqdm.notebook.tqdm(total=total_jackhmmer_chunks, bar_format=TQDM_BAR_FORMAT) as pbar:
            def jackhmmer_chunk_callback(i):
                pbar.update(n=1)

            print(f'https://storage.googleapis.com/alphafold-colab{source}/latest/uniref90_2021_03.fasta')

            pbar.set_description('Searching uniref90')
            jackhmmer_uniref90_runner = jackhmmer.Jackhmmer(
                binary_path=jackhmmer_binary_path,
                database_path=f'https://storage.googleapis.com/alphafold-colab{source}/latest/uniref90_2021_03.fasta',
                get_tblout=True,
                num_streamed_chunks=num_jackhmmer_chunks['uniref90'],
                streaming_callback=jackhmmer_chunk_callback,
                z_value=135301051)
            dbs.append(('uniref90', jackhmmer_uniref90_runner.query('target.fasta')))

            pbar.set_description('Searching smallbfd')
            jackhmmer_smallbfd_runner = jackhmmer.Jackhmmer(
                binary_path=jackhmmer_binary_path,
                database_path=f'https://storage.googleapis.com/alphafold-colab{source}/latest/bfd-first_non_consensus_sequences.fasta',
                get_tblout=True,
                num_streamed_chunks=num_jackhmmer_chunks['smallbfd'],
                streaming_callback=jackhmmer_chunk_callback,
                z_value=65984053)
            dbs.append(('smallbfd', jackhmmer_smallbfd_runner.query('target.fasta')))

            pbar.set_description('Searching mgnify')
            jackhmmer_mgnify_runner = jackhmmer.Jackhmmer(
                binary_path=jackhmmer_binary_path,
                database_path=f'https://storage.googleapis.com/alphafold-colab{source}/latest/mgy_clusters_2019_05.fasta',
                get_tblout=True,
                num_streamed_chunks=num_jackhmmer_chunks['mgnify'],
                streaming_callback=jackhmmer_chunk_callback,
                z_value=304820129)
            dbs.append(('mgnify', jackhmmer_mgnify_runner.query('target.fasta')))

        # --- Extract the MSAs and visualize ---
        # Extract the MSAs from the Stockholm files.
        # NB: deduplication happens later in pipeline.make_msa_features.

        mgnify_max_hits = 5010000
        msas = []
        deletion_matrices = []
        for db_name, db_results in dbs:
            unsorted_results = []
            for i, result in enumerate(db_results):
                msa, deletion_matrix, target_names = parsers.parse_stockholm(result['sto'])
                e_values_dict = parsers.parse_e_values_from_tblout(result['tbl'])
                e_values = [e_values_dict[t.split('/')[0]] for t in target_names]
                zipped_results = zip(msa, deletion_matrix, target_names, e_values)
                if i != 0:
                    # Only take query from the first chunk
                    zipped_results = [x for x in zipped_results if x[2] != 'query']
                unsorted_results.extend(zipped_results)
            sorted_by_evalue = sorted(unsorted_results, key=lambda x: x[3])
            db_msas, db_deletion_matrices, _, _ = zip(*sorted_by_evalue)
            if db_msas:
                if db_name == 'mgnify':
                    db_msas = db_msas[:mgnify_max_hits]
                    db_deletion_matrices = db_deletion_matrices[:mgnify_max_hits]
                msas.append(db_msas)
                deletion_matrices.append(db_deletion_matrices)
                msa_size = len(set(db_msas))
                print(f'{msa_size} Sequences Found in {db_name}')

            pickle.dump({"msas": msas, "deletion_matrices": deletion_matrices},
                        open(pickled_msa_path, "wb"))

    return msas