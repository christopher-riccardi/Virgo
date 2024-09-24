__version__="1.0.0"
__citation__="TBD"
"""
Copyright (c) 2024 Christopher Riccardi, Yuqiu Wang

Permission is hereby granted, free of charge, to any person obtaining
a copy of this software and associated documentation files (the
"Software"), to deal in the Software without restriction, including
without limitation the rights to use, copy, modify, merge, publish,
distribute, sublicense, and/or sell copies of the Software, and to
permit persons to whom the Software is furnished to do so, subject to
the following conditions:

The above copyright notice and this permission notice shall be
included in all copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,
EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF
MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND
NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE
LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION
OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION
WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
"""

from multiprocessing import Pool
import pandas as pd
import numpy as np
import subprocess
import argparse
import logging
import sys, os
import pickle
import shutil
import time
import glob

logging.basicConfig(format='%(asctime)s - %(funcName)s:%(lineno)d [%(levelname)s] %(message)s', datefmt='%d-%b-%y %H:%M:%S', level=logging.INFO)

def parse_arguments():
    parser = argparse.ArgumentParser(description='')
    parser.add_argument('-d', '--data', type=str, required=True, 
                        help='Virgo database folder')
    parser.add_argument('-i', '--input', type=str, required=True, 
                        help='Folder with input query fasta files')
    parser.add_argument('-o', '--output', type=str, required=True,
                        help='Name of output directory that will contain Virgo results')
    parser.add_argument('-t', '--num_threads', type=int, required=False,
                        help='Number of threads to use for multiprocessing')
    parser.add_argument('--with_replacement', action='store_true',
                        help='Skip database entries that have the same name as the queries. Useful for Leave-One-Out studies.')
    parser.add_argument('--no_gc', action='store_true',
                        help='Do not use G+C content to break ties. Default is to use it.')
    parser.add_argument('--version', action='version', version=f'{__version__}')
    args = parser.parse_args()
    return args
    
def CreateDirectory(directory_path):
    try:
        os.mkdir(directory_path)
    except FileExistsError:
        pass # Ignore if already present
    except: # Anything else produces error, return status 1
        logging.error(f'Cannot create directory {directory_path}')
        return 1
    logging.info(f'Successfully created/updated directory at {directory_path}')
    return 0

def Seq2Dict(file):
    d = {}
    header = None
    lines = [line.rstrip() for line in open(file)]
    ftells = []
    for i, line in enumerate(lines):
        if line.startswith('>'):
            header = line.split('>')[1] #.split()[0] Normally my Seq2Dict function would split this too.
            d[header] = ""
            ftells.append(i+1)
    ftells.append(i+2)
    for i, header in enumerate(list(d.keys())):
        d[header] = ''.join(lines[ftells[i]:ftells[i+1]-1])
    return d

def jaccard_similarity(set1, set2):
    intersection = len(set1 & set2)
    union = len(set1 | set2)
    return intersection / union if union != 0 else 0

def compute_similarity_matrix(list_a, list_b):
    num_a = len(list_a)
    num_b = len(list_b)
    
    similarity_matrix = np.zeros((num_a, num_b))
    for i, set_a in enumerate(list_a):
        for j, set_b in enumerate(list_b):
            similarity_matrix[i, j] = jaccard_similarity(set_a, set_b)
    
    return similarity_matrix

def best_match(index_set, similarity_matrix):
    return np.max(similarity_matrix[index_set])

def bidirectional_subsethood(list_a, list_b):
    similarity_matrix = compute_similarity_matrix(list_a, list_b)
    if np.sum(similarity_matrix) == 0: return 0
    coverage_a = np.sum([best_match(i, similarity_matrix) for i in range(len(list_a))])
    coverage_b = np.sum([best_match(j, similarity_matrix.T) for j in range(len(list_b))])
    return (coverage_a + coverage_b) / (len(list_a) + len(list_b))

def run_prodigal(input_fasta):
    output_fasta = io_map[input_fasta]
    prodigal_command = ['prodigal-gv', 
                        '-p', 'meta',
                        '-i', input_fasta,
                        '-a', output_fasta]
    result = subprocess.run(
        prodigal_command,
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
        text=True)
    ## We run prodigal 
    gc_contents = []
    weights = []
    try:
        seq = Seq2Dict(output_fasta)
    except:
        return None
    for key, value in seq.items():
        gc = float(key.split('gc_cont=')[1])
        size = len(value)
        gc_contents.append(gc)
        weights.append(size)
    bname = os.path.splitext(os.path.basename(output_fasta))[0]
    prot_counter = 1
    for key, value in seq.items():
        seq[key] = f'>{bname}_{prot_counter}\n{value}\n'
        prot_counter += 1
    if prot_counter == 1: return None
    return (''.join([value for value in seq.values()]), bname, np.average(gc_contents)) ## return the formatted protein sequences

def run_mmseqs(input_fasta):
    cmd = ['mmseqs', 'easy-search', 
           input_fasta, params['virus_specific_markers'], params['mmseqs_output'], params['tmp_dir'],
           '-s', '7.5',
           '-e', '1e-3',
           '-c', '0.2',
           '--cov-mode', '1',
           '--threads', str(params['num_threads'])]
    try:
        result = subprocess.run(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True, check=True)
    except subprocess.CalledProcessError as e:
        logging.error(f'{e}')
        return None
    logging.info("Markers alignment complete")
    return result

def m8_reader(m8_file):
    m8 = pd.read_csv(m8_file,
                    sep='\t',
                    header=None,
                    names=['query',
                        'target',
                        'percid',
                        'alnlen',
                        'mis',
                        'gaps',
                        'qstart',
                        'qend',
                        'tstart',
                        'tend',
                        'evalue',
                        'bitscore'])
    m8 = m8[['query', 'target']]
    return m8

def find_one_virus_replacement(virus_id):
    chosen, ties = "nan;nan;nan;nan;nan;nan;nan;nan;nan", []
    s = 0
    query_seq = queries[virus_id]
    use_gc = 1 - args.no_gc
    query_gc = gc_map[virus_id]
    
    for _, representation in database.items():
        if _ == virus_id:
            continue
        
        score = bidirectional_subsethood(query_seq, representation[0])
        if 0 == score: continue
        if score > s:
            s = score
            ties = [representation[1]]
            chosen = representation[1]
        elif score == s:
            ties.append(representation[1])
    if s == 0:
        return virus_id, None, 0, 0
    if use_gc and len(ties) > 1:
        chosen = min(ties, key=lambda elem: abs(query_gc - float(elem.split(';')[-1])))
    gc_diff = f"{abs(query_gc - float(chosen.split(';')[-1])):4.4f}"
    chosen = chosen.split(';')[:-3] ##
    chosen.append(gc_diff)
    chosen = ';'.join(chosen)
    return virus_id, chosen, s, ties

def find_one_virus(virus_id):
    chosen, ties = "nan;nan;nan;nan;nan;nan;nan;nan;nan", []
    s = 0
    query_seq = queries[virus_id]
    use_gc = 1 - args.no_gc
    query_gc = gc_map[virus_id]
    
    for _, representation in database.items():
        score = bidirectional_subsethood(query_seq, representation[0])
        if 0 == score: continue
        if score > s:
            s = score
            ties = [representation[1]]
            chosen = representation[1]
        elif score == s:
            ties.append(representation[1])
    if s == 0:
        return virus_id, None, 0, 0
    if use_gc and len(ties) > 1:
        chosen = min(ties, key=lambda elem: abs(query_gc - float(elem.split(';')[-1])))
    gc_diff = f"{abs(query_gc - float(chosen.split(';')[-1])):4.4f}"
    chosen = chosen.split(';')[:-3] ##
    chosen.append(gc_diff)
    chosen = ';'.join(chosen)
    return virus_id, chosen, s, ties

def get_tie_score(ties):
    """
    We wanted to create a straightforward way to determine the uniqueness of a virus family given the observed score.
    One family: tie score = 1
    Two families: tie score = 0.5
    and so on.
    """
    n = len(ties)
    if n < 2:
        return 1.0
    return 1 / len( set( [x.split(';')[5] for x in ties] ) )

def check_input_paths(paths):
    """
    Check if the given input paths (files or directories) exist.

    :param paths: List of paths (files or directories) provided by the user
    :return: Boolean (True if all paths exist, False if any path is missing)
    """
    missing_paths = []
    
    for path in paths:
        if not (os.path.isfile(path) or os.path.isdir(path)):
            missing_paths.append(path)

    if missing_paths:
        # Print the error message and stop the script
        logging.error(f"The following path(s) do not exist: {', '.join(missing_paths)}")
        sys.exit(1)  # Exit with an error code
    else:
        logging.info("All input files or directories exist.")
        return True

if __name__=='__main__':
    args = parse_arguments()
    sys.stdout.write(f'This is Virgo v{__version__}\n')

    params = {}
    params['input_dir'] = args.input
    params['output_dir'] = args.output
    params['database_dir'] = args.data
    params['with_replacement'] = args.with_replacement
    params['with_gc'] = 1 - args.no_gc
    params['virus_specific_markers'] = os.path.join(params['database_dir'], 'DB')
    params['database_file'] = os.path.join(params['database_dir'], 'database.pkl')
    params['tmp_dir'] = os.path.join(params['output_dir'], 'tmp_dir')
    params['mmseqs_output'] = os.path.join(params['tmp_dir'], 'output.m8')
    params['num_threads'] = args.num_threads
    params['merged_orfs'] = os.path.join(params['tmp_dir'], 'merged.faa')
    params['results_file'] = os.path.join(params['output_dir'], 'results.csv')

    logging.info(f"Parameters set:\n\
          input_dir: {params['input_dir']},\n\
          output_dir: {params['output_dir']},\n\
          database_dir: {params['database_dir']},\n\
          virus_specific_markers: {params['virus_specific_markers']},\n\
          database_file: {params['database_file']},\n\
          replacement: {params['with_replacement']},\n\
          with_gc: {params['with_gc']},\n\
          tmp_dir: {params['tmp_dir']},\n\
          mmseqs_output: {params['mmseqs_output']},\n\
          num_threads: {params['num_threads']},\n\
          merged_orfs: {params['merged_orfs']},\n\
          results_file: {params['results_file']}\n\
          "
    )
    
    check_input_paths([args.input, args.data])

    logging.info('[0]')
    if CreateDirectory(args.output) == 1:
        sys.exit(1)

    if CreateDirectory(params['tmp_dir']) == 1:
        sys.exit(1)

    allowed_extensions = {'.fa', '.fasta', '.fna', '.fas'}
    logging.info('Reading files from input directory')
    input_files = glob.glob(os.path.join(params['input_dir'], '*'))
    input_files = [file for file in input_files if os.path.splitext(file)[1] in allowed_extensions]
    logging.info(f'n={len(input_files)} have a suitable FASTA extension')
    if len(input_files) == 0:
        logging.error(f'No files had a suitable file extension. Allowed extension are {allowed_extensions}')
        shutil.rmtree(params['tmp_dir'])
        sys.exit(1)
    output_files = [os.path.join(params['tmp_dir'], os.path.basename(os.path.splitext(file)[0]) + '.faa') for file in input_files]
    io_map = {input_files[i]:output_files[i] for i in range(len(input_files))}

    logging.info('[1]')
    logging.info('Detecting vORFs in multithreading')
    with Pool(params['num_threads']) as p1:
        prodigal_results = p1.map(run_prodigal, input_files)
    
    logging.info('Performing sanity check and merging vORFs')
    with open(params['merged_orfs'], 'w') as hndl:
       n = [hndl.write(x[0]) for x in prodigal_results if x]
    gc_map = {x[1]:x[2] for x in prodigal_results if x}
    if n == 0:
        logging.error('No vORFs were found / merged')
        shutil.rmtree(params['tmp_dir'])
        sys.exit(1)
    logging.info('[2]')
    logging.info('Aligning virus-specific markers to your vORFS in multithreading (Note: This part is faster with more threads)')
    run_mmseqs(params['merged_orfs'])

    logging.info('[3]')
    logging.info('Generating queries file with the unordered collection of sets (matched virus-specific markers)')
    m8 = m8_reader(params['mmseqs_output'])
    markers_redundant = [line.split('_')[2] for line in m8['target']]
    m8['mapping'] = markers_redundant
    mapped_sets = [(x[0][:x[0].rfind('_')],set(x[1]['mapping'])) for x in m8.groupby(m8['query'])]
    queries = {genome:[] for genome in set([elem[0] for elem in mapped_sets])}
    for genome, markers_set in mapped_sets:
        queries[genome].append(set(markers_set))

    ## We comment the following two lines of code that allow to write queries to disk (for troubleshooting)
    #with open(os.path.join(params['output_dir'], 'queries.pkl'), 'wb') as f:
    #    pickle.dump(queries, f)

    logging.info('[4]')
    logging.info('Loading database, getting ready to search')
    with open(params['database_file'], 'rb') as hndl:
        database = pickle.load(hndl)

    ## We also measure the actual search wall-clock time execution
    before = time.time()

    logging.info('[5]')
    logging.info('Running virus search in multithreading')

    if params['with_replacement']:
        logging.info('Replacement option active: note that this is meaningful when the query filenames are drawn from the database sequences!')
        with Pool(params['num_threads']) as p1:
                search_results = p1.map(find_one_virus_replacement, [key for key in queries.keys()] )
    else:
        with Pool(params['num_threads']) as p1:
            search_results = p1.map(find_one_virus, [key for key in queries.keys()] )

    after = time.time()

    logging.info('[6]')
    logging.info(f'Writing results to disk at {params["results_file"]}')
    with open(params['results_file'], 'w') as hndl:
        #print('id,Realm,Kingdom,Phylum,Class,Order,Family,Genus,Species,gc_delta,score,tie_score,n_ties', end='\n', file=hndl)
        print('id,Realm,Kingdom,Phylum,Class,Order,Family,gc_delta,score,tie_score,n_ties', end='\n', file=hndl)
        for search_result in search_results:
            virus_id, chosen, score, ties = search_result
            if not chosen:
                continue
            tie_score = get_tie_score(ties)
            lineage = chosen.split(';')
            print(f"{virus_id},{','.join(lineage)},{score:3.3f},{tie_score:3.3f},{len(ties)-1}", end='\n', file=hndl)
    results = pd.read_csv(params['results_file'])
    if len(results) == 0:
        logging.error(f'Search took {after-before}s. No viruses found. Was the input correct?')
    else:
        logging.info(f'Search took {after-before:3.3}s. Taxonomy for n={len(results)} written to file. Removing temporary directory and exiting.')
    results = results.sort_values(by='id')
    results.to_csv(params['results_file'], index=False)
    shutil.rmtree(params['tmp_dir'])

    sys.stdout.write(f'\nThank you for using Virgo. If you intend to use this program in your work, please cite our paper! \n{__citation__}\n')
