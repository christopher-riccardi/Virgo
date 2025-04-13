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
import json
import time
import glob

## I'm using logging module, but also defining a plain for simpler user text at the beginning and end
logging.basicConfig(format='%(asctime)s - %(funcName)s:%(lineno)d [%(levelname)s] %(message)s', datefmt='%d-%b-%y %H:%M:%S', level=logging.INFO)
plain_logger = logging.getLogger("plain")
plain_logger.setLevel(logging.INFO)
plain_handler = logging.StreamHandler()
plain_handler.setFormatter(logging.Formatter("%(message)s"))
plain_logger.addHandler(plain_handler)
plain_logger.propagate = False

def parse_arguments():
    parser = argparse.ArgumentParser(description='')
    parser.add_argument('-d', '--data', type=str, required=True, 
                        help='Virgo database folder')
    parser.add_argument('-i', '--input', type=str, required=True, 
                        help='Folder with input query fasta files')
    parser.add_argument('-o', '--output', type=str, required=True,
                        help='Name of output directory that will contain Virgo results')
    parser.add_argument('-t', '--num_threads', type=int, required=False, default=os.cpu_count(),
                        help='Number of threads to use for multiprocessing')
    parser.add_argument('--with-replacement', action='store_true',
                        help='Skip one database entry identical to query. Useful for Leave-One-Out studies')
    parser.add_argument('--no-gc', action='store_true',
                        help='Do not use G+C content to break ties. Default is to use it')
    parser.add_argument('--min_score', type=float, required=False, default=0.05,
                        help='Show only viruses with bidirectional subsethood score > value. Must be in range (0, 1]. Default: 0.05')
    parser.add_argument('--drop-ties', action='store_true',
                        help='Do not include viruses with ties in the final result table. Default is to include them')
    parser.add_argument('--virus-by-virus', action='store_true',
                        help='Write JSON file with all query-db comparisons that yield a score > –-min_score. Default is to not write it')
    # parser.add_argument('--conf_thresh', type=float, required=False, default=0.61,
    #                     help='Queries with a confidence threshold < value do not pass confidence filter. Must be in range (0, 1]. Default: 0.61')
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
    # This is a custom function for reading .m8 files
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
    ## General parameters
    use_gc = 1 - params['no_gc']

    ## Query information
    query_seq = queries[virus_id]
    query_gc = gc_map[virus_id]
    query_vORFs_count = len(query_seq)
    if not query_vORFs_count:
        return None

    ## Scores information; {virus_id : score}
    scores_dict = {}

    ## Main loop, the costliest thing to run
    skipped = 0
    for key, value in database.items():
        if query_seq == value[0] and skipped == 0:
          skipped += 1
          continue # skip identity once
        score = bidirectional_subsethood(query_seq, value[0])
        if 0 == score: continue
        scores_dict[key] = score

    ## At this stage check sanity
    if not scores_dict:
        ## Must return
        return None

    ## Instead, at this point write files if comparisons are requested
    if params['virus_by_virus']:
        temporary_fname = os.path.join(params['tmp_dir'], f'{virus_id.replace("/", "").replace("\\", "").replace(" ", "")}.part')       
        taxonomy_scores_dict = {virus_id:{}}
        for key, value in scores_dict.items():
                database_taxon_from_id = ';'.join(database[key][1].split(';')[:6])
                if not taxonomy_scores_dict[virus_id].get(database_taxon_from_id, None):
                    taxonomy_scores_dict[virus_id][database_taxon_from_id] = []
                taxonomy_scores_dict[virus_id][database_taxon_from_id].append(value)
        with open(temporary_fname, 'w') as hndl:
            json.dump(taxonomy_scores_dict, hndl)

    max_score = max(scores_dict.values())
    if max_score <= params['min_score']:
        return None

    ## Store the keys (database virus ids) of best-scoring
    max_scoring_viruses = [key for key, value in scores_dict.items() if value==max_score]
    ## Number of viruses in the database that scored max
    database_max_scoring = len(max_scoring_viruses)
    ## Pick first candidate for now. If G+C is enabled, it gets updated later
    best_by_gc_and_score = max_scoring_viruses[0]
    gc_delta = abs(query_gc - float(database[best_by_gc_and_score][1].split(';')[-1]))

    ## Handle ties
    if database_max_scoring - 1:
        if use_gc:
            for max_scoring_virus in max_scoring_viruses:
                current_gc_diff = abs(query_gc - float(database[max_scoring_virus][1].split(';')[-1]))
                if current_gc_diff < gc_delta:
                    gc_delta = current_gc_diff
                    best_by_gc_and_score = max_scoring_virus
    
    ## Now, relative to the chosen reference from which to borrow taxonomy
    representation, information = database[best_by_gc_and_score]
    refnc_vORFs_count = len(representation)
    #refnc_vORFs_count = database_vORFs_content[best_by_gc_and_score]
    taxonomy = ';'.join(information.split(';')[:6])

    ## Relative to ties and tie scores
    number_of_ties = len(set([';'.join(database[x][1].split(';')[:6]) for x in max_scoring_viruses]))
    if number_of_ties == 0: # Always > 0
        logging.error('Number of ties equal to zero, something could be wrong in the database or database parsing. Has anything changed?')
        return None
    tie_score = 1 / number_of_ties

    ## Filter by tie_score if requester by user
    if params['drop_ties'] and tie_score < 1:
        return None

    ## Set pass confidence filter according to user input
    pass_confidence_filter = 1
    if (query_vORFs_count < 2) or (refnc_vORFs_count < 2):
        if (tie_score < 1) and (score < 0.8):
            pass_confidence_filter = 0

    return [virus_id, 
        taxonomy,
        gc_delta,
        max_score,
        tie_score,
        number_of_ties,
        query_vORFs_count,
        refnc_vORFs_count,
        pass_confidence_filter]

def find_one_virus(virus_id):
    ## General parameters
    use_gc = 1 - params['no_gc']

    ## Query information
    query_seq = queries[virus_id]
    query_gc = gc_map[virus_id]
    query_vORFs_count = len(query_seq)
    if not query_vORFs_count:
        return None

    ## Scores information; {virus_id : score}
    scores_dict = {}

    ## Main loop, the costliest thing to run
    skipped = 0
    for key, value in database.items():
        score = bidirectional_subsethood(query_seq, value[0])
        if 0 == score: continue
        scores_dict[key] = score

    ## At this stage check sanity
    if not scores_dict:
        ## Must return
        return None

    ## Instead, at this point write files if comparisons are requested
    if params['virus_by_virus']:
        temporary_fname = os.path.join(params['tmp_dir'], f'{virus_id.replace("/", "").replace("\\", "").replace(" ", "")}.part')       
        taxonomy_scores_dict = {virus_id:{}}
        for key, value in scores_dict.items():
                database_taxon_from_id = ';'.join(database[key][1].split(';')[:6])
                if not taxonomy_scores_dict[virus_id].get(database_taxon_from_id, None):
                    taxonomy_scores_dict[virus_id][database_taxon_from_id] = []
                taxonomy_scores_dict[virus_id][database_taxon_from_id].append(value)
        with open(temporary_fname, 'w') as hndl:
            json.dump(taxonomy_scores_dict, hndl)

    max_score = max(scores_dict.values())
    if max_score <= params['min_score']:
        return None

    ## Store the keys (database virus ids) of best-scoring
    max_scoring_viruses = [key for key, value in scores_dict.items() if value==max_score]
    ## Number of viruses in the database that scored max
    database_max_scoring = len(max_scoring_viruses)
    ## Pick first candidate for now. If G+C is enabled, it gets updated later
    best_by_gc_and_score = max_scoring_viruses[0]
    gc_delta = abs(query_gc - float(database[best_by_gc_and_score][1].split(';')[-1]))

    ## Handle ties
    if database_max_scoring - 1:
        if use_gc:
            for max_scoring_virus in max_scoring_viruses:
                current_gc_diff = abs(query_gc - float(database[max_scoring_virus][1].split(';')[-1]))
                if current_gc_diff < gc_delta:
                    gc_delta = current_gc_diff
                    best_by_gc_and_score = max_scoring_virus
    
    ## Now, relative to the chosen reference from which to borrow taxonomy
    representation, information = database[best_by_gc_and_score]
    refnc_vORFs_count = len(representation)
    #refnc_vORFs_count = database_vORFs_content[best_by_gc_and_score]
    taxonomy = ';'.join(information.split(';')[:6])

    ## Relative to ties and tie scores
    number_of_ties = len(set([';'.join(database[x][1].split(';')[:6]) for x in max_scoring_viruses]))
    if number_of_ties == 0: # Always > 0
        logging.error('Number of ties equal to zero, something could be wrong in the database or database parsing. Has anything changed?')
        return None
    tie_score = 1 / number_of_ties

    ## Filter by tie_score if requester by user
    if params['drop_ties'] and tie_score < 1:
        return None

    ## Set pass confidence filter according to user input
    pass_confidence_filter = 1
    if (query_vORFs_count < 2) or (refnc_vORFs_count < 2):
        if (tie_score < 1) and (score < 0.8):
            pass_confidence_filter = 0

    return [virus_id, 
        taxonomy,
        gc_delta,
        max_score,
        tie_score,
        number_of_ties,
        query_vORFs_count,
        refnc_vORFs_count,
        pass_confidence_filter]

def check_input_paths(paths):
    missing_paths = []
    for path in paths:
        if not (os.path.isfile(path) or os.path.isdir(path)):
            missing_paths.append(path)
    if missing_paths:
        logging.error(f"The following path(s) do not exist: {', '.join(missing_paths)}")
        sys.exit(1)
    else:
        return True

if __name__=='__main__':
    args = parse_arguments()
    plain_logger.info(f'This is Virgo v{__version__}\n')
    

    params = {}
    params['input_dir'] = args.input
    params['output_dir'] = args.output
    params['database_dir'] = args.data
    params['with_replacement'] = args.with_replacement
    params['no_gc'] = args.no_gc
    params['min_score'] = args.min_score
    # params['conf_thresh'] = args.conf_thresh
    params['drop_ties'] = args.drop_ties
    params['virus_by_virus'] = args.virus_by_virus
    params['virus_specific_markers'] = os.path.join(params['database_dir'], 'DB')
    params['database_file'] = os.path.join(params['database_dir'], 'database.pkl')
    params['tmp_dir'] = os.path.join(params['output_dir'], 'tmp_dir')
    params['mmseqs_output'] = os.path.join(params['tmp_dir'], 'output.m8')
    params['num_threads'] = args.num_threads
    params['merged_orfs'] = os.path.join(params['tmp_dir'], 'merged.faa')
    params['results_file'] = os.path.join(params['output_dir'], 'results.csv')
    params['virus_by_virus_file'] = os.path.join(params['output_dir'], 'virus_by_virus.json')

    logging.info(f"Parameters set:\n\
          input_dir: {params['input_dir']},\n\
          output_dir: {params['output_dir']},\n\
          database_dir: {params['database_dir']},\n\
          virus_specific_markers: {params['virus_specific_markers']},\n\
          database_file: {params['database_file']},\n\
          replacement: {params['with_replacement']},\n\
          no_gc: {params['no_gc']},\n\
          min_score: {params['min_score']},\n\
          drop_ties: {params['drop_ties']},\n\
          virus_by_virus: {params['virus_by_virus']},\n\
          tmp_dir: {params['tmp_dir']},\n\
          mmseqs_output: {params['mmseqs_output']},\n\
          num_threads: {params['num_threads']},\n\
          merged_orfs: {params['merged_orfs']},\n\
          results_file: {params['results_file']}\n\
          virus_by_virus_file: {params['virus_by_virus_file']}\n\
          "
    )
    
    ## Checking user input
    check_input_paths([params['input_dir'], params['database_dir']])

    if 1 <= params['min_score'] < 0:
        logging.error('--min_score command line argument must be in range (0, 1]')
        sys.exit(1)

    # if 1 <= params['conf_thresh'] < 0:
    #     logging.error('--conf_thresh command line argument must be in range (0, 1]')
    #     sys.exit(1)

    logging.info('[0]')
    if CreateDirectory(params['output_dir']) == 1:
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

    ##==============================================================================================##
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

    ##==============================================================================================##
    logging.info('[2]')
    logging.info('Aligning virus-specific markers to your vORFs in multithreading (Note: This part is faster with more threads)')
    run_mmseqs(params['merged_orfs'])


    ##==============================================================================================##
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


    ##==============================================================================================##
    logging.info('[4]')
    logging.info('Loading database, getting ready to search')
    with open(params['database_file'], 'rb') as hndl:
        database = pickle.load(hndl)

    ## Store reference vORFs count
    database_vORFs_content = {}
    for key, value in database.items():
        vORFs = len(value[0])
        database_vORFs_content[key] = vORFs


    ## We also measure the actual search wall-clock time execution
    before = time.time()


    ##==============================================================================================##
    logging.info('[5]')
    logging.info('Running virus search in multithreading')

    if params['with_replacement']:
        logging.info('Replacement option active: will skip the first best-scoring database entry')
        with Pool(params['num_threads']) as p1:
                search_results = p1.map(find_one_virus_replacement, [key for key in queries.keys()] )
    else:
        with Pool(params['num_threads']) as p1:
            search_results = p1.map(find_one_virus, [key for key in queries.keys()] )

    after = time.time()


    ##==============================================================================================##
    logging.info('[6]')
    logging.info(f'Pooling results, will write to file: {params["results_file"]}')
    with open(params['results_file'], 'w') as hndl:
        print('id,Realm,Kingdom,Phylum,Class,Order,Family,gc_delta,score,tie_score,n_ties,query_vORFs_count,refnc_vORFs_count,pass_confidence_filter', end='\n', file=hndl)
        for search_result in search_results:
            if not search_result:
                continue
            virus_id, taxonomy, gc_delta, max_score, tie_score, number_of_ties, query_vORFs_count, refnc_vORFs_count, pass_confidence_filter = search_result
            
            lineage = taxonomy.split(';')
            print(f"{virus_id},{','.join(lineage)},{gc_delta:3.3f},{max_score:3.3f},{tie_score:3.3f},{number_of_ties:3.3f},{query_vORFs_count},{refnc_vORFs_count},{pass_confidence_filter}",end='\n',file=hndl)
    results = pd.read_csv(params['results_file'])
    if len(results) == 0:
        logging.error(f'Search took {after-before}s. No viruses found. Was the input correct?')
    else:
        logging.info(f'Search took {after-before:3.3}s. Taxonomy for n={len(results)} written to file. Removing temporary directory and exiting.')

    results = results.sort_values(by='id')
    results.to_csv(params['results_file'], index=False)
    if params['virus_by_virus']:
        merged_part_files = {}
        part_files = glob.glob(params['tmp_dir']+'/*.part')
        for part_file in part_files:
            part_dict = json.load(open(part_file))
            key = next(iter(part_dict))
            values = part_dict[key]
            merged_part_files[key] = values
        with open(params['virus_by_virus_file'], 'w') as hndl:
            json.dump(merged_part_files, hndl, indent=4)
    shutil.rmtree(params['tmp_dir'])

    plain_logger.info(f'\nThank you for using Virgo. If you intend to use this program in your work, please cite our paper! \n{__citation__}\n')
