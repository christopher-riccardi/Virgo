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
import sqlite3
import os, sys
import pickle
import shutil
import glob

logging.basicConfig(format='%(asctime)s - %(funcName)s:%(lineno)d [%(levelname)s] %(message)s', datefmt='%d-%b-%y %H:%M:%S', level=logging.INFO)

def parse_arguments():
    parser = argparse.ArgumentParser(description='Create Virgo database for ICTV-based virus classification')
    parser.add_argument('--sql', type=str, required=True, 
                        help='SQL database created by ICTVdump')
    parser.add_argument('--taxonomy_table', type=str, required=True, 
                        help='Taxonomy table created by ICTVdump')
    parser.add_argument('--virus_markers', type=str, required=True, 
                        help='Markers dataset folder')
    parser.add_argument('-o', '--output', type=str, required=True,
                        help='Name of output directory that will contain Virgo database files')
    parser.add_argument('-t', '--num_threads', type=int, required=False,
                        help='Number of threads to use for multiprocessing')
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

def run_prodigal(input_fasta):
    sequence_name = os.path.basename(input_fasta).replace('.fa', '')
    output_fasta = input_fasta.replace('.fa', '.faa')
    prodigal_command = [
        'prodigal-gv',
        '-p', 'meta',
        '-i', input_fasta,
        '-a', output_fasta
    ]
    try:
        result = subprocess.run(
            prodigal_command,
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
            text=True
        )
        if result.returncode == 0:
            gc_contents = []
            weights = []
            try: seq = Seq2Dict(output_fasta)
            except: return None
            for key, value in seq.items():
                gc = float(key.split('gc_cont=')[1])
                size = len(value)
                gc_contents.append(gc)
                weights.append(size)
            
            orf_in_memory = [line for line in open(output_fasta) if line != '\n']
            index_i = 1
            for i, line in enumerate(orf_in_memory):
                if line.startswith('>'):
                    orf_in_memory[i] = '>' + sequence_name + '_' + str(index_i) + '\n'
                    index_i += 1
            if index_i == 1:
                logging.warning(f'No ORFs detected in sequence {input_fasta}, skipping this entry.')
                return None
            ## Update ORFs headers inside the output_fasta
            with open(output_fasta, 'w') as hndl:
                for line in orf_in_memory:
                    hndl.write(line)
            weighted_avg = np.average(gc_contents) #, weights=weights)
            return (sequence_name, weighted_avg)
        else:
            logging.error(f'Error running Prodigal: {result.stderr}.')
            return None
    except Exception as e:
        logging.error(f'An error occurred: {e}.')
    return None

def run_mmseqs(input_fasta, virus_specific_markers, mmseqs_output, tmp_dir, num_threads):
    cmd = [
        'mmseqs', 'easy-search',
        input_fasta,
        virus_specific_markers,
        mmseqs_output,
        tmp_dir,
        '-s', '7.5', #7
        '-e', '1e-3', #1e-3
        '-c', '0.2', #0.2
        '--cov-mode', '1',
        '--threads', str(num_threads)
    ]
    try:
        result = subprocess.run(
            cmd,
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
            text=True
        )
        if result.returncode == 0:
            return 0
        else:
            logging.info("Error running MMSeqs2:")
            logging.info(result.stderr)
    except Exception as e:
        logging.info(f"An error occurred: {e}.")
    return 1

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
    
#########################
## Quite procedural
#########################
args = parse_arguments()

params = {}
params['tmp_dir'] = os.path.join(args.output, 'tmp_dir')

check_input_paths([args.sql, args.taxonomy_table, args.virus_markers])

if CreateDirectory(args.output) == 1:
    sys.exit(1)

if CreateDirectory(params['tmp_dir']) == 1:
    sys.exit(1)

## (3) connect to database
logging.info('Connecting to database using sqlite3')
conn = sqlite3.connect(args.sql)
cursor = conn.cursor()
logging.info('Reading sequence information from database')
xlsx = pd.read_csv(args.taxonomy_table)
rows = xlsx.iterrows()
fasta_files = []
for i, row in rows:
    accessions = sorted(row['Virus GENBANK accession'].split(';'))
    placeholders = ','.join('?' for _ in accessions)
    cursor.execute(f'SELECT sequence FROM sequences WHERE id IN ({placeholders})', accessions)
    results = cursor.fetchall()
    sequences = [x[0] for x in results] # split the tuple
    fasta_file = os.path.join(params['tmp_dir'], row['header']+'.fa')
    with open(fasta_file, 'w') as hndl:
        hndl.write('\n'.join(sequences))
    fasta_files.append(fasta_file)
logging.info('Closing connection')
cursor.close()

logging.info('Inferring vORFs in reference genomes')
with Pool(args.num_threads) as p1:
    prodigal_results = p1.map(run_prodigal, fasta_files)

## Here prodigal returns tuples with sequence name : weighted average gc content, weighted on the vORF size
gc_map = {x[0]:x[1] for x in prodigal_results if x}

## Not all fasta files may have survived, we keep only those with at least one match to the virus-specific dataset
fasta_files = [x for x in glob.glob(params['tmp_dir'] + '/*.faa') if os.path.basename(x).replace('.faa', '') in gc_map.keys()]

logging.info('Merging ORFs for faster computation')
## Merge files
with open(os.path.join(params['tmp_dir'], 'merged_orfs.faa'), 'w') as hndl:
    for file in fasta_files:
        hndl.write(open(file).read())

logging.info('Mapping vORFs to virus-specific markers')
run_mmseqs(os.path.join(params['tmp_dir'], 'merged_orfs.faa'), 
            args.virus_markers +'/DB', 
            os.path.join(params['tmp_dir'], 'output.m8'), 
            params['tmp_dir'],
            args.num_threads)

logging.info('Generating Virgo database file')
m8 = m8_reader(os.path.join(params['tmp_dir'], 'output.m8'))
markers_redundant = [line.split('_')[2] for line in m8['target']]
m8['mapping'] = markers_redundant
mapped_sets = [(x[0][:x[0].rfind('_')],set(x[1]['mapping'])) for x in m8.groupby(m8['query'])]
database = {genome:[[]] for genome in set([elem[0] for elem in mapped_sets])}
for genome, markers_set in mapped_sets:
    database[genome][0].append(markers_set)

for genome in database.keys():
    lineage = list(xlsx.loc[(xlsx['header']==genome), ['Realm', 'Kingdom', 'Phylum', 'Class', 'Order', 'Family', 'Genus', 'Species']].values[0])
    data_string = ';'.join([str(x) for x in lineage] + [f'{gc_map[genome]:3.3f}'])
    database[genome].append(data_string)


with open(os.path.join(args.output, 'database.pkl'), 'wb') as hndl:
    pickle.dump(database, hndl)

logging.info('Cleaning up temporary files and copying virus-specific markers MMseqs database into the output folder')
shutil.rmtree(params['tmp_dir'])
for elem in glob.glob(args.virus_markers + '/*'):
    shutil.copy(elem, os.path.join(args.output, os.path.basename(elem)))
logging.info('Execution finished')
