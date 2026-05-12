# runner.py: Pipeline runner for Euclid surface photometry

import os
import logging
import logging.handlers
from os.path import join, abspath, basename
from datetime import datetime
import yaml
import pandas as pd
from tqdm import tqdm
import multiprocessing as mp
import argparse
import traceback
from astropy.table import Table
from .pipeline import *
import warnings
warnings.filterwarnings("ignore")

__version__ = '0.1.2'

# Description: Pipeline script for Euclid profile analysis
# parallelize with pool
# @pmsastro

import os
import logging
import logging.handlers
from os.path import join, abspath, basename
from datetime import datetime
import yaml
import pandas as pd
from tqdm import tqdm
import multiprocessing as mp
import argparse
import traceback
from astropy.table import Table
from euclid_pipelines import *
import warnings
warnings.filterwarnings("ignore")

__version__ = '0.1.2'

# ---------------------------------------------------------------------------
# Logging helpers
# ---------------------------------------------------------------------------

def _build_handlers(log_cfg):
    """Return a list of logging.Handler objects from a logging config dict."""
    fmt = logging.Formatter(
        fmt=log_cfg.get('format', '%(asctime)s [%(levelname)-8s] %(processName)s - %(message)s'),
        datefmt=log_cfg.get('datefmt', '%Y-%m-%d %H:%M:%S'),
    )
    handlers = []
    # Always add a stream (console) handler
    sh = logging.StreamHandler()
    sh.setFormatter(fmt)
    handlers.append(sh)
    # Optionally add a rotating file handler
    log_file = log_cfg.get('file')
    if log_file:
        fh = logging.handlers.RotatingFileHandler(
            log_file, maxBytes=10 * 1024 * 1024, backupCount=5, encoding='utf-8'
        )
        fh.setFormatter(fmt)
        handlers.append(fh)
    return handlers


def setup_logging(log_cfg):
    """Configure the root logger for the *main* process.

    Returns a ``logging.handlers.QueueListener`` that must be started and
    stopped by the caller so that worker-process log records are forwarded
    to the real handlers without file-locking conflicts.
    """
    level = getattr(logging, log_cfg.get('level', 'INFO').upper(), logging.INFO)
    # Queue used by worker processes to send records back to main process
    log_queue = mp.Queue(-1)
    handlers = _build_handlers(log_cfg)
    listener = logging.handlers.QueueListener(
        log_queue, *handlers, respect_handler_level=True
    )
    # Root logger in the main process just forwards to the queue
    root = logging.getLogger()
    root.setLevel(level)
    qh = logging.handlers.QueueHandler(log_queue)
    root.handlers = [qh]
    return listener, log_queue


def _worker_logging_init(log_queue, level):
    """Initialiser executed once per worker process.

    Redirects the root logger to the shared queue so all records are handled
    centrally in the main process.
    """
    root = logging.getLogger()
    root.handlers = [logging.handlers.QueueHandler(log_queue)]
    root.setLevel(level)

def read_yaml(file_path):
    with open(file_path, 'r') as stream:
        config = yaml.safe_load(stream)
    return config


logger = logging.getLogger(__name__)

def euclid_pipeline(data):
    """Euclid pipeline to analyse a single galaxy.

    Parameters
    ----------
    data : tuple
        ``(object_id, config)`` where *config* is the parsed YAML dict.

    Returns
    -------
    dict
        Dictionary with ``object_id`` and boolean status flags for each module
        (``M1`` – ``M5``).
    """
    log = logging.getLogger(__name__)

    # Read data to process
    object_id = data[0]
    config = data[1]

    # Define parameters of pipeline
    if 'n_radius' in config['parameters']: 
        n_radius = config['parameters']['n_radius'] 
    else:
        n_radius = 5
    
    if 'sb_threshold' in config['parameters']:
        sb_threshold = config['parameters']['sb_threshold']
    else: 
        sb_threshold = None
    if 'max_size' in config['parameters']:
        max_size = config['parameters']['max_size']
    else:
        max_size = 1000
    max_size *= 0.1/3600  # pixel --> degrees 

    outPath = join(config['absolute_paths']['output_path'],f'source_{object_id}')
    if not os.path.isdir(outPath): os.makedirs(outPath)

    # Read master catalog
    cols = ['object_id', 'segmentation_map_id', 'right_ascension', 'declination', 'kron_radius', 'phz_pp_median_redshift']
    master = pd.read_csv(config['catalogs']['objects'], usecols=cols)
    galind = np.argwhere(master['object_id'] == object_id)[0][0]
    master = master.iloc[galind]
    
    # Get object information
    ra,dec = master['right_ascension'], master['declination']
    segmentation_map_id = master['segmentation_map_id']
    tile = np.int64(str(segmentation_map_id)[:9])
    kron_radius = master['kron_radius']*0.1/3600   # pixel --> degrees
    redshift = master['phz_pp_median_redshift']
    
    # Define data path
    dataPath = join(config['absolute_paths']['data_path'],f'{tile}')
    
    # Define status of each pipeline step
    STATUSP1,STATUSP2,STATUSP3,STATUSP4,STATUSP5 = [False]*5

    log.debug('Starting pipeline for object %s (z=%.3f, tile=%s)', object_id, redshift, tile)

    try:
        # Run M1 to create cutout and mask
        if config['pipeline_steps']['M1']:
            log.debug('[%s] M1 – creating cutout and mask', object_id)
            size = n_radius*kron_radius if n_radius*kron_radius < max_size else max_size
            cutoutFile, maskFile, (ra,dec) = euclid_module1(object_id,ra,dec,
                                    size, dataPath, outPath)
            STATUSP1 = True
            log.debug('[%s] M1 done – cutout: %s', object_id, cutoutFile)
        else:
            cutoutFile = find(outPath, f'*{object_id}*VIS_cutout.fits')[0]
            maskFile = find(outPath, f'*{object_id}*mask.fits')[0]
            log.debug('[%s] M1 skipped – using existing files', object_id)

        # Run M2 to create master profile
        if config['pipeline_steps']['M2']:
            log.debug('[%s] M2 – fitting isophotes (VIS)', object_id)
            master = euclid_module2(cutoutFile, maskFile, ra, dec, object_id)
            STATUSP2 = True
            log.debug('[%s] M2 done – %d isophotes fitted', object_id, len(master))
        else:
            masterFile = find(outPath, f'*{object_id}*master-profile.csv')[0]
            master = Table.read(masterFile, format='ascii.ecsv')
            log.debug('[%s] M2 skipped – loaded %s', object_id, masterFile)

        # Run M3 to measure profiles on different bands
        if config['pipeline_steps']['M3']:
            log.debug('[%s] M3 – measuring multi-band SB profiles', object_id)
            profiles = euclid_module3(object_id, master, maskFile)
            STATUSP3 = True
            log.debug('[%s] M3 done – %d radial bins', object_id, len(profiles))
        else:
            profilesFile = find(outPath, f'*{object_id}*profiles.csv')[0]
            profiles = Table.read(profilesFile, format='ascii.ecsv')
            log.debug('[%s] M3 skipped – loaded %s', object_id, profilesFile)

        # Run M4 to measure photometry on the profiles
        if config['pipeline_steps']['M4']:
            log.debug('[%s] M4 – deriving photometric parameters', object_id)
            photometry = euclid_module4(profiles, object_id, outPath, z=redshift)
            STATUSP4 = True
            log.debug('[%s] M4 done', object_id)
        else:
            photometryFile = find(outPath, f'*{object_id}*photometry.csv')[0]
            photometry = Table.read(photometryFile, format='csv')
            log.debug('[%s] M4 skipped – loaded %s', object_id, photometryFile)

        # Run M5 to classify the breaks in the profiles
        if config['pipeline_steps']['M5']:
            log.debug('[%s] M5 – classifying disc breaks', object_id)
            breaks = euclid_module5(profiles, photometry, object_id, outPath, z=redshift, sb_threshold=sb_threshold)
            STATUSP5 = True
            log.debug('[%s] M5 done – classification: %s',
                      object_id, breaks.get('final_classification', '?') if isinstance(breaks, dict) else '?')

        log.info('[%s] completed  M1=%s M2=%s M3=%s M4=%s M5=%s',
                 object_id, STATUSP1, STATUSP2, STATUSP3, STATUSP4, STATUSP5)

    except Exception as e:
        log.error('[%s] pipeline error: %s', object_id, e)
        log.debug('[%s] traceback:\n%s', object_id, traceback.format_exc())
    
    results = {'object_id':object_id, 'M1':STATUSP1, 'M2':STATUSP2, 'M3':STATUSP3, 'M4':STATUSP4, 'M5':STATUSP5}
    return results 

def save_results(results, filename='profile_results.csv'):
    """Save intermediate or final results to file."""
    tbl = {k: [] for k in results[0].keys()}
    for res in results:
        for key in res:
            tbl[key].append(res[key])
    df = pd.DataFrame(tbl)
    df.to_csv(filename, index=False)
    return filename

def parallel_run_pipeline(dataset, num_processes, log_queue=None, log_level=logging.INFO):
    """Run ``euclid_pipeline`` in parallel using a ``multiprocessing.Pool``.

    Parameters
    ----------
    dataset : list
        List of ``(object_id, config)`` tuples.
    num_processes : int
        Number of worker processes.
    log_queue : multiprocessing.Queue or None
        Logging queue shared with the ``QueueListener`` in the main process.
        When *None* workers fall back to basic logging.
    log_level : int
        Logging level forwarded to each worker initialiser.

    Returns
    -------
    list
        List of result dicts from ``euclid_pipeline``.
    """
    initargs = (log_queue, log_level) if log_queue is not None else ()
    initfunc  = _worker_logging_init   if log_queue is not None else None
    with mp.Pool(processes=num_processes,
                 initializer=initfunc,
                 initargs=initargs) as pool:
        results = list(tqdm(
            pool.imap_unordered(euclid_pipeline, dataset, chunksize=1),
            total=len(dataset),
        ))
    return results


def make_parser():
    """Create an argument parser for the Euclid Pipeline."""
    parser = argparse.ArgumentParser(description='Run Euclid profile pipeline.')
    parser.add_argument('configFile', type=str, help='Path to the YAML configuration file.')
    parser.add_argument('-n', type=int, default=6, help='Number of processes to use for parallel execution.')
    return parser

def main():
    # Create the argument parser
    parser = make_parser().parse_args()

    # Read yaml file
    config = read_yaml(parser.configFile)

    # -----------------------------------------------------------------------
    # Logging setup
    # -----------------------------------------------------------------------
    log_cfg = config.get('logging', {})
    log_level = getattr(logging, log_cfg.get('level', 'INFO').upper(), logging.INFO)
    listener, log_queue = setup_logging(log_cfg)
    listener.start()

    # After setup_logging the root logger already has a QueueHandler, so
    # module-level `logger` will work from this point onward.
    logger.info('Euclid Pipeline v%s', __version__)
    logger.info('Config file  : %s', parser.configFile)
    logger.info('Log level    : %s', log_cfg.get('level', 'INFO').upper())
    logger.info('Log file     : %s', log_cfg.get('file', 'none'))

    # Create dataset to pass to parallel function
    table = pd.read_csv(config['catalogs']['ids'])
    objects_id = table['object_id']
    dataset = [(object_id, config) for object_id in objects_id]

    logger.info('Processors   : %d', parser.n)
    logger.info('Galaxies     : %d', len(dataset))
    time_init = datetime.now()
    logger.info('Starting at  : %s', time_init.strftime('%Y-%m-%d %H:%M'))

    results = parallel_run_pipeline(dataset, parser.n,
                                    log_queue=log_queue, log_level=log_level)

    # Save results
    resultFile = save_results(
        results,
        join(dirname(__file__),
             f'profile_results_{datetime.now().strftime("%Y-%m-%d_%H%M")}.csv')
    )
    logger.info('Results saved to: %s', resultFile)

    # Execution time
    td = datetime.now() - time_init
    hours, remainder = divmod(td.seconds, 3600)
    minutes, seconds  = divmod(remainder, 60)
    logger.info('Execution time: %d days %02dh:%02dm:%02ds',
                td.days, hours, minutes, seconds)
    logger.info('Pipeline completed successfully.')

    # Shut down the log listener cleanly
    listener.stop()

if (__name__ == '__main__'):
    main()
    


