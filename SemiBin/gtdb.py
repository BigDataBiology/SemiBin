import os
from os import path
import shutil
import sys
from time import sleep

def download_GTDB_to(logger, GTDB_dir):
    """
    Download GTDB.
    logger: logger
    GTDB_dir: where to store the data
    """
    import requests
    import tarfile
    import hashlib
    logger.info(f'Downloading GTDB to {GTDB_dir}.  It will take a while..')
    os.makedirs(GTDB_dir, exist_ok=True)

    download_path = os.path.join(GTDB_dir, 'GTDB_v95.tar.gz')

    # This allows one to bypass the actual downloading for testing (otherwise
    # tests take a few hours!)
    if path.exists('__SemiBin__internal_test_GTDB_v95.tar.gz'):
        logger.warning(f'Found __SemiBin__internal_test_GTDB_v95.tar.gz')
        logger.warning(f'Will copy it to {download_path} and use it as GTDB tarball')
        shutil.copy2('__SemiBin__internal_test_GTDB_v95.tar.gz', download_path)
        md5_ok = True
    else:
        download_url = 'https://zenodo.org/record/4751564/files/GTDB_v95.tar.gz?download=1'
        with requests.get(download_url, stream=True) as r:
            m = hashlib.md5()
            with open(download_path, 'wb') as f:
                while True:
                    data = r.raw.read(4096)
                    if not data:
                        break
                    f.write(data)
                    m.update(data)
        logger.info('Download finished. Checking MD5...')
        md5_ok = m.hexdigest() == '4a70301c54104e87d5615e3f2383c8b5'
    if md5_ok:
        try:
            tar = tarfile.open(download_path, "r:gz")
            file_names = tar.getnames()
            for file_name in file_names:
                tar.extract(file_name, GTDB_dir)
            tar.close()
        except Exception:
            sys.stderr.write(
                f"Error: cannot untar the file.")
            sys.exit(1)

        os.remove(download_path)
    else:
        os.remove(download_path)
        sys.stderr.write(
            f"Error: MD5 check failed. Downloading GTDB database failed (Please check the internet connections or storage space), removing '{download_path}'.\n")
        sys.exit(1)


def find_or_download_gtdb(logger, GTDB_reference_dir, force):
    '''Find or (lazily) download GTDB

    logger: logger
    GTDB_reference_dir: str
        where to store the database. Can be None, in which case, the default will be used
    force: bool
        If a database is found and `force` is true, then remove it and redownload it.
    '''
    using_default_dir = False
    if GTDB_reference_dir is None:
        cache_home = os.environ.get('XDG_CACHE_HOME')
        if cache_home is None:
            cache_home = path.join(os.environ['HOME'], '.cache')

        GTDB_reference_dir = os.path.join(
                cache_home,
                'SemiBin',
                'mmseqs2-GTDB')
        using_default_dir = True
    GTDB_reference = path.join(GTDB_reference_dir, 'GTDB')
    if path.exists(GTDB_reference):
        if force:
            logger.warning(f'Found existing directory {GTDB_reference_dir}, but `--force` was used and will be removed before re-downloading GTDB')
            if sys.stdin and sys.stdin.isatty():
                print("Waiting 5 seconds to remove...")
                for i in range(5):
                    print(f'{5-i}...')
                    sleep(1)
            os.unlink(GTDB_reference)
        else:
            logger.info(f'Found GTDB directory: `{GTDB_reference_dir}`.')
    if not path.exists(GTDB_reference):
        download_GTDB_to(logger, GTDB_reference_dir)
        if using_default_dir:
            with open(path.join(GTDB_reference_dir, 'CACHEDIR.TAG'), 'wt') as cachedir_tag:
                cachedir_tag.write('Signature: 8a477f597d28d172789f06886806bc55\n')
                cachedir_tag.write('# This file is a cache directory tag created SemiBin.\n')
                cachedir_tag.write('# For information about cache directory tags, see:\n')
                cachedir_tag.write('# http://www.brynosaurus.com/cachedir/\n')

    return GTDB_reference
