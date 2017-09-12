from bio_hansel.utils import exc_exists


def test_exc_exists():
    which_exists = exc_exists('which')
    assert which_exists, 'which should exist on most/all Linux systems'
    makeblastdb_exists = exc_exists('makeblastdb')
    assert makeblastdb_exists, 'makeblastdb from BLAST+ should be installed in $PATH'
    blastn_exists = exc_exists('blastn')
    assert blastn_exists, 'blastn from BLAST+ should be installed in $PATH'
    jellyfish_exists = exc_exists('jellyfish')
    assert jellyfish_exists, 'jellyfish should be installed in $PATH'
