#!/usr/bin/env python3

import os
import re
import gzip
import logging
import argparse
from pathlib import Path
from multiprocessing import Pool
from functools import partial

import skbio

from strainge.io.utils import open_compressed

logging.basicConfig()
logger = logging.getLogger()
logger.setLevel(logging.DEBUG)

special_chars = re.compile(r'[^a-zA-Z0-9\-]+')
trailing_underscore = re.compile(r'_+$')
leading_underscore = re.compile(r'^_+')


def generate_friendly_id(fpath):
    genus = fpath.parts[-4]
    species = fpath.parts[-3]
    strain = fpath.parts[-2]

    full_species = f"{genus}_{species}"

    if strain.startswith(full_species):
        strain = strain[len(full_species):]

    genus = genus[:4]
    strain = special_chars.sub('_', strain)
    strain = trailing_underscore.sub('', strain)
    strain = leading_underscore.sub('', strain)

    return f"{genus}_{species}_{strain}"

def read_assembly_report(f):
    attr = {}
    for line in f:
        if not line.startswith('#'):
            continue

        parts = line.split(':')

        if len(parts) != 2:
            continue

        key = parts[0][2:].strip()
        value = parts[1].strip()
        attr[key] = value

    return attr

def split_chrom_plasmids(path):
    chromosomes = []
    plasmids = []
    with open_compressed(path) as f:
        for entry in skbio.io.read(f, "fasta"):
            desc = entry.metadata['description']
            logger.debug('%s %s (length: %d)', entry.metadata['id'], desc,
                         len(entry))

            if "plasmid" in desc or len(entry) < 500e3:
                plasmids.append(entry)
            else:
                chromosomes.append(entry)

    logger.info('# Chrom: %d, # Plasmids: %d', len(chromosomes), len(plasmids))

    prefix = path.name.replace(".fna", "").replace(".gz", "")
    chrom_path = path.parent / f"{prefix}.chrom.fna.gz"
    with gzip.open(chrom_path, "wt") as f:
        for chrom in chromosomes:
            skbio.io.write(chrom, "fasta", into=f)

    plasmid_path = path.parent / f"{prefix}.plasmids.fna.gz"
    with gzip.open(plasmid_path, "wt") as f:
        for plasmid in plasmids:
            skbio.io.write(plasmid, "fasta", into=f)
            
def get_paths(directory):
    used_names = set()
    paths = []
    for fpath in directory.glob("**/*_genomic.fna.gz"):
        if fpath.name.endswith('_cds_from_genomic.fna.gz'):
            continue

        if fpath.name.endswith('_rna_from_genomic.fna.gz'):
            continue

        friendly_id = generate_friendly_id(fpath)
        assembly_report_path = fpath.parent / fpath.name.replace(
            "_genomic.fna.gz", "_assembly_report.txt")

        with assembly_report_path.open() as f:
            attrs = read_assembly_report(f)
        accession = attrs['RefSeq assembly accession']

        if friendly_id in used_names:
            friendly_id += f"_{accession}"

        used_names.add(friendly_id)
        
        paths.append([friendly_id,  accession, fpath])
        logger.info(f'{friendly_id} ({accession}): {fpath}')
        
    return paths

def make_symlink(symlink: Path, original: Path, overwrite = True):
    if overwrite and symlink.is_symlink():
        symlink.unlink()
    symlink.symlink_to(original)    

def prepare_db(x, output_dir, split_plasmids=False):
    friendly_id, accession, fpath = x
    target_path = output_dir / f"{friendly_id}.fa.gz"

    if split_plasmids:
        prefix = fpath.name.replace(".fna.gz", "")
        chrom_path = fpath.parent / f"{prefix}.chrom.fna.gz"

        if not chrom_path.is_file():
            logger.info("No FASTA file found with chromosomes only.")
            logger.info("Generating %s", chrom_path)
            split_chrom_plasmids(fpath)
        else:
            logger.info("Found chromosome only FASTA file at %s",
                        chrom_path)

        fpath = chrom_path

    # symlink
    make_symlink(target_path, fpath.resolve())
    return [friendly_id, str(target_path.resolve()), accession]

        
def main():
    class CustomFormatter(argparse.ArgumentDefaultsHelpFormatter,
                      argparse.RawDescriptionHelpFormatter):
        pass
    
    desc = """Search a directory recursively and collect all reference genomes,
then symlink each reference genome with a human-readable filename to the
given output directory. This directory then serves as input for a StrainGST
database.

The directory structure should follow the `ncbi-genome-download` human readable
output structure.

This script can optionally split chromosomes and plasmids too.
"""   
    parser = argparse.ArgumentParser(
        description=desc,
        formatter_class=CustomFormatter
    )

    parser.add_argument(
        'directories', type=Path, nargs='+',
        help="Directories to search for references"
    )

    parser.add_argument(
        '-s', '--split-plasmids', action="store_true", default=False,
        help="Only use chromosome scaffolds for the database."
    )

    parser.add_argument(
        '-o', '--output-dir', type=Path, required=True,
        help="Output directory. Will contain symlinks to the all references "
             "and a TSV file with some metadata on the included references."
    )
    parser.add_argument(
        '-t', '--threads', type=int, default=1,
        help="Number of parallel processes."
    )
    
    args = parser.parse_args()

    args.output_dir.mkdir(parents=True, exist_ok=True)
    
    # getting directory info
    if args.threads > 1:
        pool = Pool(args.threads)
        paths = pool.map(get_paths, args.directories)
    else:
        paths = map(get_paths, args.directories)        
    paths = [y for x in paths for y in x]
    
    # preparing db
    func = partial(prepare_db,
                   output_dir=args.output_dir,
                   split_plasmids=args.split_plasmids)
    if args.threads > 1:
        res = pool.map(func, paths)
    else:
        res = map(func, paths)

    # writing metadata
    metadata_tsv = args.output_dir / "references_meta.tsv"
    with metadata_tsv.open('w') as ofile:
        for x in res:
            print(x, sep='\t', file=ofile)


if __name__ == '__main__':
    main()
