#!/usr/bin/env python

import os
import urllib
import csv
from subprocess import call

def load_assemblies():
    """fetch bacteria assembly_summary.txt file from NCBI refseq"""
    summary = "assembly_summary.txt"
    if not os.path.exists(summary):
        print("Fetching assembly file")
        urllib.urlretrieve("ftp://ftp.ncbi.nlm.nih.gov/genomes/refseq/bacteria/" + summary, summary)
    with open(summary, 'r') as sum:
        sum.readline()
        assemblies =[row for row in csv.DictReader(sum, delimiter='\t')]
    return assemblies

def assembly_filter(ass):
    """Do we want this assembly?"""
    if (ass['version_status'] == "latest" and ass['assembly_level'] == 'Complete Genome') \
        or ass['refseq_category'] == "reference genome":
        name = ass['organism_name'].lower()
        if name.startswith("escherichia coli") or name.startswith("shigella"):
            return True
    return False

def rsync_assembly(ass):
    uri = ass['ftp_path'].replace("ftp:", "rsync:")
    dirname = assembly_dir(ass)
    if os.path.exists(dirname):
        return
    command = ["rsync", "-av", uri, "."]
    print "Running", ' '.join(command)
    call(command)
    return dirname

def link_assembly(ass):
    """Create a link to assembly using descriptive strain name"""
    dirname = assembly_dir(ass)
    source = os.path.join(dirname, dirname + "_genomic.fna.gz")
    if os.path.exists(source):
        dest = os.path.join("assemblies", assembly_name(ass) + ".fasta.gz")
        if not os.path.exists("assemblies"):
            os.mkdir("assemblies")
        if not os.path.exists(dest):
            print 'Link from', source, 'to', dest
            os.symlink(os.path.join("..", source), dest)

def assembly_dir(ass):
    return ass['ftp_path'].split('/')[-1]

assembly_names = set()

def assembly_name(ass):
    """Generate a name representing the name and suitable for use as a filename"""
    org = ass['organism_name']
    strain = ass['infraspecific_name']
    isolate = ass['isolate']

    org = org.replace("Escherichia", "E")
    org = org.replace("Shigella", "S")
    strain = strain.replace("strain=", "")
    name = org
    if strain and name.find(strain) < 0:
        name += "_" + strain
    if isolate and name.find(isolate) < 0:
        name += "_" + isolate
    name = name.replace(".", "")
    name = name.replace("/", "-")
    name = name.replace("(", "")
    name = name.replace(")", "")
    name = name.replace("'", "")
    name = name.replace(";", "-")
    name = name.replace(":", "-")
    name = name.replace(" ", "_")
    name = name.replace("K-12_K-12", "K-12")
    if name in assembly_names:
        name += "_" + ass['# assembly_accession'].split('.')[0]
    assembly_names.add(name)
    # print (org, strain, isolate), name
    return name

###################################
# MAIN
###################################

assemblies = load_assemblies()
print len(assemblies), "assemblies"

keepers = [assembly for assembly in assemblies if assembly_filter(assembly)]
print len(keepers), "keepers"

for keeper in keepers:
    rsync_assembly(keeper)
    link_assembly(keeper)