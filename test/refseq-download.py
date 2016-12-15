#!/usr/bin/env python

import sys
import urllib
import csv
from subprocess import call

def load_assemblies():
    # fetch bacteria assembly_summary.txt file
    summary = "assembly_summary.txt"
    #urllib.urlretrieve("ftp://ftp.ncbi.nlm.nih.gov/genomes/refseq/bacteria/" + summary, summary)
    with open(summary, 'r') as sum:
        sum.readline()
        assemblies =[row for row in csv.DictReader(sum, delimiter='\t')]
    return assemblies

def assembly_filter(ass):
    """Do we want this assembly?"""
    if (ass['version_status'] == "latest" and ass['assembly_level'] == 'Complete Genome') \
        or ass['refseq_category'] == "reference genome":
        name = ass['organism_name'].lower()
        if name.startswith("escherichia") or name.startswith("shigella"):
            return True
    return False

def rsync_assembly(ass):
    uri = ass['ftp_path'].replace("ftp:", "rsync:")
    command = ["rsync", "-av", uri, "."]
    print "Running", ' '.join(command)
    call(command)

def link_assembly(ass):
    """Create a link to assembly using descriptive strain name"""

###################################
# MAIN
###################################

print("Fetching assembly file")
assemblies = load_assemblies()
print len(assemblies), "assemblies"

keepers = [assembly for assembly in assemblies if assembly_filter(assembly)]
print len(keepers), "keepers"

for keeper in keepers:
    rsync_assembly(keeper)
    link_assembly(keeper)
