#!/usr/bin/env python3
import argparse
import glob
import os
import re
import time
from io import StringIO

import pandas as pd
from Bio import SeqIO

# ---------------------------------------------------------------------------- #
# Classes
# ---------------------------------------------------------------------------- #

# Store individual fasta information


class fasta:
    def __init__(self, ogid, seqDict):
        self.ogid = ogid
        self.dict = seqDict


# Store information for many fasta class objects
class fastaCollection:
    def __init__(self):
        self.collection = {}
        self.subset = {}

    def updateCollection(self, classFasta=None):
        if classFasta is not None:
            self.collection.update(classFasta)

    def subsetCollection(self, ids=None):
        if ids is not None:
            keep = {key: value for key, value in self.collection.items() if key in ids}
            self.subset.update(keep)


# ---------------------------------------------------------------------------- #
# Functions
# ---------------------------------------------------------------------------- #

# Create fasta object for each file
def buildOGFasta(path):
    ogid = os.path.splitext(os.path.basename(path))[0]
    seqDict = SeqIO.to_dict(
        SeqIO.parse(path, "fasta"), key_function=lambda rec: rec.description
    )

    # Clean sequence headers (if they contain ':') - e.g. Trinity/seqtk
    clean_dict = {}
    for header in seqDict.keys():
        if ":" in header:
            rec = seqDict[header]
            h = header.split(":", 1)[0]  # Identifer
            s = header.split(" ", 1)[1]  # Species
            new_header = " ".join([h, s])

            # Update seq-record with new naming
            rec.id = h
            rec.name = h
            rec.description = new_header

            clean_dict.update({new_header: rec})
        else:
            clean_dict.update({header: seqDict[header]})

    # Create dict object: { ogid: class fasta }
    return {ogid: fasta(ogid, clean_dict)}


# Read Phyletic Profile table - quick to subset
def readPhyleticProfileOMAGroups(path, col_subset=False):
    ppg = pd.read_table(path, delimiter="\t", skiprows=4, na_values=0)

    # Subset for user selection
    if col_subset is False:
        col_subset = ppg.iloc[:, 1:].columns
        col_subset_group = col_subset.insert(0, "Group")
    else:
        col_subset_group = ["Group"] + col_subset

    # Get complete cases
    idx = ppg.dropna(thresh=len(col_subset), subset=col_subset).index

    # subset for data
    sub = ppg.iloc[idx, :]
    sub = sub.loc[:, col_subset_group]

    # Convert OMA header to OG header (match file names)
    groupids = sub["Group"].tolist()
    groupids = [re.sub("OMA", "", i) for i in groupids]
    groupids = [re.sub("(?<![0-9])0+", "", i) for i in groupids]
    groupids = ["OG" + i for i in groupids]

    # Return list: [ species, [ogids, ...] ]
    return [sub.columns[1:], groupids]


# Get headers from OrthologGroups text file - use OGIDs from PPOMAgroups
def getOrthologousGroups(path, subset):
    # Read lines as list - skip first 4 lines
    with open(path, "r") as file:
        lines = file.readlines()[4:]

    # Remove species from string
    dict_seqs = {}
    for l in lines:

        # Split on tab
        split = l.rstrip("\n").split("\t")

        # Convert OMA header to OG header (match file names)
        ogid = re.sub("OMA", "", split[0])
        ogid = re.sub("(?<![0-9])0+", "", ogid)
        ogid = "OG" + ogid

        # Iterate over present sequences and build a dict: { species: header }
        dict_local = {}
        for i in split[1:]:
            sp = i.split(":", 1)[0]
            hd = i.split(":", 1)[1]
            if (
                ":" in hd
            ):  # Used for seqtk/Trinity subset files (seqtk appends ':1-len(seq)' to header)
                hd = hd.split(":")[0]
            hd = hd + " [" + sp + "]"
            dict_local.update({sp: hd})

        # Update master dict
        dict_seqs.update({ogid: dict_local})

    # Subset data for user specified
    keep = {key: value for key, value in dict_seqs.items() if key in subset[1]}
    return keep


# Write OMA subset to directory
def writeOMAfasta(fastaCollection, dict_orthGroups, phyGroupsObj, outdir):

    # Make output directory
    outdir = os.path.join(outdir, "protein")
    try:
        if not os.path.exists(outdir):
            os.makedirs(outdir)
            # print("Outdir:", outdir)
    except OSError:
        print("Can't create output directory %s" % outdir)

    for key, value in fastaCollection.subset.items():

        headers_dict = dict_orthGroups[key]  # Dict of structure { species: header }
        species = phyGroupsObj[0]  # List of user specified species
        fn = os.path.join(outdir, key + ".fa")  # output file

        # Iterate over species, species specific headers
        headers = []
        for s in species:
            headers.append(headers_dict[s])

        # Iterate over headers + append to seqList -> headers are now species IDs
        seqrecord_lst = []
        for h in headers:
            sp = h.split(" [")[1].replace("]", "")
            seq_rec = value.dict[h]
            seq_rec.description = sp
            seq_rec.id = sp
            seqrecord_lst.append(seq_rec)

        # Write subset of sequences to file
        SeqIO.write(seqrecord_lst, fn, "fasta")


# Get list of Nucleotide fasta files - complementary to peptide sequences
def getNucleotideFiles(path, ext):
    files = glob.glob(path + "/*" + ext)

    file_dict = {}
    for f in files:
        bn = os.path.splitext(os.path.basename(f))[0]

        if "-" in bn:
            bn = bn.split("-")[0]

        seqDict = SeqIO.to_dict(
            SeqIO.parse(f, "fasta"), key_function=lambda rec: rec.description
        )

        # Clean sequence headers (if they contain ':') - e.g. Trinity/seqtk
        clean = {}
        for header in seqDict.keys():
            if ":" in header:
                rec = seqDict[header]
                h = header.split(":", 1)[0]  # Identifer

                # Update seq-record with new naming
                rec.id = h
                rec.name = h
                rec.description = h

                clean.update({h: rec})
            else:
                clean.update({header: seqDict[header]})

        # Create dict object: { species: SeqIO dict }
        file_dict.update({bn: clean})

    return file_dict


# Write nucleotide multi-fasta files that are complementary to OMA peptide
def writeNucFasta(nuc_dict, orthGroup_dict, phyGroup_lst, outdir):

    # Make output directory
    outdir = os.path.join(outdir, "nucleotide")
    try:
        if not os.path.exists(outdir):
            os.makedirs(outdir)
            # print("Outdir:", outdir)
    except OSError:
        print("Can't create output directory %s" % outdir)

    # Species + OG groups of interest
    species = phyGroup_lst[0]
    ogids = phyGroup_lst[1]

    for o in ogids:
        species_dict = orthGroup_dict[
            o
        ]  # Dictionary of { (species: header), (species: header) }
        subset_dict = {
            k: v for k, v in species_dict.items() if k in species
        }  # Subset of 'species_dict' above
        subset_dict.update(
            (k, [v, v.split(" [")[0]]) for k, v in subset_dict.items()
        )  # Remove the '[species]' string from header (put there by OMA)
        fn = os.path.join(outdir, o + ".fa")

        # Extract nucleotide sequences, save to list and write to file
        nuc_seqs = []
        for s, h in subset_dict.items():

            try:
                seq_dict = nuc_dict[s][h[1]]
            except KeyError:
                print(f"\nERROR: Can't find species {s} or header {h[1]}")
                print()

            # # seq_dict = nuc_dict[s][h[1]]
            # seq_dict.description = h[0]
            seq_dict.description = s
            seq_dict.id = s
            nuc_seqs.append(seq_dict)

        SeqIO.write(nuc_seqs, fn, "fasta")


# Convert OrthologousGroups.txt to a CSV file with samples as columns (headers as values) for each OG group
def writeOrthologousGroupsCsv(path, outdir):
    # Read lines as list - skip first 4 lines
    with open(path, "r") as file:
        lines = file.readlines()[4:]

    # Remove species from string
    df_list = []
    for l in lines:

        # Split on tab
        split = l.rstrip("\n").split("\t")

        # Convert OMA id to OG id
        ogid = re.sub("OMA", "", split[0])
        ogid = re.sub("(?<![0-9])0+", "", ogid)
        ogid = "OG" + ogid

        # Get species' as list
        colnames = ["Group"]
        for s in split[1:]:
            colnames.append(s.split(":")[0])

        # clean headers without species
        colvalues = [ogid]
        for h in split[1:]:
            colvalues.append(h.split(":")[1])

        # Create table for each line
        df_list.append(
            pd.read_csv(StringIO(",".join(colvalues)), sep=",", names=colnames)
        )

        # Subset for columns of interest
    df = pd.concat(df_list)
    df.to_csv(os.path.join(outdir, "OrthologousGroups.csv"), index=False)


# Main function
def main():

    # ---------------------------------------------------------------------------- #
    # Parser object
    # ---------------------------------------------------------------------------- #
    parser = argparse.ArgumentParser(
        description="""Filter OMA output for complete Orthologue groups using
        species filters. Write complementary peptide + nucleotide files to fasta
        files."""
    )

    parser.add_argument(
        "-oma", "--omadir", help="Directory path to OMA output directory"
    )
    parser.add_argument(
        "-s",
        "--species",
        default=False,
        help="String of species to subset OMA results by",
        type=str,
    )
    parser.add_argument(
        "-n",
        "--nucleotideDir",
        help="""Directory path to coding sequence fasta files. Must be complementary 
                                                           to protein sequences used in the OMA analysis""",
    )
    parser.add_argument(
        "-e", "--nucleotideExtension", help="File extension for CDS files"
    )
    parser.add_argument(
        "-o",
        "--outdir",
        help="Directory path to output location. 'protein' and 'nucleotide' directories will be created within this location",
    )

    args = parser.parse_args()

    # ---------------------------------------------------------------------------- #
    # OMA fasta files + species into list
    # ---------------------------------------------------------------------------- #

    # Pipeline start
    tic = time.perf_counter()

    # OMA specific files
    oma_fasta = glob.glob(os.path.join(args.omadir, "OrthologousGroupsFasta") + "/*.fa")
    phyleticProfileFile = os.path.join(args.omadir, "PhyleticProfileOMAGroups.txt")
    orthologousGroupsFile = os.path.join(args.omadir, "OrthologousGroups.txt")

    # Species string to list
    if args.species is not False:
        args.species = args.species.split(" ")

    # ---------------------------------------------------------------------------- #
    # Pipeline
    # ---------------------------------------------------------------------------- #
    print("Subsetting phyletic Profile")
    phyleticProfile = readPhyleticProfileOMAGroups(
        path=phyleticProfileFile, col_subset=args.species
    )

    print("Working on Orthologous Groups file")
    orthologousGroups = getOrthologousGroups(
        path=orthologousGroupsFile, subset=phyleticProfile
    )

    print("Importing nucleotide files")
    nucleotideFiles = getNucleotideFiles(
        path=args.nucleotideDir, ext=args.nucleotideExtension
    )

    # Initialise fastaCollection class
    fasta_collection = fastaCollection()

    # Iterate through OMA output fasta files
    for fastaPath in oma_fasta:
        fasta_collection.updateCollection(buildOGFasta(fastaPath))

    # Subset fastaCollection by user specified species filters
    print("Subsetting OMA fasta files for complete orthologues")
    fasta_collection.subsetCollection(ids=orthologousGroups.keys())

    # Write fasta files - peptide and nucleotide (complementary to each other)
    print("Writing peptide and nucleotide files")
    writeOMAfasta(
        fastaCollection=fasta_collection,
        dict_orthGroups=orthologousGroups,
        phyGroupsObj=phyleticProfile,
        outdir=args.outdir,
    )

    writeNucFasta(
        nuc_dict=nucleotideFiles,
        orthGroup_dict=orthologousGroups,
        phyGroup_lst=phyleticProfile,
        outdir=args.outdir,
    )

    # Write orthologous groups CSV file
    writeOrthologousGroupsCsv(path=orthologousGroupsFile, outdir=args.omadir)

    # Pipeline end
    toc = time.perf_counter()
    print(f"Time: {toc-tic:0.4f} seconds")


# ---------------------------------------------------------------------------- #
# Run the pipeline
# ---------------------------------------------------------------------------- #
if __name__ == "__main__":
    main()
