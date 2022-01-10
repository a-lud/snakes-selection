#!/Users/alastair.ludington/miniconda3/bin/python3.8

import subprocess
import argparse
import glob
import logging
import os
import textwrap
from ete3.treeview.main import NodeStyle
import pandas as pd
from functools import reduce
from pathlib import Path

from Bio import SeqIO

import pyham

# Logging information
logging.basicConfig(
    level=logging.INFO, format="%(asctime)s %(name)-12s %(levelname)-8s %(message)s"
)

# --------------------------------------------------------------------------- #
# Functions


def getArgs():
    desc = """\
    # -------------------------------------------------------- #
    #                   Explorer HOG results                   #
    # -------------------------------------------------------- #

    A simple script to obtain genes that have been 'gained' or
    'lost' based on their hierarchical-orthology. That is, by
    using software such as OMA, we can identify sets of genes
    that have descended from a common ancestral sequence in an
    ancestral species. These orthology relationships form a
    graph, with groups within this graph between species forming
    the HOG.

    The hierarchical nature of these orthology groups comes from
    the fact that 'groupings' are based on the taxanomic level.
    A group defined at a recent node in a tree is likely to be
    a sub-group of an older clade.

    This script aims to extract HOGs at a user provided level.
    It will output the following files:
        - Gained: Top-level-hogs shared by all samples at user
                  specified taxonomic level
        - Gained: Fasta files belonging to the HOGs above
        - Lost: Top-level-hogs shared by all samples at a user
                specified taxonomic level
        - Lost: Fasta files beloning to the HOGs above
        - Summary: Summary of gained/lost/duplicated and
                   retained genes within each sample.
    ------------------------------------------------------------
    """

    epi = """\
    Code written by Alastair J. Ludington
    University of Adelaide
    July, 2021
    """

    parser = argparse.ArgumentParser(
        formatter_class=argparse.RawDescriptionHelpFormatter,
        description=textwrap.dedent(desc),
        epilog=textwrap.dedent(epi),
    )

    # Required, positional input file arguments
    parser.add_argument(
        "omadir",
        help="Path to input directory containing OMA results",
        metavar="/path/to/omadir",
    )
    parser.add_argument(
        "outdir", help="Pipeline output directory", metavar="/path/to/outdir"
    )

    parser.add_argument(
        "-pn",
        "--printNodes",
        help="Print the names of all nodes within the tree",
        action="store_true",
    )

    parser.add_argument(
        "-n",
        "--nodesFile",
        help="Text file containing ancestral nodes to get gains/losses",
        type=str,
    )

    parser.add_argument(
        "-tr",
        "--tree",
        help="Custom tree file",
        default=None,
    )

    args = parser.parse_args()
    return args


# Get Tree/OrthoXML file - dict: {id: File}
def getInputFiles(omapath, treepath):
    logging.info("getInputFiles: Getting tree and OrthoXML files")
    files = {}
    if treepath is not None:
        if not Path(treepath).is_file():
            logging.error("getInputFiles: Provided tree file doesn't exist")
            print(f"\nFile: {treepath}")
            exit()
        else:
            files.update({"tree": treepath})
    else:
        files.update({"tree": os.path.join(omapath, "EstimatedSpeciesTree.nwk")})

    files.update({"orthoxml": os.path.join(omapath, "HierarchicalGroups.orthoxml")})
    files.update({"phyleticProfile": os.path.join(omapath, "PhyleticProfileHOGs.txt")})

    return files


# Get HOG fasta files - dict: {HOGID: SeqIO dict}
def getHogFastas(path):
    logging.info("getHogFastas: Reading in HOG Fasta files using SeqIO ")
    files = glob.glob(os.path.join(path, "HOGFasta", "*.fa"))

    file_dict = {}
    for f in files:

        bn = os.path.splitext(os.path.basename(f))[0].replace("HOG", "")

        seqDict = SeqIO.to_dict(
            SeqIO.parse(f, "fasta"), key_function=lambda rec: rec.description
        )

        # Create dict object: { species: SeqIO dict }
        file_dict.update({bn: seqDict})

    return file_dict


def getHam(tree, xml):
    logging.info("getHam: Reading in OrthoXML as HAM object")
    ham = pyham.Ham(tree, xml, use_internal_name=False)
    return ham


# Print Nodes to screen
def printNodes(ham, outdir):
    internal_nodes = os.path.join(outdir, "nodes.txt")
    logging.info("printNodes: Printing tree and nodes")
    print("\nTree Visualisaion:", ham.taxonomy.tree)
    print("\nInternal Nodes to Choose From:")

    try:
        os.remove(internal_nodes)
    except OSError:
        pass

    with open(internal_nodes, "w") as file:
        for i in ham.taxonomy.internal_nodes:
            print("\t- {}".format(i.name))
            print(i.name + "\n", file=file)

    exit()


# Read nodes file and return list of ancestral nodes
def getNodes(outdir):
    node_file = os.path.join(outdir, "nodes.txt")
    with open(node_file, "r") as file:
        lines = list(filter(None, (line.rstrip() for line in file)))

    return lines


# Return samples not apart of the selected node
def buildGroups(node, ham):
    logging.info(
        "buildGroups: Splitting tree into samples within and outside of selected node"
    )
    selected = node.split("/")
    leaves = ham.taxonomy.tree.get_leaf_names()

    dif = list(set(leaves) - set(selected))

    return {"samples": selected, "background": dif}


# Get the 'root' name
def getRootName(ham):

    return ham.taxonomy.tree.get_tree_root().name


# Get ancestral genome objects - group 1 and top level to compare to
def getAncestralGenomes(ham, node, root):
    anc_samples = ham.get_ancestral_genome_by_name(node)
    anc_root = ham.get_ancestral_genome_by_name(root)

    return {"anc_sample": anc_samples, "anc_root": anc_root}


# Iterate over each snake and get: gain/loss/dup/retained
def doVerticalAnalysis(ham, groups_dict, anc_root, outdir):
    groups_out = {}  # Dict for all output to be stored

    # Iterate over groups
    for group, lst_snakes in groups_dict.items():
        logging.info(f"doVerticalAnalysis: Current group - {group}")
        current_snake = {}

        # Iterate over the snakes within the group
        for snake in lst_snakes:
            logging.info(f"doVerticalAnalysis: Current snake - {snake}")

            # Get ancestral genome object for current snake
            anc_snake = ham.get_extant_genome_by_name(snake)

            # Vertical analysis - Gene evo. from root to current snake
            vertical = ham.compare_genomes_vertically(anc_snake, anc_root)

            retained = vertical.get_retained()
            duplicated = vertical.get_duplicated()
            gained = vertical.get_gained()
            lost = vertical.get_lost()

            # Count number of HOGs/genes (?) belonging to each
            retained_len = len(retained)
            duplicated_len = len(duplicated)
            gained_len = len(gained)
            lost_len = len(lost)

            # Append current snakes results to dictionary
            current_snake.update(
                {
                    snake: {
                        "vertical": vertical,
                        "retained": [retained, retained_len],
                        "duplicated": [duplicated, duplicated_len],
                        "gained": [gained, gained_len],
                        "lost": [lost, lost_len],
                    }
                }
            )
        groups_out.update({group: current_snake})

    # Write summary file
    summary = os.path.join(outdir, "summary-by-sample.csv")
    try:
        os.remove(summary)
    except OSError:
        pass

    # Column names
    with open(file=summary, mode="w") as file:
        print("Sample,Retained,Duplicated,Gained,Lost", file=file)

    # Iterate over dictionary and write file information
    for group, results in groups_out.items():
        for snake, data in results.items():
            with open(summary, mode="a") as file:
                rl = data["retained"][1]
                dl = data["duplicated"][1]
                gl = data["gained"][1]
                ll = data["lost"][1]
                print(f"{snake},{rl},{dl},{gl},{ll}", file=file)

    return groups_out


# Get the gained and lost genes within each snake - BY SNAKE CURRENTLY
def getGainLostByGroup(verticaldict, root):
    # Get outgroup
    outgroup = root.split("/")[-1]

    # Dictionaries for each group
    gained_samples = {}
    gained_background = {}
    lost_samples = {}
    lost_background = {}

    logging.info("getGainLostByGroup: Getting gained/lost HOGs by group")

    # Iterate over verticaldic's groups and samples
    for group, dict_snakes in verticaldict.items():

        # iterate over each snake
        for snake, analyses in dict_snakes.items():
            gained = analyses["gained"][0]
            lost = analyses["lost"][0]

            if group == "samples":
                gained_samples.update({snake: gained})
                lost_samples.update({snake: lost})
            else:
                gained_background.update({snake: gained})
                lost_background.update({snake: lost})

    # Remove outgroup leaf from lost - being top level,
    # it won't be missing any genes.
    if outgroup in lost_samples.keys():
        lost_samples.pop(outgroup)
    else:
        lost_background.pop(outgroup)

    return {
        "gained": {"samples": gained_samples, "background": gained_background},
        "lost": {"samples": lost_samples, "background": lost_background},
    }


def parsePhyleticProfile(group_dict, phyleticProfile, scriptDir, outdir):
    rscript = os.path.join(scriptDir, "parsePhyleticProfile.R")

    # R-script command
    cmd = [
        rscript,
        "-f",
        str(phyleticProfile),
        "-s",
        str(" ".join(group_dict["samples"])),
        "-b",
        str(" ".join(group_dict["background"])),
        "-o",
        str(outdir),
    ]

    # Run the R-script
    try:
        subprocess.run(
            cmd,
            check=True,
            stderr=subprocess.PIPE,
            universal_newlines=True,
        )
        logging.info(
            "parsePhyleticProfile: Obtaining gained/lost HOGs for sample group"
        )
    except subprocess.CalledProcessError as error:
        logging.error("parsePhyleticProfile: Issue parsing HOG Phyletic-Profile")
        print(
            f"\nCommand: {error.cmd}\nExit code: {error.returncode}\nSTDERR: {error.stderr}"
        )
        exit()

    # Import the HOG ids as lists
    with open(os.path.join(outdir, "gained-samples-id.txt")) as file:
        gained_hogs = file.readlines()

    with open(os.path.join(outdir, "lost-samples-id.txt")) as file:
        lost_hogs = file.readlines()

    with open(os.path.join(outdir, "gained-samples-threshold-id.txt")) as file:
        gain_threshold = file.readlines()

    # Files containing reference samples the gain/lost genes need to be
    # validated against
    refs_gain = str(" ".join(group_dict["background"]))
    refs_lost = str(" ".join(group_dict["samples"]))

    with open(os.path.join(outdir, "refs-validate-gains.txt"), mode="w") as file:
        print(refs_gain, file=file)

    with open(os.path.join(outdir, "refs-validate-lost.txt"), mode="w") as file:
        print(refs_lost, file=file)

    return {
        "gained": gained_hogs,
        "lost": lost_hogs,
        "gain_threshold": gain_threshold,
    }


# Copy relevant HOGs to ouput locations
def copyFasta(
    gainlost,
    hogfastas,
    lostfastadir,
    gainedfastadir,
    gainedthresholdfastadir,
):

    logging.info("copyFasta: Writing lost HOGs to file")
    for hog in gainlost["lost"]:
        hog = hog.rstrip()
        filename = os.path.join(lostfastadir, hog + ".fa")
        seqDict = hogfastas[hog[3:]]
        seqrecords_list = []
        for headers in seqDict.keys():
            seqrecords_list.append(seqDict[headers])

        SeqIO.write(seqrecords_list, filename, "fasta")

    for hog in gainlost["gained"]:
        hog = hog.rstrip()
        filename = os.path.join(gainedfastadir, hog + ".fa")
        seqDict = hogfastas[hog[3:]]
        seqrecords_list = []
        for headers in seqDict.keys():
            seqrecords_list.append(seqDict[headers])

        SeqIO.write(seqrecords_list, filename, "fasta")

    for hog in gainlost["gain_threshold"]:
        hog = hog.rstrip()
        filename = os.path.join(gainedthresholdfastadir, hog + ".fa")
        seqDict = hogfastas[hog[3:]]
        seqrecords_list = []
        for headers in seqDict.keys():
            seqrecords_list.append(seqDict[headers])

        SeqIO.write(seqrecords_list, filename, "fasta")


# TODO: Create output directories and files once I know what I want
def createOutputs(arguments, nodeNum=None, node=None):
    logging.info(
        "createOutputs: Create output directory and filepaths",
    )

    if all(n is not None for n in [nodeNum, node]):
        # Make the output directory - if can't do this, exit
        try:
            outdir = os.path.join(arguments.outdir, "ancestral-node-" + nodeNum)
            Path(outdir).mkdir(parents=True, exist_ok=True)
        except OSError:
            logging.error("Couldn't create output directory")
            exit()

        # Output directories
        lost_samples_fasta_dir = os.path.join(outdir, "samples", "lost_fasta")
        ganied_samples_fasta_dir = os.path.join(outdir, "samples", "gained_fasta")
        gained_threshold_samples_fasta_dir = os.path.join(
            outdir, "samples", "gained_threshold_fasta"
        )
        Path(lost_samples_fasta_dir).mkdir(parents=True, exist_ok=True)
        Path(ganied_samples_fasta_dir).mkdir(parents=True, exist_ok=True)
        Path(gained_threshold_samples_fasta_dir).mkdir(parents=True, exist_ok=True)

        # Output files
        gained_hogs = os.path.join(outdir, "gained-samples-id.txt")
        lost_hogs = os.path.join(outdir, "lost-samples-id.txt")
        node_file = os.path.join(outdir, "analysis-node.txt")
        # gained_samples_singleton = os.path.join(outdir, "gained-samples-singleton.csv")

        # Remove them if they exist already
        try:
            os.remove(gained_hogs)
            os.remove(lost_hogs)
            # os.remove(gained_samples_singleton)
        except OSError:
            pass

        # with open(file=gained_samples_singleton, mode="w") as file:
        #     print("Sample,id,protID", file=file)

        with open(file=node_file, mode="w") as file:
            print(node, file=file)

        # Return a dict
        outputs = {
            "outdir": outdir,
            "lost_samples_fasta_dir": lost_samples_fasta_dir,
            # "gained_samples_singleton": gained_samples_singleton,
            "gained_samples_fasta_dir": ganied_samples_fasta_dir,
            "gained_threshold_samples_fasta_dir": gained_threshold_samples_fasta_dir,
        }

        return outputs
    else:
        # Make the output directory - if can't do this, exit
        try:
            outdir = arguments.outdir
            Path(outdir).mkdir(parents=True, exist_ok=True)
        except OSError:
            logging.error("Couldn't create output directory")
            exit()


if __name__ == "__main__":

    # Get arguments
    args = getArgs()

    # Script directory
    script_dir = os.path.abspath(os.path.dirname(__file__))

    # Create output directory + output files
    out_files = createOutputs(args)

    # Get the paths to required files for analysis
    oma_files = getInputFiles(args.omadir, args.tree)

    # Create the HAM object
    ham = getHam(oma_files["tree"], oma_files["orthoxml"])

    # Print the tree +  ancestral node names and exit if requested
    if args.printNodes:
        printNodes(ham, args.outdir)

    # Get ancestral nodes to test
    ancestral_node_list = getNodes(args.outdir)

    # Can take a while - do after checking print-nodes argument
    hog_fastas = getHogFastas(args.omadir)

    # Top level ancestralGenome name (I take the root)
    root = getRootName(ham)

    # Iterate over ancestral nodes generating output
    for count, node in enumerate(ancestral_node_list):
        # Create outputs for each node
        out_files = createOutputs(args, str(count), node)

        # Create groups: samples within selected ancNode and outside it
        groups_dict = buildGroups(node, ham)

        # Get ancestral genome objects for selected ancNode and root
        anc_dict = getAncestralGenomes(ham, node, root)

        # Get the gain/loss profiles for the sample-group of interest
        gained_lost_hogs = parsePhyleticProfile(
            groups_dict,
            oma_files["phyleticProfile"],
            script_dir,
            out_files["outdir"],
        )

        copyFasta(
            gained_lost_hogs,
            hog_fastas,
            out_files["lost_samples_fasta_dir"],
            out_files["gained_samples_fasta_dir"],
            out_files["gained_threshold_samples_fasta_dir"],
        )

        # Veritcal analysis: Snakes to ancestral genome
        vertical = doVerticalAnalysis(
            ham,
            groups_dict,
            anc_dict["anc_root"],
            args.outdir,
        )

    # Get dictionaries of gained/lost genes by snake
    # gained_lost = getGainLostByGroup(vertical, root)
