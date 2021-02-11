#!/usr/bin/env python3
"""
snapp_prep.py was ported to python by Carlos Alonso Maya-Lastra. 
snapp_prep.rb originally was written by Michael Matschiner (michaelmatschiner@mac.com), 2018-05-21 
https://github.com/mmatschiner/snapp_prep

Original Snapp_prep.rb description:
This script prepares XML format input files for the software SNAPP (http://beast2.org/snapp/),
given a phylip format SNP matrix, a table linking species IDs and specimen IDs, and a file
specifying age constraints. Optionally, a starting tree can be provided (recommended) and the
number of MCMC iterations used by SNAPP can be specified. Files prepared by this script use
a particular combination of priors and operators optimized for SNAPP analyses that aim to
estimate species divergence times (note that all population sizes are linked in these analyses).
Detailed annotation explaining the choice of priors and operators will be included in the
XML file (unless turned off with option '-n').

This script can be executed from command line as:
python3 snapp_prep.py -p example.phylip -s example.tre -x snapp.xml

Or as a module:
import snapp_prep
snapp_prep.run(phylip="example.phylip", tree="example.tre", xml="snapp.xml")

"""

# Load required libraries.
import argparse
import sys
import random 
import re




#Main function that prepares the XML for SNAPP
def run(phylip=None, vcf=None, table="example.spc.txt", constraints="example.con.txt", tree=None, length=500000, weight=1.0, max_snps=None, min_dist=None, transversions=False, transitions=False, xml="snapp.xml", out="snapp", no_annotation=False):
    """
    Prepares the XML file to run SNAPP
    """

    # Feedback.
    print ("\nsnapp_prep.py (a port of snapp_prep.rb)\n----------------------------------------------------------------------------------------")

    #Define Unique function without arbitrary reorder
    def unique(sequence):
        seen = set()
        return [x for x in sequence if not (x in seen or seen.add(x))]

    # Make sure that input is provided in either phylip or vcf format.
    if phylip is None and vcf is None:
        print("ERROR: An input file must be provided, either in phylip format with option '-p' or in vcf format with option '-v'!")
        sys.exit(1)
    elif phylip and vcf:
        print("ERROR: Only one of the two options '-p' and '-v' can be used!")
        sys.exit(1)

    # Make sure that the -r and -i options are not used jointly.
    if transversions and transitions:
        print("ERROR: Only one of the two options '-r' and '-i' can be used!")
        sys.exit(1)


    # Make sure that the specified weight is not negative.
    if weight < 0:
        print("ERROR: The specified relative weight of the topology operator must not be negative!")
        sys.exit(1)

    # Initiate a warn string and counts for excluded sites.
    warn_string = ""
    number_of_excluded_sites_missing = 0
    number_of_excluded_sites_monomorphic = 0
    number_of_excluded_sites_triallelic = 0
    number_of_excluded_sites_tetraallelic = 0
    number_of_excluded_sites_indel = 0
    number_of_excluded_sites_transition = 0
    number_of_excluded_sites_transversion = 0
    number_of_sites_with_half_call = 0

    # Initiate arrays for specimen ids and sequences.
    specimen_ids = []
    seqs = []

    # Define various characters.
    binary_chars = ["0","1","2"]
    nucleotide_chars = ["A","C","G","T","R","Y","S","W","K","M"]
    ambiguous_chars = ["R","Y","S","W","K","M"]
    missing_chars = ["-","?","N"]


    # Read the phylip file.
    if phylip:
        with open(phylip) as phylip_file: 
            phylip_lines = phylip_file.readlines()
            for l in phylip_lines[1:]:
                if l is not "":
                    specimen_ids.append(l.split()[0])
                    seqs.append(l.split()[1].upper())

        # Recognize the sequence format (nucleotides or binary).
        sequence_format_is_nucleotide = False 
        sequence_format_is_binary =  False 
        unique_seq_chars = unique("".join(seqs)) #
        if set(unique_seq_chars).intersection(set(nucleotide_chars)):
            sequence_format_is_nucleotide = True
        if set(unique_seq_chars).intersection(set(binary_chars)):
            sequence_format_is_binary = True
        if sequence_format_is_binary == False and sequence_format_is_nucleotide == False:
            print("ERROR: Sequence format could not be recognized as either 'nucleotide' or 'binary'!")
            sys.exit(1)
        if sequence_format_is_binary == True and sequence_format_is_nucleotide == True:
            print("ERROR: Sequence format was recognized as a mixture of 'nucleotide' and 'binary'!")
            sys.exit(1)
        if sequence_format_is_binary == True:
            print("SORRY: Binary data is not supported by this port.")
            sys.exit(1)


    # Alternatively, read the vcf file.
    ##TODO not implemented in this port


    # If necessary, translate the sequences into SNAPP's "0", "1", "2" code, where "1" is heterozygous.
    binary_seqs = ["" for i in seqs]


    # if sequence_format_is_binary: .......TODO.....

    if sequence_format_is_nucleotide:
        for pos in range(len(seqs[0])):
            # Collect all bases at this position.
            bases_at_this_pos = []
            for s in seqs:
                if s[pos] == "A":
                    bases_at_this_pos.append("A")
                    bases_at_this_pos.append("A")
                elif s[pos] == "C":
                    bases_at_this_pos.append("C")
                    bases_at_this_pos.append("C")
                elif s[pos] == "G":
                    bases_at_this_pos.append("G")
                    bases_at_this_pos.append("G")
                elif s[pos] == "T":
                    bases_at_this_pos.append("T")
                    bases_at_this_pos.append("T")
                elif s[pos] == "R":
                    bases_at_this_pos.append("A")
                    bases_at_this_pos.append("G")
                elif s[pos] == "Y":
                    bases_at_this_pos.append("C")
                    bases_at_this_pos.append("T")
                elif s[pos] == "S":
                    bases_at_this_pos.append("G")
                    bases_at_this_pos.append("C")
                elif s[pos] == "W":
                    bases_at_this_pos.append("A")
                    bases_at_this_pos.append("T")
                elif s[pos] == "K":
                    bases_at_this_pos.append("G")
                    bases_at_this_pos.append("T")
                elif s[pos] == "M":
                    bases_at_this_pos.append("A")
                    bases_at_this_pos.append("C")
                elif s[pos] not in missing_chars:
                    print(f"ERROR: Found unexpected base at position {pos+1}: {s[pos]}!")
                    sys.exit(1)
            
            uniq_bases_at_this_pos = sorted(unique(bases_at_this_pos)) #
            # Issue a warning if non-bi-allelic sites are excluded.
            if len(uniq_bases_at_this_pos) == 2:
                # Determine if this is a transition or transversion site.
                transversion_site = False
                if uniq_bases_at_this_pos == ["A","C"]:
                    transversion_site = True
                elif uniq_bases_at_this_pos == ["A","G"]:
                    transversion_site = False
                elif uniq_bases_at_this_pos == ["A","T"]:
                    transversion_site = True
                elif uniq_bases_at_this_pos == ["C","G"]:
                    transversion_site = True
                elif uniq_bases_at_this_pos == ["C","T"]:
                    transversion_site = False
                elif uniq_bases_at_this_pos == ["G","T"]:
                    transversion_site = True
                else:
                    print (f"ERROR: Unexpected combination of unique bases at position {pos+1}: {uniq_bases_at_this_pos[0]}, {uniq_bases_at_this_pos[1]}")
                    sys.exit(1)
                
                transition_site = True
                if transversion_site == True: 
                    transition_site = False 
                if transversions == True and transversion_site == False: # If the site is a transition and only transversions are allowed.
                    number_of_excluded_sites_transition += 1
                elif transitions == True and transition_site == False: # If the site is a transversion and only transitions are allowed.
                    number_of_excluded_sites_transversion += 1
                else:
                    # Randomly define what's "0" and "2".
                    random.shuffle(uniq_bases_at_this_pos)
                    for x in range(len(seqs)):
                        if seqs[x][pos] == uniq_bases_at_this_pos[0]:
                            binary_seqs[x] += "0"
                        elif seqs[x][pos] == uniq_bases_at_this_pos[1]:
                            binary_seqs[x] += "2"
                        elif seqs[x][pos] in missing_chars:
                            binary_seqs[x] += "-"
                        elif seqs[x][pos] in ambiguous_chars:
                            binary_seqs[x] += "1"
                        else:
                            print (f"ERROR: Found unexpected base at position {pos+1}: {seqs[x][pos]}!")
                            sys.exit(1)
                    
                    #if options[:vcf] != nil #Not implemented in this port yet
            elif len(uniq_bases_at_this_pos) == 0:
                number_of_excluded_sites_missing += 1
            elif len(uniq_bases_at_this_pos) == 1:
                number_of_excluded_sites_monomorphic += 1
            elif len(uniq_bases_at_this_pos) == 3:
                number_of_excluded_sites_triallelic += 1
            elif len(uniq_bases_at_this_pos) == 4:
                number_of_excluded_sites_tetraallelic += 1
            else:
                print (f"ERROR: Found unexpected number of alleles at position {pos+1}!")
                sys.exit(1)


    #if options[:vcf] != nil


    # Read the file with a table linking species and specimens.

    with open(table) as table_file: 
        table_lines = table_file.readlines()
        table_species = []
        table_specimens = []
        for l in table_lines:
            if l is not "":
                line_ary = l.split()
                header_line = False
                if line_ary[0].lower() == "species" and line_ary[1].lower() in ["specimen", "specimens", "sample", "samples"]:
                    header_line = True
                if header_line == False:
                    table_species.append(line_ary[0])
                    table_specimens.append(line_ary[1])

    # Make sure that the arrays table_specimens and specimen_ids are identical when sorted.
    if sorted(table_specimens) != sorted(specimen_ids):
        print (f"ERROR: The specimens listed in file {table} and those included in the input file are not identical! Check the following names: {', '.join(set(table_specimens).difference(set(specimen_ids)))}") #
        sys.exit(1)


    # Remove sites at which one or more species have only missing data; these could not be used by SNAPP anyway.
    binary_seqs_filtered = ["" for i in binary_seqs]
    # if options[:vcf] != nil #Not implemented in this port

    for pos in range(len(binary_seqs[0])):
        one_or_more_species_have_only_missing_data_at_this_pos = False
        for spc in unique(table_species):
            specimens_for_this_species = []
            for x in range(len(table_specimens)):
                if table_species[x] == spc:
                    specimens_for_this_species.append(table_specimens[x])

            alleles_for_this_species_at_this_pos = []
            for x in range(len(specimen_ids)):
                if specimen_ids[x] in specimens_for_this_species:
                    alleles_for_this_species_at_this_pos.append(binary_seqs[x][pos])
            if unique(alleles_for_this_species_at_this_pos) == ["-"]:
                one_or_more_species_have_only_missing_data_at_this_pos =  True

        # Set all alleles at this position to nil if one species had only missing data.
        if one_or_more_species_have_only_missing_data_at_this_pos:
            number_of_excluded_sites_missing += 1
        else:
            for x in range(len(binary_seqs)):
                binary_seqs_filtered[x] += binary_seqs[x][pos]

        # 	# if options[:vcf] != nil #Not implemented in this port
    binary_seqs = binary_seqs_filtered

        # # if options[:vcf] != nil #Not implemented in this port

    # Make sure that information on linkage groups and positions on these have the same length.
    # if options[:vcf] != nil #Not implemented in this port

    # If a minimum distance between SNPs has been set, thin the data set according to this number.
    # if options[:vcf] != nil and options[:min_dist] != nil and options[:min_dist] > 1
    ##vcf is not implemented in this port



    # If a maximum number of SNPs has been set, reduce the data set to this number.
    number_of_sites_before_excluding_due_to_max = len(binary_seqs[0])
    number_of_excluded_sites_due_to_max = 0

    if max_snps != None:
        if max_snps < 1:
            print ( "ERROR: The specified maximum number of SNPs must be positive!")
            sys.exit(1)
        elif max_snps < len(binary_seqs[0]):
            seq_indices = []
            for x in range(len(binary_seqs[0])):
                seq_indices.append(x)
            selected_seq_indices = sorted(random.sample(seq_indices, max_snps))
            binary_seqs_red = []
            for s in binary_seqs:
                binary_seq_red = ""
                for i in selected_seq_indices:
                    binary_seq_red += s[i]
                binary_seqs_red.append(binary_seq_red)
            binary_seqs = binary_seqs_red
            number_of_excluded_sites_due_to_max = number_of_sites_before_excluding_due_to_max - max_snps
        else:
            warn_string += f"WARNING: The maximum number of SNPs has been set to {max_snps}, which is greater\n"
            warn_string += f"   than the number of bi-allelic SNPs with sufficient information ({len(binary_seqs[0])}) for SNAPP.\n"



    # Compose the warn string if necessary.
    if number_of_sites_with_half_call > 0:
        warn_string += f"WARNING: Found {number_of_sites_with_half_call} site"
        if number_of_sites_with_half_call > 1: warn_string += "s" 
        warn_string += " with genotypes that were half missing. These genotypes were ignored.\n"
    if number_of_excluded_sites_missing > 0:
        warn_string += f"WARNING: Excluded {number_of_excluded_sites_missing} site"
        if number_of_excluded_sites_missing > 1: warn_string += "s" 
        warn_string += " with only missing data in one or more species.\n"
    if number_of_excluded_sites_monomorphic > 0:
        warn_string += f"WARNING: Excluded {number_of_excluded_sites_monomorphic} monomorphic site"
        if number_of_excluded_sites_monomorphic > 1: warn_string += "s" 
        warn_string += ".\n"
    if number_of_excluded_sites_transition > 0:
        warn_string += f"WARNING: Excluded {number_of_excluded_sites_transition} transition site"
        if number_of_excluded_sites_transition > 1: warn_string += "s" 
        warn_string += ".\n"
    if number_of_excluded_sites_transversion > 0:
        warn_string += f"WARNING: Excluded {number_of_excluded_sites_transversion} transversion site"
        if number_of_excluded_sites_transversion > 1: warn_string += "s" 
        warn_string += ".\n"
    if number_of_excluded_sites_triallelic > 0:
        warn_string += f"WARNING: Excluded {number_of_excluded_sites_triallelic} tri-allelic site"
        if number_of_excluded_sites_triallelic > 1: warn_string += "s" 
        warn_string += ".\n"
    if number_of_excluded_sites_tetraallelic > 0:
        warn_string += f"WARNING: Excluded {number_of_excluded_sites_tetraallelic} tetra-allelic site"
        if number_of_excluded_sites_tetraallelic > 1: warn_string += "s" 
        warn_string += ".\n"

    # If there were any warning, print them.
    if warn_string:
        # warn_string += "\n"
        print(warn_string)


    # Compose the info string if necessary.
    info_string = ""
    if vcf != None and option.min_dist != None and min_dist > 1:
        info_string += f"INFO: Removed {number_of_excluded_sites_due_to_min_dist} bi-allelic sites due to specified minimum distance between sites of {min_dist} bp.\n"
    if max_snps != None:
        info_string += f"INFO: Removed {number_of_excluded_sites_due_to_max} bi-allelic sites due to specified maximum number of {max_snps} sites.\n"
    if transversions:
        info_string = f"INFO: Retained {len(binary_seqs[0])} bi-allelic transversion sites.\n"
    elif transitions:
        info_string = f"INFO: Retained {len(binary_seqs[0])} bi-allelic transition sites.\n"
    else:
        info_string = f"INFO: Retained {len(binary_seqs[0])} bi-allelic sites.\n"

    # If there were any infos, print them.
    if info_string:
        # info_string += "\n"
        print(info_string)


    # Read the file with age constraint information.
    with open(constraints) as constraint_file: 
        constraint_lines = constraint_file.readlines()
        constraint_strings = []
        cladeage_constraints_used = False
        for l in constraint_lines:
            if l[0:6].lower() in ["normal","lognor","unifor","cladea","monoph"]:
                if l[0:8].lower() == "cladeage":
                    cladeage_constraints_used = True
                constraint_strings.append(l)
    if len(constraint_strings) == 0:
        print (f"WARNING: No age constraints could be found in file {constraints}. You will have to manually modify the XML file to include age constraints.")
    if cladeage_constraints_used:
        info_string = f"INFO: CladeAge constraints are specified in file {constraints}.\n"
        info_string += "    To use these in SNAPP, make sure that the CladeAge package for BEAST2 is installed.\n"
        info_string += "    Installation instructions can be found at http://evoinformatics.eu/cladeage.pdf.\n"
        # info_string += "\n"
        print (info_string)

    # Read the starting tree input file if specified.
    if tree:
        with open(tree) as tree_file: 
            tree_lines = tree_file.readlines()
            tree_string = re.sub("\[.+?\]", "", tree_lines[0].strip().strip(";"))
    else:
        warn_string = "WARNING: As no starting tree has been specified, a random starting tree will be used by\n"
        warn_string += "    SNAPP. If the random starting tree is in conflict with specified constraints, SNAPP\n"
        warn_string += "    may not be able to find a suitable state to initiate the MCMC chain. This will lead\n"
        warn_string += "    to an error message such as 'Could not find a proper state to initialise'. If such\n"
        warn_string += "    a problem is encountered, it can be solved by providing a starting tree in which\n"
        warn_string += "    the ages of constrained clades agree with the constraints placed on these clades.\n"
        # warn_string += "\n"
        print (warn_string)


    # Set run parameters.
    store_frequency = length/2000
    screen_frequency = 50
    snapp_log_file_name = f"{out}.log"
    snapp_trees_file_name = f"{out}.trees"


    # Prepare SNAPP input string.
    snapp_string = ""
    snapp_string += "<?xml version=\"1.0\" encoding=\"UTF-8\" standalone=\"no\"?>\n"
    snapp_string += "<beast beautitemplate='SNAPP' beautistatus='' namespace=\"beast.core:beast.evolution.alignment:beast.evolution.tree.coalescent:beast.core.util:beast.evolution.nuc:beast.evolution.operators:beast.evolution.sitemodel:beast.evolution.substitutionmodel:beast.evolution.likelihood\" version=\"2.0\">\n"
    snapp_string += "\n"

    snapp_string += "<!-- Data -->\n"
    if not no_annotation:
        snapp_string += "<!--\n"
        if sequence_format_is_binary:
            snapp_string += f"The SNP data matrix, taken from file {phylip}.\n"
        else:
            snapp_string += f"The SNP data matrix, converted to binary format from file {phylip}.\n"
        snapp_string += "-->\n"
    snapp_string += "<data id=\"snps\" dataType=\"integer\" name=\"rawdata\">\n"
    for x in range(len(specimen_ids)):
        snapp_string += f"    <sequence id=\"seq_{specimen_ids[x]}\" taxon=\"{specimen_ids[x]}\" totalcount=\"3\" value=\"{binary_seqs[x]}\"/>\n"
    snapp_string += "</data>\n"
    snapp_string += "\n"

    snapp_string += "<!-- Maps -->\n"
    snapp_string += "<map name=\"Uniform\" >beast.math.distributions.Uniform</map>\n"
    snapp_string += "<map name=\"Exponential\" >beast.math.distributions.Exponential</map>\n"
    snapp_string += "<map name=\"LogNormal\" >beast.math.distributions.LogNormalDistributionModel</map>\n"
    snapp_string += "<map name=\"Normal\" >beast.math.distributions.Normal</map>\n"
    snapp_string += "<map name=\"Gamma\" >beast.math.distributions.Gamma</map>\n"
    snapp_string += "<map name=\"OneOnX\" >beast.math.distributions.OneOnX</map>\n"
    snapp_string += "<map name=\"prior\" >beast.math.distributions.Prior</map>\n"
    snapp_string += "\n"
    snapp_string += f"<run id=\"mcmc\" spec=\"MCMC\" chainLength=\"{length}\" storeEvery=\"{int(store_frequency)}\">\n"
    snapp_string += "\n"

    snapp_string += "    <!-- State -->\n"
    snapp_string += f"    <state id=\"state\" storeEvery=\"{int(store_frequency)}\">\n"
    if tree:
        snapp_string += f"        <stateNode id=\"tree\" spec=\"beast.util.TreeParser\" IsLabelledNewick=\"true\" nodetype=\"snap.NodeData\" newick=\"{tree_string};\">\n"
    else:
        snapp_string += "        <stateNode id=\"tree\" spec=\"beast.util.ClusterTree\" clusterType=\"upgma\" nodetype=\"snap.NodeData\">\n"
    snapp_string += "            <taxa id=\"data\" spec=\"snap.Data\" dataType=\"integerdata\">\n"
    snapp_string += "                <rawdata idref=\"snps\"/>\n"
    for s in unique(table_species):
        snapp_string += f"                <taxonset id=\"{s}\" spec=\"TaxonSet\">\n"
        for x in range(len(table_species)):
            if table_species[x] == s:
                snapp_string += f"                    <taxon id=\"{table_specimens[x]}\" spec=\"Taxon\"/>\n"
        snapp_string += "                </taxonset>\n"
    snapp_string += "            </taxa>\n"
    snapp_string += "        </stateNode>\n"
    snapp_string += "        <!-- Parameter starting values -->\n"
    snapp_string += "        <parameter id=\"lambda\" lower=\"0.0\" name=\"stateNode\">0.1</parameter>\n"
    snapp_string += "        <parameter id=\"coalescenceRate\" lower=\"0.0\" name=\"stateNode\">200.0</parameter>\n"
    snapp_string += "        <parameter id=\"clockRate\" lower=\"0.0\" name=\"stateNode\">0.01</parameter>\n"
    snapp_string += "    </state>\n"
    snapp_string += "\n"
    snapp_string += "    <!-- Posterior -->\n"
    snapp_string += "    <distribution id=\"posterior\" spec=\"util.CompoundDistribution\">\n"
    snapp_string += "        <distribution id=\"prior\" spec=\"util.CompoundDistribution\">\n"
    snapp_string += "\n"

    snapp_string += "            <!-- Divergence age priors -->\n"
    constraint_count = 0
    for c in constraint_strings:
        constraint_count += 1
        constraint_ary = c.split()
        if len(constraint_ary) != 3:
            print ("ERROR: Expected three character strings per line for each constraint specification, but found")
            print (f"    '{c.strip()}',")
            sys.exit(1)
        constraint_distribution = constraint_ary[0].lower()
        constraint_placement = constraint_ary[1].lower()
        constraint_clade = constraint_ary[2]
        constraint_clade_ary = constraint_clade.split(",")
        for s in constraint_clade_ary:
            if s not in table_species:
                print (f"ERROR: A constraint has been specified for a clade that includes species {s}; however,")
                print (f"    this species could not be found in file {table}!")
                sys.exit(1)
        if constraint_placement not in ["stem","crown","na"]:
            print ("ERROR: Expected 'stem', 'crown', or 'NA' (only for monophyly constraints without calibration)")
            print (f"    but found '{constraint_placement}' as the second character string in '{c}'!")
            sys.exit(1)
        if "(" not in constraint_distribution and constraint_distribution != "monophyletic":
            print ("ERROR: Expected parameters in parentheses as part of the first character string in ")
            print (f"    '{c.strip()}',")
            print (f"    but found '{constraint_distribution}'!")
            sys.exit(1)
        if constraint_distribution == "monophyletic":
            constraint_type = "monophyletic"
        else:
            constraint_type = constraint_distribution.split("(")[0].strip()
        #This detection is not working in this port neither the original script, to this point any line with a different type is not passed to constraint_strings (check when constraint file is loaded)
        if constraint_type not in ["normal","lognormal","uniform","cladeage","monophyletic"]:
            print ("ERROR: Expected 'normal', 'lognormal', 'uniform', 'cladeage', or 'monophyletic' as part")
            print (f"    of the first character string in '{c}' but found '{constraint_type}'!")
            sys.exit(1)
        if constraint_type != "monophyletic":
            constraint_parameters = constraint_distribution.split("(")[1].strip().replace(")","").split(",")
    ##constraint_parameter in mnophyletic has last memory check for possible bugs

        constraint_parameters_size = len(constraint_parameters)
        if constraint_type == "normal":
            if constraint_parameters_size != 3:
                print (f"ERROR: Expected 3 parameters for normal distribution, but found {constraint_parameters_size}!")
                sys.exit(1)
        elif constraint_type == "lognormal":
            if constraint_parameters_size != 3:
                print (f"ERROR: Expected 3 parameters for lognormal distribution, but found {constraint_parameters_size}!")
                sys.exit(1)
        elif constraint_type == "uniform":
            if constraint_parameters_size != 2:
                print (f"ERROR: Expected 2 parameters for uniform distribution, but found {constraint_parameters_size}!")
                sys.exit(1)
        elif constraint_type == "cladeage":
            if constraint_parameters_size != 8:
                print (f"ERROR: Expected 8 parameters for cladeage, but found {constraint_parameters_size}!")
                sys.exit(1)
        else:
            if constraint_type != "monophyletic":
                print (f"ERROR: Unexpected constraint type '{constraint_type}'!")
                sys.exit(1)

        constraint_id = str(constraint_count).zfill(4)

        if constraint_type == "cladeage":
            snapp_string += f"            <distribution id=\"Clade{constraint_id}\" spec=\"beast.math.distributions.FossilPrior\" monophyletic=\"true\"  tree=\"@tree\">\n"
        else:
            snapp_string += f"            <distribution id=\"Clade{constraint_id}\" spec=\"beast.math.distributions.MRCAPrior\" "
            if constraint_placement == "stem":
                if len(constraint_clade_ary) == len(unique(table_species)):
                    print (f"ERROR: It seems that a {constraint_type} constraint should be placed on the stem of the root of the phylogeny!")
                    sys.exit(1)
                else:
                    snapp_string += "useOriginate=\"true\" "
            else:
                snapp_string += "useOriginate=\"false\" "
            snapp_string += "monophyletic=\"true\"  tree=\"@tree\">\n"
        snapp_string += f"                <taxonset id=\"Constraint{constraint_id}\" spec=\"TaxonSet\">\n"
        for s in constraint_clade_ary:
            snapp_string += f"                    <taxon idref=\"{s}\"/>\n"
        snapp_string += "                </taxonset>\n"

        if constraint_type == "normal":
            snapp_string += f"                <Normal name=\"distr\" offset=\"{constraint_parameters[0]}\">\n"
            snapp_string += f"                    <parameter estimate=\"false\" lower=\"0.0\" name=\"mean\">{constraint_parameters[1]}</parameter>\n"
            snapp_string += f"                    <parameter estimate=\"false\" lower=\"0.0\" name=\"sigma\">{constraint_parameters[2]}</parameter>\n"
            snapp_string += "                </Normal>\n"
        elif constraint_type == "lognormal":
            snapp_string += f"                <LogNormal meanInRealSpace=\"true\" name=\"distr\" offset=\"{constraint_parameters[0]}\">\n"
            snapp_string += f"                    <parameter estimate=\"false\" lower=\"0.0\" name=\"M\">{constraint_parameters[1]}</parameter>\n"
            snapp_string += f"                    <parameter estimate=\"false\" lower=\"0.0\" name=\"S\">{constraint_parameters[2]}</parameter>\n"
            snapp_string += "                </LogNormal>\n"
        elif constraint_type == "uniform":
            snapp_string += f"                <Uniform name=\"distr\" lower=\"{constraint_parameters[0]}\" upper=\"{constraint_parameters[1]}\"/>\n"
        elif constraint_type == "cladeage":
            snapp_string += "                <fossilDistr\n"
            snapp_string += f"                    id=\"{constraint_id}\"\n"
            snapp_string += f"                    minOccuranceAge=\"{constraint_parameters[0]}\"\n"
            snapp_string += f"                    maxOccuranceAge=\"{constraint_parameters[1]}\"\n"
            snapp_string += f"                    minDivRate=\"{constraint_parameters[2]}\"\n"
            snapp_string += f"                    maxDivRate=\"{constraint_parameters[3]}\"\n"
            snapp_string += f"                    minTurnoverRate=\"{constraint_parameters[4]}\"\n"
            snapp_string += f"                    maxTurnoverRate=\"{constraint_parameters[5]}\"\n"
            snapp_string += f"                    minSamplingRate=\"{constraint_parameters[6]}\"\n"
            snapp_string += f"                    maxSamplingRate=\"{constraint_parameters[7]}\"\n"
            snapp_string += "                    minSamplingGap=\"0\"\n"
            snapp_string += "                    maxSamplingGap=\"0\"\n"
            snapp_string += "                    spec=\"beast.math.distributions.FossilCalibration\"/>\n"
        snapp_string += "            </distribution>\n"

    if not no_annotation:
        snapp_string += "            <!--\n"
        snapp_string += "            The clock rate affects the model only by scaling branches before likelihood calculations.\n"
        snapp_string += "            A one-on-x prior distribution is used, which has the advantage that the shape of its\n"
        snapp_string += "            distribution is independent of the time scales used, and can therefore equally be applied\n"
        snapp_string += "            in phylogenetic analyses of very recent or old groups. However, note that the one-on-x\n"
        snapp_string += "            prior distribution is not a proper probability distribution as its total probability mass\n"
        snapp_string += "            does not sum to 1. For this reason, this prior distribution can not be used for Bayes Factor\n"
        snapp_string += "            comparison based on Path Sampling or Stepping Stone analyses. If such analyses are to be\n"
        snapp_string += "            performed, a different type of prior distribution (uniform, lognormal, gamma,...) will need\n"
        snapp_string += "            to be chosen.\n"
        snapp_string += "            -->\n"
    snapp_string += "            <prior name=\"distribution\" x=\"@clockRate\">\n"
    snapp_string += "                <OneOnX name=\"distr\"/>\n"
    snapp_string += "            </prior>\n"

    if not no_annotation:
        snapp_string += "            <!--\n"
        snapp_string += "            The scaling of branch lengths based on the clock rate does not affect the evaluation of the\n"
        snapp_string += "            likelihood of the species tree given the speciation rate lambda. Thus, lambda is measured in\n"
        snapp_string += "            the same time units as the unscaled species tree. As for the clock rate, a one-on-x prior\n"
        snapp_string += "            distribution is used, and an alternative prior distribution will need to be chosen if Bayes\n"
        snapp_string += "            Factor comparisons are to be performed.\n"
        snapp_string += "            -->\n"
    snapp_string += "            <prior name=\"distribution\" x=\"@lambda\">\n"
    snapp_string += "                <OneOnX name=\"distr\"/>\n"
    snapp_string += "            </prior>\n"

    if not no_annotation:
        snapp_string += "            <!--\n"
        snapp_string += "            The below distribution defines the prior probability for the population mutation rate Theta.\n"
        snapp_string += "            In standard SNAPP analyses, a gamma distribution is commonly used to define this probability,\n"
        snapp_string += "            with a mean according to expectations based on the assumed mean effective population size\n"
        snapp_string += "            (across all branches) and the assumed mutation rate per site and generation. However, with\n"
        snapp_string += "            SNP matrices, the mutation rate of selected SNPs usually differs strongly from the genome-wide\n"
        snapp_string += "            mutation rate, and the degree of this difference depends on the way in which SNPs were selected\n"
        snapp_string += "            for the analysis. The SNP matrices used for this analysis are assumed to all be bi-allelic,\n"
        snapp_string += "            excluding invariable sites, and are thus subject to ascertainment bias that will affect the\n"
        snapp_string += "            mutation rate estimate. If only a single species would be considered, this ascertainment bias\n"
        snapp_string += "            could be accounted for (Kuhner et al. 2000; Genetics 156: 439–447). However, with multiple\n"
        snapp_string += "            species, the proportion of SNPs that are variable among individuals of the same species,\n"
        snapp_string += "            compared to the overall number of SNPs variable among all species, depends on relationships\n"
        snapp_string += "            between species and the age of the phylogeny, parameters that are to be inferred in the analysis.\n"
        snapp_string += "            Thus, SNP ascertainment bias can not be accounted for before the analysis, and the prior\n"
        snapp_string += "            expectation of Theta is extremely vague. Therefore, a uniform prior probability distribution\n"
        snapp_string += "            is here used for this parameter, instead of the commonly used gamma distribution. By default,\n"
        snapp_string += "            SNAPP uses a lower boundary of 0 and an upper boundary of 10000 when a uniform prior probability\n"
        snapp_string += "            distribution is chosen for Theta, and these lower and upper boundaries can not be changed without\n"
        snapp_string += "            editing the SNAPP source code. Regardless of the wide boundaries, the uniform distribution works\n"
        snapp_string += "            well in practice, at least when the Theta parameter is constrained to be identical on all branches\n"
        snapp_string += "            (see below).\n"
        snapp_string += "            -->\n"
    snapp_string += "            <distribution spec=\"snap.likelihood.SnAPPrior\" coalescenceRate=\"@coalescenceRate\" lambda=\"@lambda\" rateprior=\"uniform\" tree=\"@tree\">\n"

    if not no_annotation:
        snapp_string += "                <!--\n"
        snapp_string += "                SNAPP requires input for parameters alpha and beta regardless of the chosen type of prior,\n"
        snapp_string += "                however, the values of these two parameters are ignored when a uniform prior is selected.\n"
        snapp_string += "                Thus, they are both set arbitrarily to 1.0.\n"
        snapp_string += "                -->\n"

    snapp_string += "                <parameter estimate=\"false\" lower=\"0.0\" name=\"alpha\">1.0</parameter>\n"
    snapp_string += "                <parameter estimate=\"false\" lower=\"0.0\" name=\"beta\">1.0</parameter>\n"
    snapp_string += "            </distribution>\n"
    snapp_string += "        </distribution>\n"
    snapp_string += "        <distribution id=\"likelihood\" spec=\"util.CompoundDistribution\">\n"
    snapp_string += "            <distribution spec=\"snap.likelihood.SnAPTreeLikelihood\" data=\"@data\" non-polymorphic=\"false\" pattern=\"coalescenceRate\" tree=\"@tree\">\n"
    snapp_string += "                <siteModel spec=\"SiteModel\">\n"
    snapp_string += "                    <substModel spec=\"snap.likelihood.SnapSubstitutionModel\" coalescenceRate=\"@coalescenceRate\">\n"

    if not no_annotation:
        snapp_string += "                        <!--\n"
        snapp_string += "                        The forward and backward mutation rates are fixed so that the number of expected mutations\n"
        snapp_string += "                        per time unit (after scaling branch lengths with the clock rate) is 1.0. This is done to\n"
        snapp_string += "                        avoid non-identifability of rates, given that the clock rate is estimated. Both parameters\n"
        snapp_string += "                        are fixed at the same values, since it is assumed that alleles were translated to binary\n"
        snapp_string += "                        code by random assignment of '0' and '2' to homozygous alleles, at each site individually.\n"
        snapp_string += "                        Thus, the probabilities for '0' and '2' are identical and the resulting frequencies of '0'\n"
        snapp_string += "                        and '2' in the data matrix should be very similar.\n"
        snapp_string += "                        -->\n"
    snapp_string += "                        <parameter estimate=\"false\" lower=\"0.0\" name=\"mutationRateU\">1.0</parameter>\n"
    snapp_string += "                        <parameter estimate=\"false\" lower=\"0.0\" name=\"mutationRateV\">1.0</parameter>\n"
    snapp_string += "                    </substModel>\n"
    snapp_string += "                </siteModel>\n"

    if not no_annotation:
        snapp_string += "                <!--\n"
        snapp_string += "                A strict clock rate is used, assuming that only closely related species are used in SNAPP\n"
        snapp_string += "                analyses and that branch rate variation among closely related species is negligible.\n"
        snapp_string += "                The use of a relaxed clock is not supported in SNAPP.\n"
        snapp_string += "                -->\n"
    snapp_string += "                <branchRateModel spec=\"beast.evolution.branchratemodel.StrictClockModel\" clock.rate=\"@clockRate\"/>\n"
    snapp_string += "            </distribution>\n"
    snapp_string += "        </distribution>\n"
    snapp_string += "    </distribution>\n"
    snapp_string += "\n"
    snapp_string += "    <!-- Operators -->\n"

    if weight > 0:
        if not no_annotation:
            snapp_string += "    <!--\n"
            snapp_string += "    The treeNodeSwapper operator is the only operator on the tree topology. Thus if the tree topology\n"
            snapp_string += "    should be fixed, simply remove this operator by deleting the second-next line.\n"
            snapp_string += "    -->\n"
        snapp_string += f"    <operator id=\"treeNodeSwapper\" spec=\"snap.operators.NodeSwapper\" tree=\"@tree\" weight=\"{weight}\"/>\n"
    if not no_annotation:
        snapp_string += "    <!--\n"
        snapp_string += "    The treeNodeBudger and treeScaler operators modify node heights of the tree.\n"
        snapp_string += "    -->\n"
    snapp_string += "    <operator id=\"treeNodeBudger\" spec=\"snap.operators.NodeBudger\" size=\"0.75\" tree=\"@tree\" weight=\"1.0\"/>\n"
    snapp_string += "    <operator id=\"treeScaler\" spec=\"snap.operators.ScaleOperator\" scaleFactor=\"0.75\" tree=\"@tree\" weight=\"1.0\"/>\n"
    if not no_annotation:
        snapp_string += "    <!--\n"
        snapp_string += "    To constrain the Theta parameter so that all branches always share the same value (and thus\n"
        snapp_string += "    the same population size estimates), a single operator is used to modify Theta values by scaling\n"
        snapp_string += "    the Thetas of all branches up or down by the same factor. Instead, SNAPP's default Theta operator\n"
        snapp_string += "    types 'GammaMover' and 'RateMixer' are not not used.\n"
        snapp_string += "    -->\n"
    snapp_string += "    <operator id=\"thetaScaler\" spec=\"snap.operators.ScaleOperator\" parameter=\"@coalescenceRate\" scaleFactor=\"0.75\" scaleAll=\"true\" weight=\"1.0\"/>\n"
    if not no_annotation:
        snapp_string += "    <!--\n"
        snapp_string += "    The lamdaScaler and clockScaler operators modify the speciation and clock rates.\n"
        snapp_string += "    -->\n"
    snapp_string += "    <operator id=\"lamdaScaler\" spec=\"snap.operators.ScaleOperator\" parameter=\"@lambda\" scaleFactor=\"0.75\" weight=\"1.0\"/>\n"
    snapp_string += "    <operator id=\"clockScaler\" spec=\"snap.operators.ScaleOperator\" parameter=\"@clockRate\" scaleFactor=\"0.75\" weight=\"1.0\"/>\n"
    snapp_string += "\n"
    snapp_string += "    <!-- Loggers -->\n"
    snapp_string += f"    <logger fileName=\"{snapp_log_file_name}\" logEvery=\"{int(store_frequency)}\">\n"
    snapp_string += "        <log idref=\"posterior\"/>\n"
    snapp_string += "        <log idref=\"likelihood\"/>\n"
    snapp_string += "        <log idref=\"prior\"/>\n"
    snapp_string += "        <log idref=\"lambda\"/>\n"
    snapp_string += "        <log id=\"treeHeightLogger\" spec=\"beast.evolution.tree.TreeHeightLogger\" tree=\"@tree\"/>\n"
    snapp_string += "        <log idref=\"clockRate\"/>\n"
    snapp_string += "    </logger>\n"
    snapp_string += f"    <logger logEvery=\"{screen_frequency}\">\n"
    snapp_string += "        <log idref=\"posterior\"/>\n"
    snapp_string += "        <log spec=\"util.ESS\" arg=\"@posterior\"/>\n"
    snapp_string += "        <log idref=\"likelihood\"/>\n"
    snapp_string += "        <log idref=\"prior\"/>\n"
    snapp_string += "        <log idref=\"treeHeightLogger\"/>\n"
    snapp_string += "        <log idref=\"clockRate\"/>\n"
    snapp_string += "    </logger>\n"
    snapp_string += f"    <logger fileName=\"{snapp_trees_file_name}\" logEvery=\"{int(store_frequency)}\" mode=\"tree\">\n"
    snapp_string += "        <log spec=\"beast.evolution.tree.TreeWithMetaDataLogger\" tree=\"@tree\">\n"
    snapp_string += "            <metadata id=\"theta\" spec=\"snap.RateToTheta\" coalescenceRate=\"@coalescenceRate\"/>\n"
    snapp_string += "        </log>\n"
    snapp_string += "    </logger>\n"
    snapp_string += "\n"
    snapp_string += "</run>\n"
    snapp_string += "\n"
    snapp_string += "</beast>\n"


    # Write the SNAPP input file.
    with open(xml, "w") as snapp_file:
        snapp_file.write(snapp_string)
    print (f"Wrote SNAPP input in XML format to file {xml}.\n")



if __name__ == '__main__':
    # Get the command line options and define defaults if script is used in CLI
    opt_parser = argparse.ArgumentParser(description="This script prepares XML format input files for the software SNAPP (http://beast2.org/snapp/), originally wrote by Michael Matschiner (michaelmatschiner@mac.com) 2018, and ported to Python3 by Carlos Alonso Maya-Lastra 2020", formatter_class=argparse.ArgumentDefaultsHelpFormatter)

    opt_parser.add_argument("-p", "--phylip", metavar="FILENAME", help="File with SNP data in phylip format.", default=None, dest="phylip")
    opt_parser.add_argument("-v","--vcf", metavar="FILENAME", help="File with SNP data in vcf format. [NOT SUPPORTED IN THIS PORT]", default=None, dest="vcf")
    opt_parser.add_argument("-t","--table", metavar="FILENAME", help="File with table linking species and specimens.", default="example.spc.txt", dest="table")
    opt_parser.add_argument("-c","--constraints", metavar="FILENAME", help="File with age constraint information.", default="example.con.txt", dest="constraints")
    opt_parser.add_argument("-s","--starting-tree", metavar="FILENAME", help="File with starting tree in Nexus or Newick format.", default=None, dest="tree")
    opt_parser.add_argument("-l","--length", metavar="LENGTH", type=int, help="Number of MCMC generations.", default=500000, dest="length")
    opt_parser.add_argument("-w","--weight", metavar="WEIGHT", type=float, help="Relative weight of topology operator.", default=1.0, dest="weight")
    opt_parser.add_argument("-m","--max-snps", metavar="NUMBER", type=int, help="Maximum number of SNPs to be used.", default=None, dest="max_snps")
    opt_parser.add_argument("-q","--min-dist", metavar="NUMBER", type=int, help="Minimum distance between SNPs to be used.", default=None, dest="min_dist")
    opt_parser.add_argument("-r","--transversions", action='store_true', help="Use transversions only.", default=False, dest="transversions")
    opt_parser.add_argument("-i","--transitions", action='store_true', help="Use transitions only.", default=False, dest="transitions")
    opt_parser.add_argument("-x","--xml", metavar="FILENAME", help="Output file in XML format.", default="snapp.xml", dest="xml")
    opt_parser.add_argument("-o","--out", metavar="PREFIX", help="Prefix for SNAPP's.log and.trees output files.", default="snapp", dest="out")
    opt_parser.add_argument("-n","--no-annotation", action='store_true', help="Do not add explanatory annotation to XML file.", default=False, dest="no_annotation")

    options = opt_parser.parse_args()

    #Print help if no argument is passed
    if len(sys.argv)==1:
        opt_parser.print_help(sys.stderr)
        sys.exit(1)

    run(phylip=options.phylip, vcf=options.vcf,table=options.table,constraints=options.constraints,tree=options.tree,length=options.length,weight=options.weight,max_snps=options.max_snps,min_dist=options.min_dist,transversions=options.transversions,transitions=options.transitions,xml=options.xml,out=options.out,no_annotation=options.no_annotation)
