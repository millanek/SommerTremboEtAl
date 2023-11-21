# Authors: Michael Matschiner and Milan Malinsky 
# Last edit: 21st Nov 2023

# Import libraries.
import msprime
import sys
import newick
import random
from random import randint
import datetime
import re
import os
import numpy
import argparse

def get_generations_per_branch_length_unit(
        branch_length_units=None,
        generation_time=None):
    """
    Method to calculate the number of generations per branch length
    unit, given the branch length unit and a generation time.
    """
    if branch_length_units == "gen":
        generations_per_branch_length_unit = 1
    elif branch_length_units == "myr":
        generations_per_branch_length_unit = 10**6/generation_time
    else:
        generations_per_branch_length_unit = 1/generation_time
    return generations_per_branch_length_unit


def parse_species_tree(
        species_tree=None,
        branch_length_units="gen",
        Ne=None,
        generation_time=None,
        migration_matrix=None,
        geneFlowPeriod=None):
    """
    Method to parse species trees in Newick
    (https://en.wikipedia.org/wiki/Newick_format) format.

    Trees are assumed to be rooted and ultrametric and branch lengths
    must be included and correspond to time, either in units of millions
    of years ("myr"), years ("yr"), or generations ("gen"; default).
    Leafs must be named. An example for an accepted tree string in
    Newick format is:
    (((human:5.6,chimpanzee:5.6):3.0,gorilla:8.6):9.4,orangutan:18.0)
    The tree string can end with a semi-colon, but this is not required.

    - An estimate of the effective population size Ne should be
        specified.
    - If and only if the branch lengths are not in units of
        generations, the generation time should be specified.
    """

    # Make sure a species tree is specified.
    if type(species_tree) is not str:
        raise ValueError("A species tree must be specified.")

    # Make sure that branch length units are either "myr", "yr", or "gen".
    allowed_branch_lenth_units = ["myr", "yr", "gen"]
    if branch_length_units not in allowed_branch_lenth_units:
        err = 'The specified units for branch lengths ('
        err += '"{}") are not accepted. '.format(branch_length_units)
        err += 'Accepted units are "myr" (millions of years), "yr" (years), '
        err += 'and "gen" (generations).'
        raise ValueError(err)

    # Make sure that the population size is either None or positive.
    if Ne is not None:
        try:
            Ne = float(Ne)
        except ValueError:
            raise ValueError("Population size Ne must be numeric.")
        if Ne <= 0:
            raise ValueError("Population size Ne must be > 0.")

    # Make sure that the generation time is either None or positive.
    if generation_time is not None:
        try:
            generation_time = float(generation_time)
        except ValueError:
            raise ValueError("Generation time must be numeric.")
        if generation_time <= 0:
            raise ValueError("Generation time must be > 0.")

    # Make sure that the generation time is specified if and only if
    # branch lengths are not in units of generations.
    if branch_length_units == "gen":
        if generation_time is not None:
            err = 'With branch lengths in units of generations ("gen"), '
            err += 'a generation time should not be specified additionally.'
            raise ValueError(err)
    else:
        if generation_time is None:
            err = 'With branch lengths in units of '
            err += '"{}", a generation time must be '.format(branch_length_units)
            err += 'specified additionally.'
            raise ValueError(err)

    # Make sure that a population size is specified.
    if Ne is None:
        raise ValueError("Ne should be specified.")

    # Get the number of generations per branch length unit.
    generations_per_branch_length_unit = get_generations_per_branch_length_unit(
        branch_length_units, generation_time
        )
    # Get the gene flow period limit parameter in the time units of the tree
    geneFlowPeriodInTreeUnits = geneFlowPeriod / generations_per_branch_length_unit

    # Read the input file.
    species_tree_lines = species_tree.splitlines(False)
    if len(species_tree_lines) > 1:
        raise ValueError("The species tree has multiple lines.")
    tree_patterns = re.search('\\(.+\\)', species_tree_lines[0])
    if tree_patterns is None:
        raise ValueError("No species tree string found.")
    species_tree_string = tree_patterns.group(0)

    # Parse the newick tree string.
    root = newick.loads(species_tree_string)[0]

    # Set node depths (distances from root).
    max_depth = 0
    for node in root.walk():
        depth = 0
        moving_node = node
        while moving_node.ancestor:
            depth += moving_node.length
            moving_node = moving_node.ancestor
        node.depth = depth
        if depth > max_depth:
            max_depth = depth

    # Set node heights (distances from present).
    for node in root.walk():
        node.height = max_depth - node.depth

    # Get a list of species IDs.
    species_ids = sorted(root.get_leaf_names())
    
    # Create a new matrix which will hold the migration rates for time 0.0 - to be supplied as a param to msprime.simulate
    timeZero_migration_matrix = [[0.0 for col in range(len(species_ids))] for row in range(len(species_ids))]

    # Determine at which time which populations should merge.
    sources = []
    destinations = []
    divergence_times = []
    allDescendantsOfDestinations = []
    allDescendantsOfSources = []
    for node in root.walk():
        if node.is_leaf is False:
            name_indices = []
            leaf_names_of_nodeDescendants = []
            for descendants in node.descendants:
                leaf_names = (descendants.get_leaf_names())
                name_indices.append(species_ids.index(sorted(leaf_names)[0]))
                leaf_names_of_nodeDescendants.append(sorted(leaf_names))
            
            s = sorted(zip(name_indices, leaf_names_of_nodeDescendants))
            name_indices, leaf_names_of_nodeDescendants = map(list, zip(*s))
            new_destination = name_indices[0]
            name_indices.remove(name_indices[0])
            newDescendantsOfDestination = leaf_names_of_nodeDescendants[0]
            leaf_names_of_nodeDescendants.remove(leaf_names_of_nodeDescendants[0])
            for x in range(len(name_indices)):
                sources.append(name_indices[x])
                allDescendantsOfSources.append(leaf_names_of_nodeDescendants[x])
                destinations.append(new_destination)
                allDescendantsOfDestinations.append(newDescendantsOfDestination)
                divergence_times.append(node.height)

    # Sort the lists source_sets, destinations, and divergence_times
    # according to divergence_time.
    s = sorted(zip(divergence_times, sources, destinations, allDescendantsOfSources, allDescendantsOfDestinations))
    divergence_times, sources, destinations, allDescendantsOfSources, allDescendantsOfDestinations = map(list, zip(*s))
    
   
   #DEBUG STUFF:
#    print("Batmin: " + str(species_ids.index("Batmin")))
#    print("Batleo: " + str(species_ids.index("Batleo")))
#    print("Batgra: " + str(species_ids.index("Batgra")))
#    print("Batfer: " + str(species_ids.index("Batfer")))
#    print("Bathor: " + str(species_ids.index("Bathor")))
#    print("Batfas: " + str(species_ids.index("Batfas")))
#    print("Batvit: " + str(species_ids.index("Batvit")))
    
    
    for x in range(len(divergence_times)):
        # Start of gene flow ('Start' in the sense of looking backwards in coalescent terms)
        geneFlowStart =  divergence_times[x] - geneFlowPeriodInTreeUnits
        
        if geneFlowStart < 0.0:
            geneFlowStart = 0.0
        else:
            descendantsToRemove = [] # These should not have any gene-flow, because they do not exist within the geneFlowPeriod for this node
            for s in allDescendantsOfSources[x]:
                speciesID = species_ids.index(s)
                for y in range(len(divergence_times)):
                    if speciesID == sources[y]:
                        if divergence_times[y] <= geneFlowStart:
                            descendantsToRemove.append(s)
            
            for r in descendantsToRemove:
                allDescendantsOfSources[x].remove(r)
            
            descendantsToRemove = []
            
            for s in allDescendantsOfDestinations[x]:
                speciesID = species_ids.index(s)
                for y in range(len(divergence_times)):
                    if speciesID == sources[y]:
                        if divergence_times[y] <= geneFlowStart:
                            descendantsToRemove.append(s)
                            
            for r in descendantsToRemove:
                allDescendantsOfDestinations[x].remove(r)
    

    # Define the species/population tree for msprime.
    population_configurations = []
    for _ in range(len(root.get_leaves())):
        population_configurations.append(
            msprime.PopulationConfiguration(
                initial_size=Ne))
    demographic_events = []
    demographic_event_times = []
    for x in range(len(divergence_times)):
        demographic_events.append(
            msprime.MassMigration(
                time=divergence_times[x]*generations_per_branch_length_unit,
                source=sources[x],
                destination=destinations[x]))
        demographic_event_times.append(divergence_times[x])
        
        
        # Start of gene flow ('Start' in the sense of looking backwards in coalescent terms)
        geneFlowStart =  divergence_times[x] - geneFlowPeriodInTreeUnits
    #DEBUG STUFF:
#        print("geneFlowPeriodInTreeUnits: ")
#        print(geneFlowPeriodInTreeUnits)
#        print("divergence_times[x]: ")
#        print(divergence_times[x])
#        print("geneFlowStart: ")
#        print(geneFlowStart)
        
        if geneFlowStart < 0.0:
            geneFlowStart = 0.0
            for d in allDescendantsOfDestinations[x]:
                for s in allDescendantsOfSources[x]:
                    timeZero_migration_matrix[species_ids.index(s)][species_ids.index(d)] = migration_matrix[species_ids.index(s)][species_ids.index(d)]
                    timeZero_migration_matrix[species_ids.index(d)][species_ids.index(s)] = migration_matrix[species_ids.index(d)][species_ids.index(s)]
            
 #DEBUG STUFF:
#        print("geneFlowStart: ")
#        print(geneFlowStart)
#
#        print("")
        else:
        # Now loop over all descendants of this node which should have gene-flow among them and set that
            for d in allDescendantsOfDestinations[x]:
                for s in allDescendantsOfSources[x]:
        #DEBUG STUFF:
    #               if divergence_times[x] > 1.5:
    #                    print(d)
    #                    print(s)
    #                    print(species_ids[8])
    #                    print(species_ids[67])
    #                    print(species_ids[68])
    #                    sys.exit(1)
                    if migration_matrix[species_ids.index(s)][species_ids.index(d)] > 0.0:
                        demographic_events.append(
                        msprime.MigrationRateChange(
                            time=geneFlowStart*generations_per_branch_length_unit,
                            rate=migration_matrix[species_ids.index(s)][species_ids.index(d)],
                            matrix_index=(species_ids.index(s), species_ids.index(d))))
                        demographic_event_times.append(geneFlowStart)
                        
                    if migration_matrix[species_ids.index(d)][species_ids.index(s)] > 0.0:
                        demographic_events.append(
                        msprime.MigrationRateChange(
                            time=geneFlowStart*generations_per_branch_length_unit,
                            rate=migration_matrix[species_ids.index(d)][species_ids.index(s)],
                            matrix_index=(species_ids.index(d), species_ids.index(s))))
                            
                        demographic_event_times.append(geneFlowStart)
                    
                    # Of course better if the rates in internal branches were averages, not just randomply selected 'surviving' ones
        
        for y in range(len(root.get_leaves())):
            if y != sources[x]:
                demographic_events.append(
                    msprime.MigrationRateChange(
                        time=divergence_times[x]*generations_per_branch_length_unit,
                        rate=0,
                        matrix_index=(sources[x], y)))
                demographic_event_times.append(divergence_times[x])
                
                demographic_events.append(
                    msprime.MigrationRateChange(
                        time=divergence_times[x]*generations_per_branch_length_unit,
                        rate=0,
                        matrix_index=(y, sources[x])))
                demographic_event_times.append(divergence_times[x])
                
    
    # Sort the list demographic_events according to demographic_event_times (the "zip" approach does't work on the type in demographic_events)
    if len(demographic_event_times) == len(demographic_events):
        desiredOrder = numpy.argsort(demographic_event_times)
        sorted_demographic_events = []
        for x in range(len(desiredOrder)):
            sorted_demographic_events.append(demographic_events[desiredOrder[x]])
    else:
        print("Something went horribly wrong: len(demographic_event_times) != len(demographic_events)")
        sys.exit(1)

    # Return a tuple of population_configurations, demographic_events, and the initial migration_matrix at time 0.0
    return population_configurations, sorted_demographic_events, timeZero_migration_matrix

# Get the command line arguments/parameters
argparser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter,
                                    description="Build the msprime model and run simulations given a tree and a migration matrix based on f4-ratio statistics produced by Dsuite. ", add_help=True)
argparser.add_argument("tree", type = str, help="Path to a .newick tree file")
argparser.add_argument("f4ratios", type=str, help="Path to file containing f4-ratio statitics matrix derived from Dsuite results")
argparser.add_argument("vcf", type=str, help="Path to the output file which will contain the simulation results")
argparser.add_argument("-g", "--gen_time", type=float, help="Average generation time in years.", default=3.0)
argparser.add_argument("-u", "--mut_rate", type=float, help="Mutation rate per bp per generation",default=3.5e-9)
argparser.add_argument("-r", "--rho", type=float, help="Recombination rate per bp per generation (in the future could use a recombination map); set to 0 to simulate a non-recombining region", default=2.2e-8)
argparser.add_argument("-N", "--Ne", type=int, help="Effective population size (for now this is the same throughout the tree.", default=20000)
argparser.add_argument("-l","--chr_length", type=int, help="Length of the 'chromosome' to simulate", default=100)
argparser.add_argument("-p", "--gene_flow_period", help="How long gene-flow persist after species/populations split (in generations).", type=int, default=1000000)
argparser.add_argument("-s", "--scaling_factor", help="A scaling factor between the f4-ratios and the msprime migration units", type=float, default=1.0e-5)
argparser.add_argument("-n", "--num_indiv", help="The number of (diploid) individuals to simulate from each population/species", type=int, default=1)
argparser.add_argument("-d", "--debugger", help="Save the msprime.DemographyDebugger output to FILENAME", type=str, default=None,metavar="FILENAME")
args = argparser.parse_args()

# Read the species tree string.
with open(args.tree, "r") as f:
    species_tree_string = f.read()

# Parse the species tree.
root = newick.loads(species_tree_string)[0]
species_ids = sorted(root.get_leaf_names())


# Check if the specified migration rate string is a number or a file name.
print('Preparing the migration matrix...', end='', flush=True)
if os.path.isfile(args.f4ratios):
    migrations = []
    # Read the migration rate matrix.
    with open(args.f4ratios) as f:
        dsuite_file_lines = f.readlines()
        dsuite_file_species_ids = dsuite_file_lines[0].split()
        migration_matrix = []
        for species1 in species_ids:
            row = []
            if species1 in dsuite_file_species_ids:
                dsuite_species1_index = dsuite_file_species_ids.index(species1)
            else:
                dsuite_species1_index = None
          
            for species2 in species_ids:
                if species2 in dsuite_file_species_ids:
                    dsuite_species2_index = dsuite_file_species_ids.index(species2)
                else:
                    dsuite_species2_index = None
   
                if species1 == species2:
                    row.append(0)
                else:
                    if dsuite_species1_index is None or dsuite_species2_index is None:
                        row.append(0)
                    else:
                        row_list = dsuite_file_lines[dsuite_species1_index+1].split()
                        migration = float(row_list[dsuite_species2_index+1]) * args.scaling_factor
                        row.append(migration)
                        migrations.append(migration)

            migration_matrix.append(row)
    print(" done. Mean migration rate is " + str(sum(migrations)/len(migrations)) + ".")

    
#print(migration_matrix)
    
# Parse the species tree with msprime and generate population configurations and demographic evens.
parsed_tuple = parse_species_tree(
    species_tree=species_tree_string,
    branch_length_units="myr",
    Ne=args.Ne,
    generation_time=args.gen_time,
    migration_matrix=migration_matrix,
    geneFlowPeriod=args.gene_flow_period
    )
population_configurations = parsed_tuple[0]
demographic_events = parsed_tuple[1]
timeZero_migration_matrix = parsed_tuple[2]
for n in population_configurations:
    n.sample_size = 2 * args.num_indiv

if not args.debugger is None:
    dd = msprime.DemographyDebugger(
            population_configurations=parsed_tuple[0],
            demographic_events=parsed_tuple[1],
            migration_matrix=timeZero_migration_matrix)
    dd_outfile = open(args.debugger, 'w')
    dd.print_history(dd_outfile)
#sys.exit(0)

# Write the vcf file.
print('Simulating with msprime...', end='', flush=True)
new_tree_sequence_obj = msprime.simulate(
        population_configurations=population_configurations,
        migration_matrix=timeZero_migration_matrix,
        demographic_events=demographic_events,
        mutation_rate=args.mut_rate,
        length=args.chr_length,
        recombination_rate=args.rho,
        random_seed=random.randint(1, 10000000)
        )
print(" done.")

# Prepare the vcf header.
print('Preparing the vcf...', end='', flush=True)
vcf_string = '##fileformat=VCFv4.2\n'
now = datetime.datetime.now()
vcf_string += '##fileDate={}\n'.format(now.strftime("%Y%m%d"))
vcf_string += '##source=c-genie\n'
vcf_string += '##contig=<ID=1,length={}>\n'.format(args.chr_length)
vcf_string += '##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">\n'
vcf_string += '#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT'
for sp in range(0,len(species_ids)):
    for i in range(args.num_indiv):
        vcf_string += '\t{}_{}'.format(species_ids[sp],i)
vcf_string += '\n'

# Prepare vcf file.
DNA_alphabet = ['A','C','G','T']
n_variants = 0
in_gene_position = 0
for variant in new_tree_sequence_obj.variants():
    vp = int(variant.site.position)
    if vp > in_gene_position:  # exclude multiple mutations at the same site (i.e. at the moment we generate only biallelic sites) TO DO: make multiallelics in these cases
        in_gene_position = vp
        n_variants = n_variants + 1
        ra = randint(0, 3)
        ancestralBase = DNA_alphabet[ra]
        derivedPossibilities = DNA_alphabet[:ra] +  DNA_alphabet[ra+1:]
        rd = randint(0, 2)
        derivedBase = derivedPossibilities[rd]
        vcf_string += '1\t{}\t.\t{}\t{}\t.\t.\t.\tGT'.format(vp,ancestralBase,derivedBase)
        for sp in range(0,len(species_ids)):
            for i in range(args.num_indiv):
                if i == 1:
                    vcf_string += '\t{}|{}'.format(variant.genotypes[2*args.num_indiv*sp],variant.genotypes[(2*args.num_indiv*sp)+1])
                else:
                    vcf_string += '\t{}|{}'.format(variant.genotypes[(2*args.num_indiv*sp)+(i*2)],variant.genotypes[(2*args.num_indiv*sp)+(i*2)+1])
        vcf_string += '\n'

# Write the vcf output file.
vcf_outfile = open(args.vcf, 'w')
vcf_outfile.write(vcf_string)
print(" done.")
