import algebra as va
import networkx as nx
import warnings
from .data import api_get
import re
import copy

all_functions = ['unknown function', 'uncertain function', 'normal function', 'decreased function', 'no function']

def sort_function(f):
    """Sort function annotation based on severity."""
    # QUESTION: what is the difference between uncertain and unknown?
    if f not in all_functions:
        raise Exception("Unknown function: " + f)
    functions = all_functions[3:]
    if f not in functions: 
        return 0 
    return functions.index(f) + 1

def sort_types(v):
    """Sort variant types in order of specificity (somewhat arbitrary)."""
    if '*' in v:
        if '.' in v:
            return 2 # Suballele
        else:
            return 1 # Core allele
    elif v[:2] in ('HG', 'NA'):
        return 4 # Sample
    else:
        if ':' in v:
            return 3 # Variant
        return 5 # Personal variant
    
def find_contained_variants(start, cont_graph, eq_graph, matches, visited, find, stop=None):
    """Recursively find the equivalent and contained variants from a given start node."""
    if start in visited: return # Already visited
    visited.add(start) # Needed to avoid equivalence loop and to avoid doubles
    if sort_types(start) in find: # Need to find this type
        matches.append(start)
    if stop is not None and sort_types(start) == stop: # Don't look further here
        return 
    # Find equivalents
    if start in eq_graph.nodes():
        for match in eq_graph[start]: 
            find_contained_variants(match, cont_graph, eq_graph, matches, visited, find, stop) # Add equivalents and maybe traverse
    # Find contained
    if start in cont_graph.nodes():
        for match, _ in cont_graph.in_edges(start):
            find_contained_variants(match, cont_graph, eq_graph, matches, visited, find, stop) # Add contained  and maybe traverse
        
def find_best_match(sample, matches, functions, phased):
    """Find the best matches from a list of matches.
    
    Also return a measure of certainty.
    """
    # TODO change to list of sets (since there is no order within a position)
    sorted_matches = []
    # Equivalent alleles are the best matches
    if len(matches["equivalent"]) > 0:
        if len(matches["equivalent"]) > 1: # The same as some allele
            raise Exception(f"{sample}: This should not happen, multiple equivalent matches found: {matches['equivalent']}.")
        sorted_matches.append(matches["equivalent"])
    # Next step is the contained alleles
    if len(matches["contained"]) > 0:
        # In the case of multiple contained corealleles, a choice is made based on functional annotation
        annotated = {f: [m for m in matches["contained"] if functions[m] == f] for f in all_functions}
        # If the annotation is no function, this allele is most disruptive and is chosen
        if len(annotated["no function"]) > 0:	
            sorted_matches.append(annotated["no function"])
        # If uncertain or unknown functions are present, the ordering after no function cannot be determined
        # An unknown function could be no function and thus no allele should be preferred.
        if len(annotated["unknown function"]) > 0 or len(annotated["uncertain function"]) > 0:
            # Put all in the same position of ordering
            sorted_matches.append(
                annotated["decreased function"] +
                annotated["normal function"] +
                annotated["unknown function"] +
                annotated["uncertain function"]
            )
        else: # Ordering can be determined between normal and decreased function
            # More disruptive alleles are chosen first
            if len(annotated["decreased function"]) > 0:
                sorted_matches.append(annotated["decreased function"])
            if len(annotated["normal function"]) > 0:
                sorted_matches.append(annotated["normal function"])
    # Check certainty and consistency of matches
    for match in sorted_matches:
        n_cores = 0
        for allele in match:
            if sort_types(allele) == 1: # Coreallele
                n_cores += 1
            # check certainty
            if len(matches['variants'][allele]['uncertain']) > 0:
                # TODO how to express this to the user?
                # warnings.warn(f"{sample}={allele}: Uncertain variants found: {matches['variants'][allele]['uncertain']}.")
                pass
        if n_cores > 1 and phased or n_cores > 2 and not phased: # Ambiguous phased/unphased
            raise Exception(f"{sample}: Multiple corealleles found with the same priority: {match}.")
        elif n_cores == 0: # No found yet
            continue
        else: # Unambiguous
            break # No need to check further
    else: # No cores found
        # Catch-all is the default allele
        # *1 will never be in contained so it is not prioritized over others
        # TODO also find uncertain variants here
        sorted_matches.append(["CYP2D6*1"]) # QUESTION is this valid?
    return sorted_matches

def star_allele_calling(sample, eq_graph, cont_graph, functions, supremals, phased):
    """Determine star allele calling for a sample based on va relations.
    
    Find based on pruned graph containing relations between samples and alleles and between themselves.
    """
    # QUESTION: is it needed to look at suballeles for calling?
    # QUESTION: is it needed to look at individual variants for calling?
    # All information of a calling
    matches = {
        "equivalent": [], # Equivalent alleles
        "contained": [], # Contained alleles
        "variants": {
            # extra variants for each allele
        }
    }
    # Find equivalent alleles
    contained = []
    if sample in eq_graph.nodes(): matches["equivalent"] = [match for match in eq_graph[sample] if sort_types(match) != 5]
    # Find contained alleles
    find_contained_variants(sample, cont_graph, eq_graph, contained, set(), (1,2), 1)
    # Check if any contain another, keep most specific
    # TODO add to initial getting?
    for node1 in set(contained):
        for node2 in contained:
            if node1 == node2: # Skip self 
                continue
            if sort_types(node1) == 1 and sort_types(node2) == 2: # Skip containment of core in suballeles as this is expected
                continue
            if node1 in nx.ancestors(cont_graph, node2): # Node1 is contained in node2
                break # Skip this allele
        else: # No issues 
            if node1 in matches["equivalent"]: # Skip equivalent alleles
                continue
            matches["contained"].append(node1) # Keep match
    # Check certainty of calling for each matched allele
    vs = [] # All variants indirectly or directly contained in the sample 
    find_contained_variants(sample, cont_graph, eq_graph, vs, set(), (3, 5)) # TODO get from PharmVar?
    for match in matches["contained"] + matches["equivalent"]: # Equivalent added but these will be empty
        matches["variants"][match] = {"noise": [], "uncertain": []}
        # Find 'extra' variants: in the sample but not in the called allele
        ws = [] # All variants indirectly or directly contained in the sample
        find_contained_variants(match, cont_graph, eq_graph, ws, set(), (3, 5)) # TODO get from PharmVar?
        extra_variants = []
        for v in vs: # In sample
            if v in ws: # But not in this allele
                continue
            extra_variants.append(v) 
        # Check if the extra variants can influence the calling
        for variant in extra_variants:
            # Check what the impact of the variant is
            noise, function = is_noise(variant, functions, supremals[variant])
            if noise: # No evidence of relevance
                matches["variants"][match]["noise"].append(variant) # TODO annotate noise
                continue # TODO need to check interference here?
            # Evidence of possible protein impact
            # Check if the impact could interfere with the calling in some way
            if "splice defect" in function: # Splice defect, outcome uncertain
                pass
            elif "fs" == function[-2:]: # Frameshift
                # TODO check if effect downstream
                pass
            elif "X" == function[-1]: # Stop codon
                # TODO check if effect downstream
                pass
            else: # Missense
                # TODO check if interference with impactful variant (function dependent?)
                pass
            matches["variants"][match]["uncertain"].append((variant, function))
    # Filter and return
    return find_best_match(sample, matches, functions, phased) 

def unpack_unphased_calling(unphased_calling, cont_graph):
    """Unpack a calling of a unphased sample into two phased callings."""
    # TODO why not stable?
    samples = sorted([s for s in unphased_calling.keys() if s[-1] != "+"])
    phased_calling = {sample: {'A': [], 'B': []} for sample in sorted(samples)} 
    for sample in samples:                
        # Group matches to find phasing TODO do this at step star_allele_calling
        # Assumes that every sample has one optimal core match QUESTION is this always valid?
        all_alleles = [allele for alleles in unphased_calling[sample] for allele in alleles]
        connected = nx.Graph()
        connected.add_nodes_from(all_alleles)
        for allele in all_alleles:
            for allele2 in all_alleles:
                if allele == allele2:
                    continue
                if cont_graph.has_edge(allele, allele2):
                    connected.add_edge(allele, allele2)
        groups = [c for c in nx.connected_components(connected)]
        # Split phasing by groups TODO do this neater
        n_cores = 0
        designate = {}
        for alleles in unphased_calling[sample]:
            n_cores_sub = 0
            for allele in alleles:
                if allele in designate:
                    continue
                if sort_types(allele) == 1: # Keep one core in each phasing
                    # Add core and grouped to only this phase
                    n_cores_sub += 1
                    if n_cores >= 2:
                        continue
                    phasing = "AB"[n_cores] 
                    designate[allele] = phasing
                    for group in groups: # Mark grouped as being in phase
                        if allele not in group:
                            continue
                        for other in group:
                            designate[other] = phasing
                    n_cores += 1
        # Add alleles to correct phasing based on grouping
        # Unknown groupings have to be placed in both
        # Keeps priority ranks determined earlier
        first = {'A': True, 'B': True}
        for alleles in unphased_calling[sample]:
            phased_calling[sample]['A'].append(set()) # Keep priority ranks
            phased_calling[sample]['B'].append(set())
            for allele in alleles:
                if allele in designate:
                    phased_calling[sample][designate[allele]][-1].add(allele)
                else:
                    phased_calling[sample]['A'][-1].add(allele)
                    phased_calling[sample]['B'][-1].add(allele)
            # Test ambiguity
            for phase in "AB":
                test_n_cores = [sort_types(a) for a in phased_calling[sample][phase][-1]].count(1) # Number of cores in this position
                if first[phase] and test_n_cores > 1: # Ambiguous first occurrence
                    # TODO how to handle   
                    warnings.warn(f"{sample}{phase} ambiguous calling: {phased_calling[sample][phase][-1]}")
                if test_n_cores > 0: # First occurrence of core
                    first[phase] = False

        if n_cores == 1: # One match found, possibly homozygous or one wildtype
            # Use double variants to find homozygous
            doubles = sample + '+'
            if doubles in unphased_calling.keys():
                for alleles in unphased_calling[doubles]:
                    for allele in alleles:
                        if sort_types(allele) != 1: # Must contain a core
                            continue
                        if allele == "CYP2D6*1": # Must not be the default TODO don't return this 
                            continue
                        # Found core in doubles, make this the second phasing
                        phased_calling[sample]['B'] = [set(al) for al in unphased_calling[doubles]] # TODO append/extend?
                        n_cores += 1
                        break
                    else: # Not found yet
                        continue
                    break # Found
        if n_cores == 1: # Not homozygous
                phased_calling[sample]['B'].append({"CYP2D6*1",}) # resort to wildtype
        if n_cores == 0:
            raise ValueError("No core alleles found in sample: {}".format(sample))
        # if sample == "HG00276":
        #     print(groups)
        #     print(unphased_calling[sample])
        #     print(phased_calling[sample])
        #     print(designate)
        #     exit()
    return phased_calling

def star_allele_calling_all(samples, nodes, edges, functions, supremals, phased=True):
    eq_graph = nx.Graph()
    eq_graph.add_nodes_from(nodes)
    eq_graph.add_edges_from([(left, right) for left, right, relation in edges if relation == va.Relation.EQUIVALENT])
    cont_graph = nx.DiGraph()
    cont_graph.add_nodes_from(nodes)
    cont_graph.add_edges_from([(left, right) for left, right, relation in edges if relation == va.Relation.IS_CONTAINED])
    """Iterate over samples and call star alleles for each."""
    if phased: # Alleles are treated separately
        callings = {sample[:-1]: {'A': None, 'B': None} for sample in sorted(samples)} 
        for sample in samples:
            calling = star_allele_calling(sample, eq_graph, cont_graph, functions, supremals, phased)
            if phased: # For phased samples each allele has a single calling
                sample_source, phasing = sample[:-1], sample[-1]
                callings[sample_source][phasing] = calling
    else: # Unphased samples have a single calling that has to be separated into two 
        callings = {sample: star_allele_calling(sample, eq_graph, cont_graph, functions, supremals, phased) for sample in samples} 
        # Split unphased into two callings (how/where?)  
        callings = unpack_unphased_calling(callings, cont_graph)         
    return callings

def matches_core(match):
    """Print the core allele from a match."""
    # WARNING not in use since getting via traversal is better
    if sort_types(match) == 1:
        return match
    elif sort_types(match) == 2:
        return match[:-4] # QUESTION is this the best way to get the core?
    else:
        raise Exception(f"Unexpected match type: {match}")
        # QUESTION needed to also handle variants?

def impact_position(supremal):
    """Check if a variant is silent based on positions of variants.
    
    This is a different approach to the other methods since it doesn't rely on online sources.
    Therefore it can be used on novel variants.
    Instead it uses the position of a variant to see if it can influence the protein.
    This is based on the intron-exon borders.
    """
    # TODO do earlier for all personal variants (and more?)
    # TODO find exons dynamically for any gene
    # Exon borders (one-based; closed end) 
    # Source? https://www.ncbi.nlm.nih.gov/genome/gdv/browser/gene/?id=1565
    # QUESTION are there more exons sometimes?
    # QUESTION Does alternative splicing occur?
    # QUESTION are these fixed?
    exons = [ 
        (42126499, 42126752),
        (42126851, 42126992),
        (42127447, 42127634),
        (42127842, 42127983),
        (42128174, 42128350),
        (42128784, 42128944),
        (42129033, 42129185),
        (42129738, 42129909),
        (42130612, 42130810)
    ]    
    # Length of splice site
    # Part of intron that is used for recognition by spliceosome
    # Source: https://www.ncbi.nlm.nih.gov/gene/1565
    # TODO get dynamically
    # QUESTION is this correct?
    splice_sites = ("GT", "AG") 
    # Functions for checking exon influence
    # Check if a range is entirely in an exon
    in_exon = lambda s, e: any(
        start <= s and e <= end # Interval is contained in exon
        for start, end in exons)
    # Check if a range overlaps with an exon 
    overlap_exon = lambda s, e: any(
        max(s, start) <= min(e, end) # Interval overlaps with exon (highest begin is lowe than lowest end)
        for start, end in exons) 
    # Check if a range overlaps with a splice site
    overlap_splice_site = lambda s, e: any(
        s <= (start - len(splice_sites[0]) <= e) or # interval contains left splice site
        s <= (end + len(splice_sites[1])) <= e
        for start, end in exons)
    # Check if a mutation disturbs triplets by itself
    # difference (deleted or inserted) between area of influence and sequence is not a multiple of 3 
    frameshift = lambda v: (abs(v.end - v.start - len(v.sequence)) % 3) != 0 # TODO is this correct?
    # Check if supremal representation of variant (covers different placements) can influence the exon
    start = supremal.start + 1 # One-based position
    end = supremal.end # Closed end position
    if overlap_splice_site(start, end):
        return "possible splice defect" # QUESTION is this correct or can splice sites occur in a different way as well?
    if in_exon(start, end) or overlap_exon(start, end): 
        if frameshift(supremal): 
            return "possible frameshift" # QUESTION is this correct
        # TODO split further into synonymous, early stop and missense?
        return "possible missense" 
    return None # Intronic and thus certainly silent TODO correct value?

def is_noise(variant, functions, supremal): # TODO rename to relevance?
    """Determine if a variant is noise.
    
    Noise is defined as not being relevant for calling.
    The variant is noise if it has no impact on the protein through non-synonymous mutations or splicing effects.

    This approach uses the online annotations but falls back on a sequence based approach.
    """
    if variant in functions: function = functions[variant]  # PharmVar variant 
    else: function = impact_position(supremal) # Personal variant
    if function == '': function = None # QUESTION what is the correct interpretation of this (no change or unknown)?
    # Check the impact of the variant
    if function == None: # QUESTION what is the correct interpretation (no change or unknown)?
        return True, None
    elif 'splice defect' in function: # Change in expression
        return False, function
    else: # Explicit change on protein level
        return False, function # Missense mutation

def print_classification(classifications, detail_level=0):
    """Print the classification of samples.
    
    Different detail levels are available.
    0: Only print best core match
    1: Simplify to core matches
    2: Print all matches
    """
    for sample, classification in classifications.items():
        selected_alleles = {'A': [], 'B': []}
        for phase, ordered_alleles in classification.items():
            for alleles in ordered_alleles:
                for allele in alleles: # TODO do this nicer
                    if detail_level in (0, 1) and sort_types(allele) != 1: # Only core
                        continue
                    selected_alleles[phase].append(allele)
                    if detail_level == 0: # Only display best
                        break
                else: # No break 
                    continue
                break
        # Print
        print(f"{sample}: {','.join(selected_alleles['A'])}/{','.join(selected_alleles['B'])}")
