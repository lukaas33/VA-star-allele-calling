import algebra as va
import networkx as nx
import warnings
from .data import api_get
import re
import copy

all_functions = ['unknown function', 'uncertain function', 'normal function', 'decreased function', 'no function']

def sort_function(f):
    """Sort function annotation based on severity."""
    # TODO use enums
    # QUESTION: what is the difference between uncertain and unknown?
    if f not in all_functions:
        raise Exception("Unknown function: " + f)
    functions = all_functions[3:]
    if f not in functions: 
        return 0 
    return functions.index(f) + 1

def sort_types(v):
    """Sort variant types in order of specificity (somewhat arbitrary)."""
    # TODO use enums or types for the data
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
    """Recursively find the contained (and equivalent) variants from a given start node."""
    if start in visited: return # Already visited
    visited.add(start) # Needed to avoid equivalence loop and to avoid doubles
    if sort_types(start) in find: # Need to find this type
        matches.add(start)
    if stop is not None and sort_types(start) == stop: # Don't look further here (terminal nodes)
        return 
    # Find equivalent (needed for finding all suballeles or all variants)
    if start in eq_graph.nodes():
        for match in eq_graph[start]: 
            if sort_types(match) not in (1, 2, 3): # TODO allow iteration over own personal variants?
                continue
            find_contained_variants(match, cont_graph, eq_graph, matches, visited, find, stop) # Add equivalents and maybe traverse
    # Find contained
    if start in cont_graph.nodes(): 
        for match, _ in cont_graph.in_edges(start):
            if sort_types(match) not in (1, 2, 3): # TODO allow iteration over own personal variants?
                continue
            find_contained_variants(match, cont_graph, eq_graph, matches, visited, find, stop) # Add contained and maybe traverse

def find_overlapping_variants(start, cont_graph, eq_graph, overlap_graph, matches, visited, find, stop=None):
    """Recursively find the overlapping variants from a given start node."""
    if start in visited: return # Already visited
    visited.add(start) # Needed to avoid equivalence loop and to avoid doubles
    if stop is not None and sort_types(start) == stop: # Don't look further here (terminal nodes)
        return 
    # Find equivalent not needed (since everything is connected to core)
    # Find contained
    if start in cont_graph.nodes(): 
        for match, _ in cont_graph.in_edges(start):
            if sort_types(match) not in (1, 2, 3): # TODO allow iteration over own personal variants?
                continue
            find_overlapping_variants(match, cont_graph, eq_graph, overlap_graph, matches, visited, find, stop) # Don't add but traverse
    # Find overlapping
    if start in overlap_graph.nodes():
        for match in overlap_graph[start]:
            if sort_types(match) not in (1, 2, 3): # TODO allow iteration over own personal variants?
                continue
            matches.add(match) # Add overlapping, don't traverse further

def find_most_specific(matches, cont_graph):
    """Find the most specific alleles from a list of contained alleles."""  
    reduced = set()
    for node1 in matches:
        for node2 in matches:
            if node1 == node2: # Skip self 
                continue
            if node2 not in cont_graph.nodes(): # Skip alleles not in graph (in practice only *1)
                continue
            if sort_types(node1) == 1 and sort_types(node2) == 2: # Skip containment of core in suballeles as this is expected
                continue
            if node1 not in nx.ancestors(cont_graph, node2): # Node1 is contained in node2
                continue
            # No issues 
            break
        else: # No break
            reduced.add(node1) # Keep match
    return reduced

def prioritize_calling(sample, matches, functions, phased):
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
        annotated = {f: set([m for m in matches["contained"] if functions[m] == f]) for f in all_functions}
        # If the annotation is no function, this allele is most disruptive and is chosen
        if len(annotated["no function"]) > 0:	
            sorted_matches.append(annotated["no function"])
        # If uncertain or unknown functions are present, the ordering after no function cannot be determined
        # An unknown function could be no function and thus no allele should be preferred.
        if len(annotated["unknown function"]) > 0 or len(annotated["uncertain function"]) > 0:
            # Put all in the same position of ordering
            sorted_matches.append(
                annotated["decreased function"] |
                annotated["normal function"] |
                annotated["unknown function"] |
                annotated["uncertain function"]
            )
        else: # Ordering can be determined between normal and decreased function
            # More disruptive alleles are chosen first
            if len(annotated["decreased function"]) > 0:
                sorted_matches.append(annotated["decreased function"])
            if len(annotated["normal function"]) > 0:
                sorted_matches.append(annotated["normal function"])
    # Check certainty and consistency of matches
    # for match in sorted_matches:
    #     n_cores = 0
    #     for allele in match:
    #         if sort_types(allele) == 1: # Coreallele
    #             n_cores += 1
    #         # check certainty
    #         if len(matches['variants'][allele]['uncertain']) > 0:
    #             # TODO how to express this to the user?
    #             # warnings.warn(f"{sample}={allele}: Uncertain variants found: {matches['variants'][allele]['uncertain']}.")
    #             pass
    #     if n_cores > 1 and phased or n_cores > 2 and not phased: # Ambiguous phased/unphased
    #         warnings.warn(f"{sample}: Multiple corealleles found with the same priority: {match}.")
    #     elif n_cores == 0: # No found yet
    #         continue
    #     else: # Unambiguous
    #         break # No need to check further
    # else: # No cores found
    # Catch-all is the default allele and is in theory always present
    # *1 will never be in contained so it is not prioritized over others
    # TODO do this somewhere else?
    sorted_matches.append({"CYP2D6*1",})        
    return sorted_matches

def star_allele_calling(sample, eq_graph, cont_graph, overlap_graph, functions, supremals, phased):
    """Determine star allele calling for a sample based on va relations.
    
    Find based on pruned graph containing relations between samples and alleles and between themselves.
    """
    # QUESTION: is it needed to look at suballeles for calling?
    # QUESTION: is it needed to look at individual variants for calling?
    # All information of a calling
    matches = {
        "equivalent": set(), # Equivalent alleles
        "contained": set(), # Contained alleles
        "overlap": set(), # Overlapping alleles
        "variants": {
            # extra variants for each allele
        }
    }
    # Find equivalent alleles
    if sample in eq_graph.nodes(): matches["equivalent"] = set([match for match in eq_graph[sample] if sort_types(match) in (1, 2)])
    # Find contained
    find_contained_variants(sample, cont_graph, eq_graph, matches["contained"], set(), (1,2))
    matches["contained"] = set([m for m in matches["contained"] if m not in matches["equivalent"]]) # Remove equivalent from here
    # Find overlapping alleles
    find_overlapping_variants(sample, cont_graph, eq_graph, overlap_graph, matches["overlap"], set(), (1,2))
    matches["contained"] |= matches["overlap"] # Add overlap to contained TODO change to lower priority?

    # Check certainty of calling for each matched allele
    # TODO change extra variant definition to be more strict
    # vs = set() # All variants indirectly or directly contained in the sample 
    # find_contained_variants(sample, cont_graph, eq_graph, vs, set(), (3, 5)) # TODO get from PharmVar?
    # for match in matches["contained"] | matches["equivalent"]: # Equivalent added but these will be empty
    #     matches["variants"][match] = {"noise": [], "uncertain": []}
    #     # Find 'extra' variants: in the sample but not in the called allele
    #     ws = set() # All variants indirectly or directly contained in the sample
    #     find_contained_variants(match, cont_graph, eq_graph, ws, set(), (3, 5)) # TODO get from PharmVar?
    #     extra_variants = []
    #     for v in vs: # In sample
    #         if v in ws: # But not in this allele
    #             continue
    #         extra_variants.append(v) 
    #     # Check if the extra variants can influence the calling
    #     for variant in extra_variants:
    #         # Check what the impact of the variant is
    #         noise, function = is_noise(variant, functions, supremals[variant])
    #         if noise: # No evidence of relevance
    #             matches["variants"][match]["noise"].append(variant) # TODO annotate noise
    #             continue # TODO need to check interference here?
    #         # Evidence of possible protein impact
    #         # Check if the impact could interfere with the calling in some way
    #         if "splice defect" in function: # Splice defect, outcome uncertain
    #             pass
    #         elif "fs" == function[-2:]: # Frameshift
    #             # TODO check if effect downstream
    #             pass
    #         elif "X" == function[-1]: # Stop codon
    #             # TODO check if effect downstream
    #             pass
    #         else: # Missense
    #             # TODO check if interference with impactful variant (function dependent?)
    #             pass
    #         matches["variants"][match]["uncertain"].append((variant, function))
    # Filter and return
    return prioritize_calling(sample, matches, functions, phased) 

def unpack_unphased_calling(unphased_calling, cont_graph):
    """Unpack a calling of a unphased sample into two phased callings."""
    raise Exception("Not implemented yet")
    samples = sorted([s for s in unphased_calling.keys() if s[-1] != "+"])
    phased_calling = {sample: {'A': [], 'B': []} for sample in sorted(samples)} 
    for sample in samples:                
        n_cores = 0
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
        # Add alleles to phasing based on grouping
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
            # TODO what is the ordering between these steps?
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
        # if sample == "NA19207":
        #     print(unphased_calling[sample])
        #     print(phased_calling[sample])
        #     print(groups)
        #     print(designate)
        #     exit()
    return phased_calling

def star_allele_calling_all(samples, nodes, edges, functions, supremals, phased=True, detail_level=0):
    eq_graph = nx.Graph()
    eq_graph.add_nodes_from(nodes)
    eq_graph.add_edges_from([(left, right) for left, right, relation in edges if relation == va.Relation.EQUIVALENT])
    cont_graph = nx.DiGraph()
    cont_graph.add_nodes_from(nodes)
    cont_graph.add_edges_from([(left, right) for left, right, relation in edges if relation == va.Relation.IS_CONTAINED])
    overlap_graph = nx.DiGraph()
    overlap_graph.add_nodes_from(nodes)
    overlap_graph.add_edges_from([(left, right) for left, right, relation in edges if relation == va.Relation.OVERLAP])
    """Iterate over samples and call star alleles for each."""
    if phased: # Alleles are treated separately
        callings = {sample[:-1]: {'A': None, 'B': None} for sample in sorted(samples)} 
        for sample in samples:
            calling = star_allele_calling(sample, eq_graph, cont_graph, overlap_graph, functions, supremals, phased)
            sample_source, phasing = sample[:-1], sample[-1]
            callings[sample_source][phasing] = calling
    else: # Unphased samples have a single calling that has to be separated into two 
        raise NotImplementedError()
        callings = {sample: star_allele_calling(sample, eq_graph, cont_graph, functions, supremals, phased) for sample in samples} 
        # Split unphased into two callings (how/where?)  
        callings = unpack_unphased_calling(callings, cont_graph)         
    return calling_to_repr(callings, cont_graph, detail_level)

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
    # TODO retire this function
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
    # TODO use new method
    if variant in functions: function = functions[variant]  # PharmVar variant 
    else: function = impact_position(supremal) # Personal variant TODO don't do this
    if function == '': function = None # QUESTION what is the correct interpretation of this (no change or unknown)?
    # Check the impact of the variant
    if function == None: # QUESTION what is the correct interpretation (no change or unknown)?
        return True, None
    elif 'splice defect' in function: # Change in expression
        return False, function
    else: # Explicit change on protein level
        return False, function # Missense mutation

def calling_to_repr(callings, cont_graph, detail_level=0):
    """Change the calling to a representation of a certain amount of detail.
    
    Different detail levels are available.
    0: Only print best core match(es)
    1: Simplify to core matches only
    2: Print core and suballele matches
    3: Print all matches and their less specific contained alleles
    """
    # TODO only print *1 if other present
    if detail_level != 3: # Remove matches contained in others
        for sample in callings:
            for phase in "AB":
                all_alleles = set([a for al in callings[sample][phase] for a in al if al != None])
                reduced_alleles = find_most_specific(all_alleles, cont_graph)
                # print(sample, phase, callings[sample][phase])
                for al in callings[sample][phase]:
                    al &= reduced_alleles
                # print(sample, phase, callings[sample][phase])
                # print()
        # exit()
    selected_calling = {}
    for sample, calling in callings.items():
        selected_alleles = {'A': [], 'B': []} # Flatten into ordered list
        for phase, ordered_alleles in calling.items():
            for alleles in ordered_alleles:
                for allele in alleles:
                    if detail_level in (0, 1) and sort_types(allele) != 1: # Only core
                        continue
                    selected_alleles[phase].append(allele)
                    if detail_level == 0: # Only display best
                        break
                else: # No break 
                    continue
                break
        selected_calling[sample] = selected_alleles
    return selected_calling
        
def find_path(s, t, cont_graph, eq_graph, overlap_graph, path=None, visited=None):
    """Find a path from s to t in the graphs"""
    # TODO integrate in web interface
    if path is None: path = [s]
    if visited is None: visited = set()
    if s in visited: return None
    visited.add(s)
    if s == t:
        return path
    for g in (cont_graph, eq_graph, overlap_graph):
        if s in g.nodes():
            for n in g[s]:
                if sort_types(n) in (4, 5): # Don't use samples for iteration
                    continue
                path.append(n)
                result = find_path(n, t, cont_graph, eq_graph, overlap_graph, path, visited)
                if result is not None:
                    return result
                path.pop()
    return None