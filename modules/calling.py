import algebra as va
import networkx as nx
import warnings
from .data import api_get
from .other_sources import find_id_hgvs, get_annotation_entrez
import re
import copy

all_functions = ['unknown function', 'uncertain function', 'normal function', 'decreased function', 'no function']

def sort_function(f):
    """Sort function annotation based on severity."""
    raise DeprecationWarning("This function is redundant now")
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
    raise DeprecationWarning("This function is not used anymore")
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

def find_overlapping_variants(current, cont_graph, overlap_graph, find):
    """Recursively find the overlapping variants from a given current node."""
    # Find overlapping
    matches = set()
    if current in overlap_graph.nodes():
        for match in overlap_graph[current]:
            if sort_types(match) not in find: # Looking for this type
                continue
            matches.add(match) # Add overlapping, don't traverse this node further (direct overlap)
    # Iterate over equivalent not needed (since everything is connected to core)
    # Iterate over contained to find most specific overlaps
    if current in cont_graph.nodes(): 
        for node, _ in cont_graph.in_edges(current): # Go over variants contained in this
            if sort_types(node) == 4: # skip samples
                continue
            # Find more specific overlaps
            for n in find_overlapping_variants(node, cont_graph, overlap_graph, find):
                matches.add(n)
    return matches

def prioritize_calling(matches, functions):
    """Prioritize matches with the same rank based on functional annotation.

    Returns list of sets where each index is the priority level and may contain several matches.
    """
    # Find annotation for all matches
    prioritized_matching = []
    annotated = {f: set([m for m in matches if functions[m] == f]) for f in all_functions}
    # If the annotation is no function, this allele is most disruptive and is prioritized over everything
    if len(annotated["no function"]) > 0:	
        prioritized_matching.append(annotated["no function"])
    # If uncertain or unknown functions are present, the ordering after no function cannot be determined
    # An unknown function could be as bad as no function.
    if len(annotated["unknown function"]) > 0 or len(annotated["uncertain function"]) > 0:
        # Put all in the same position of ordering
        prioritized_matching.append(
            annotated["decreased function"] |
            annotated["normal function"] |
            annotated["unknown function"] |
            annotated["uncertain function"]
        )
    else: # No uncertain or unknown functions
        # Ordering can be determined between normal and decreased function when no uncertain or unknown functions are present
        # More disruptive alleles are chosen first
        if len(annotated["decreased function"]) > 0:
            prioritized_matching.append(annotated["decreased function"])
        if len(annotated["normal function"]) > 0:
            prioritized_matching.append(annotated["normal function"])
    return prioritized_matching

def star_allele_calling(sample, eq_graph, cont_graph, overlap_graph, functions, supremals, reference):
    """Determine star allele calling for a sample based on va relations.
    
    Find matches based on pruned graph.
    Graph must contain relations between all samples and PharmVar alleles.
    Graph must also contain relations between all PharmVar alleles.

    Returns a list of allele matches ordered by the strength of the relation.
    """
    # QUESTION does this method work directly for unphased data?
    # QUESTION: is it needed to look at suballeles for calling?
    # QUESTION: is it needed to look at individual variants for calling?
    # Don't predict unparsable samples
    # Would be predicted as *1 otherwise
    # Empty alleles are in the supremal dictionary with None as a value
    if sample not in supremals: 
        return [{'CYP2D6*?',}]
    # All information of a calling
    start = [sample]
    matches = []
    # Find equivalent allele if present
    if sample in eq_graph.nodes() and len(eq_graph[sample]) != 0:
        m_eq = [match for match in eq_graph[sample] if sort_types(match) in (1, 2)]
        if len(m_eq) > 1: # More equivalents found, not possible for correct dataset
            raise Exception(f"{sample}: This should not happen, multiple equivalent matches found: {matches['equivalent']}.")
        if len(m_eq) == 1: 
            matches.append(set(m_eq)) # Strongest match
            start.append(m_eq[0]) # Also look for contained/overlapping from equivalent allele
    # Find directly contained alleles
    m_cont = set()
    for eq in start:
        if eq in cont_graph.nodes() and len(cont_graph.in_edges(eq)) != 0:
            m_cont |= set([m for m, _ in cont_graph.in_edges(eq) if sort_types(m) in (1, 2)])
    if len(m_cont) > 0: matches.append(m_cont) # Less strong match
    # Find directly overlapping alleles
    # also need to check contained alleles because of most specific pruning
    m_ov = set()
    for eq in start:
        m_ov |= find_overlapping_variants(eq, cont_graph, overlap_graph, (1, 2))
    if len(m_ov) > 0: matches.append(m_ov) # Least strong match (can treat the same as contained?)
    # Add default allele, will have the lowest priority
    matches.append({"CYP2D6*1",})

    # # TODO move elsewhere
    # # Find extra variants to check if they can influence the calling
    # # TODO finalise based on notes
    # # TODO visualise
    # # Extra variants are variants that are in the sample but not in any of the called alleles
    # variants = []
    # if sample in cont_graph.nodes(): variants += [m for m, _ in cont_graph.in_edges(sample) if sort_types(m) in (3, 5)]
    # if sample in overlap_graph.nodes(): variants += [m for m in overlap_graph[sample] if sort_types(m) in (3, 5)]
    # # Check relevance of variants
    # alleles = matches["equivalent"] | matches["contained"] | matches["overlapping"]
    # for variant in variants:
    #     if is_relevant(variant, functions, supremals, alleles, reference):
    #         matches["variants"].add(variant)

    # Return all matches to be filtered later (based on detail level)
    return matches

def separate_callings(unphased_calling):
    """Unpack a calling of a unphased sample into two phased callings."""
    samples = list(unphased_calling.keys())
    phased_calling = {sample: {'A': [], 'B': []} for sample in sorted(samples)} 
    # Add phasing based on homozygosity
    for sample in samples:
        for alleles in unphased_calling[sample]['all']:
            phased_calling[sample]['A'].append(set()) # maintain priority ranks
            phased_calling[sample]['B'].append(set())
            for allele in alleles:
                # Check if match was found using homozygous variants
                if all([allele not in alls for alls in unphased_calling[sample]['hom']]):
                    continue
                # Match is phased
                phased_calling[sample]['A'][-1].add(allele)
                phased_calling[sample]['B'][-1].add(allele)
        if sample == "HG00373":
            print(phased_calling[sample])
            print(unphased_calling[sample])
            exit()
    raise NotImplementedError("Not implemented yet")
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
    return phased_calling

def star_allele_calling_all(samples, nodes, edges, functions, supremals, reference, phased=True, detail_level=0):
    """Iterate over samples and call star alleles for each."""
    eq_graph = nx.Graph()
    eq_graph.add_nodes_from(nodes)
    eq_graph.add_edges_from([(left, right) for left, right, relation in edges if relation == va.Relation.EQUIVALENT])
    cont_graph = nx.DiGraph()
    cont_graph.add_nodes_from(nodes)
    cont_graph.add_edges_from([(left, right) for left, right, relation in edges if relation == va.Relation.IS_CONTAINED])
    overlap_graph = nx.Graph()
    overlap_graph.add_nodes_from(nodes)
    overlap_graph.add_edges_from([(left, right) for left, right, relation in edges if relation == va.Relation.OVERLAP])
    callings = {sample.split('_')[0]: {} for sample in sorted(samples)} 
    for sample in samples:
        calling = star_allele_calling(sample, eq_graph, cont_graph, overlap_graph, functions, supremals, reference)
        sample_source, phasing = sample.split('_')
        callings[sample_source][phasing] = calling
    if not phased: # Single calling contains multiple matches which must be separated
        callings = separate_callings(callings)
    # Create a textual representation of the calling based on the amount of detail needed
    representation = calling_to_repr(callings, cont_graph, functions, **detail_from_level(detail_level))
    return representation

def find_core_string(match):
    """Get the core allele from a match string wise."""
    raise DeprecationWarning("not in use since getting via traversal is better")
    if sort_types(match) == 1:
        return match
    elif sort_types(match) == 2:
        return match[:-4] # QUESTION is this the best way to get the core?
    else:
        raise Exception(f"Unexpected match type: {match}")
        # QUESTION needed to also handle variants?

def find_core_traversal(match, cont_graph):
    """Get the coreallele(s) from a match based on traversal.
    
    Returns a set since a suballele can contain more corealleles.
    Will not return indirectly contained corealleles.
    This won't find the corealleles that are not contained in the suballele by design.
    """
    cores = []
    if match in cont_graph.nodes(): # Only have to check contained since everything is connected to the core
        for node, _ in cont_graph.in_edges(match):
            if sort_types(node) == 1: # Found core contained in match, don't traverse further
                cores.append(node)
            elif sort_types(node) == 2: # Traverse further
                cores.extend(find_core_traversal(node, cont_graph))
    return cores

def impact_position(supremal):
    """Check if a variant is silent based on positions of variants.
    
    This is a different approach to the other methods since it doesn't rely on online sources.
    Therefore it can be used on novel variants.
    Instead it uses the position of a variant to see if it can influence the protein.
    This is based on the intron-exon borders.
    """
    raise DeprecationWarning("Will not be finished since annotations are used")
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

def is_relevant(variant, functions, supremals, alleles, reference): 
    """Determine if a variant is relevant.
    
    Relevance is defined as being relevant for calling.
    It is not complete but tries to filter out variants that are definitely not relevant.
    These variants can then be ignored in visualisations
    This approach uses the online annotations of variants.
    """
    raise NotImplementedError("Not finished")
    # QUESTION is this method valid?
    # QUESTION why does id lookup at 42132027>T result in unparsable HGVS like NC_000022.11:g.42132044_42132047T[4]CTTTT[1]? 
    # Check if variant possibly interferes with any allele (overlaps with supremal)
    interference = False
    for allele in alleles:
        if max(supremals[allele].start, supremals[variant].start) <= min(supremals[allele].end, supremals[variant].end): # Overlap between variant and allele
            # Variant can alter the calling since it can disturb the allele's definition.
            interference = True
    # Check impact of variant individually
    # Annotations of PharmVar variants are retrieved from the PharmVar database
    # Annotations of personal variants retrieved from online sources
    impact = functions[variant]
    # Handle impact annotations
    if impact is None: # Also includes None
        # TODO check if this is correct
        return True # No annotation available, assume relevant
    if "fs" in impact or "frameshift" in impact: # Frameshift is always relevant
        return True
    if "splice" in impact: # Splice defect is always relevant
        return True
    if "stop_gained" in impact or re.match(r"[A-Z][0-9]*X", impact): # Stop gained is always relevant
        return True
    # Otherwise only relevant if it interferes with an allele
    # Other impact annotations include: (non-)synonymous, missense/[A-Z][0-9]*[A-Z], inframe deletion/insertion, intron variant and coding sequence variant
    # Since we are not predicting the protein sequence we assume the worse for them if they interfere with a allele
    # QUESTION is this correct, can this be more precise?
    # TODO handle all cases?
    return interference 

def detail_from_level(level):
    """Return a dictionary of options for a specific detail level
    
    The output corresponds to the options of calling_to_repr.

    Different detail levels are available.
    0: Only print best core match(es)
    1: Print all direct core matches
    TODO implement more levels
    """
    kwargs = {}
    kwargs["find_cores"] = True 
    kwargs["prioritize_function"] = level <= 0
    kwargs["prioritize_strength"] = level <= 0
    kwargs["suballeles"] = not (level <= 1)
    kwargs["default"] = not (level <= 1)
    return kwargs

def calling_to_repr(callings, cont_graph, functions, find_cores, suballeles, default, prioritize_function, prioritize_strength):
    """Change the calling to a representation of a certain amount of detail.
    
    Each calling contains all direct matches which can be suballeles or core alleles.
    The amount of detail in a representation can be made more or less specific.
    The representation will show an list of matches ordered on start allele number.

    These are the options:
    find_cores: find the corealleles of each suballele
    suballeles: show suballeles or not
    default: show the default allele always or only if no other allele is found
    prioritize_function: prioritize alleles based on functional annotation.
    prioritize_strength: prioritize alleles based on relation strength.
    """
    def star_num(match):
        if '?' in match:
            return 0
        else:
            return int(match.split('*')[1])
    def remove_contained(matches):
        """Filter out cores that are contained in other cores"""
        # TODO do this by not adding some cores in the first place
        for i1, match1 in enumerate(matches):  
            for i2, match2 in enumerate(matches):
                if i1 == i2:
                    continue
                if match1 == match2: # Duplicates
                    matches.remove(match1)
                if sort_types(match2) == 2: # Cores can be contained in subs
                    continue
                if match1 in nx.ancestors(cont_graph, match2): # match1 is contained in match2
                    matches.remove(match1)
                    break
        return matches
    # TODO add indirectly found alleles
    # TODO add alternative descriptions
    representation = {}
    for sample in callings:
        representation[sample] = {} # Keep structure
        for phase in callings[sample]:
            representation[sample][phase] = [] # Flatten into single list
            for alleles in callings[sample][phase]:
                for allele in alleles:
                    # Add detail
                    if find_cores and sort_types(allele) == 2: # Find cores of suballeles
                        cores = find_core_traversal(allele, cont_graph) # Cores of suballele
                        representation[sample][phase].extend(cores)
                    # Remove detail
                    if not suballeles and sort_types(allele) == 2: # Don't represent suballeles
                        continue
                    if not default and allele == "CYP2D6*1": # Don't represent default allele if there are other matches
                        if len(representation[sample][phase]) > 0: # Works since default is always last
                            continue
                    representation[sample][phase].append(allele) # Add match
                if prioritize_strength and len(representation[sample][phase]) > 0: # Only keep strongest related alleles
                    break
            # Remove cores contained in other cores 
            # These can occur because of get_core_traversal
            representation[sample][phase] = remove_contained(representation[sample][phase])
            # Prioritize alleles of equal rank (after filtering)
            if prioritize_function and len(representation[sample][phase]) > 1:
                prioritized = prioritize_calling(representation[sample][phase], functions)
                representation[sample][phase] = list(prioritized[0]) # Only keep the first priority
            # Sort alleles based on star allele number
            representation[sample][phase].sort(key=star_num)
        # Sort phases based on start allele number
        low_A = min([star_num(m) for m in representation[sample]['A']])
        low_B = min([star_num(m) for m in representation[sample]['B']])
        if low_A > low_B:
            representation[sample]['A'], representation[sample]['B'] = representation[sample]['B'], representation[sample]['A']
    return representation
        
def find_path(s, t, cont_graph, eq_graph, overlap_graph, path=None, visited=None):
    """Find a path from s to t in the graphs"""
    raise NotImplementedError("Implementation contains bugs")
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