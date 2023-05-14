import algebra as va
import networkx as nx
import warnings
from .data import api_get
from .other_sources import find_id_hgvs, get_annotation_entrez, severity_GO, severity_pharmvar
import re
import copy

all_functions = ['unknown function', 'uncertain function', 'normal function', 'decreased function', 'no function']
def sort_function(f):
    """Sort function annotation based on severity."""
    raise DeprecationWarning("This function is redundant now")
    # TODO handle function not assigned (doesn't occur in calling)
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
    raise DeprecationWarning("This function is not used any more")
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
    raise DeprecationWarning("This function is not used any more")
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
    # Don't predict unparsable samples
    # Would be predicted as *1 otherwise
    # Empty alleles are in the supremal dictionary with None as a value
    if sample not in supremals: 
        return [{'CYP2D6*?',}]
    # All information of a calling
    matches = []
    # Find equivalent allele if present
    if sample in eq_graph.nodes():
        m_eq = set([match for match in eq_graph[sample] if sort_types(match) in (1, 2)])
        if len(m_eq) > 1: # More equivalents found, not possible for correct dataset
            raise Exception(f"{sample}: multiple equivalent matches found: {matches['equivalent']}.")
        if len(m_eq) == 1: 
            matches.append(m_eq) # Strongest match
    # Find directly contained alleles
    if sample in cont_graph.nodes():
        m_cont = set([m for m, _ in cont_graph.in_edges(sample) if sort_types(m) in (1, 2)])
        if len(m_cont) > 0: 
            matches.append(m_cont) # Less strong match than equivalent
    # Find directly overlapping alleles
    # not looking at overlaps that are a result of contained alleles
    if sample in overlap_graph.nodes():
        m_ov = set([m for m in overlap_graph[sample] if sort_types(m) in (1, 2)])
        if len(m_ov) > 0: 
            if len(matches) > 0:
                # Overlap can be treated the same since there are no occurrences
                matches[-1] |= m_ov # Least strong match 
            else: # No equivalent or contained alleles
                matches.append(m_ov)
    # Add default allele, will have the lowest priority
    matches.append({"CYP2D6*1",})
    # Return all matches to be filtered later (based on detail level)
    return matches

def separate_callings(unphased_calling, cont_graph, functions):
    """Unpack a calling of a unphased sample into two phased callings."""
    samples = list(unphased_calling.keys())
    phased_calling = {sample: {'A': [], 'B': []} for sample in sorted(samples)} 
    # Infer phasing for some samples
    for sample in samples:
        for alleles in unphased_calling[sample]['all']:
            phased_calling[sample]['A'].append(set()) # Maintain relation strength rank
            phased_calling[sample]['B'].append(set())
            # Handle unparsable samples
            if alleles == {'CYP2D6*?',}:
                phased_calling[sample]['A'][-1].add("CYP2D6*?")
                phased_calling[sample]['B'][-1].add("CYP2D6*?")
                break
            # Add largest homozygous matches to both callings
            # Don't guess the heterozygous/mixed matches but add them all to A
            # This will be handled by alternative callings later
            for allele in alleles: # is largest instead of overlap
                if find_core_string(allele) == "CYP2D6*1": # Skip default, is added last
                    continue
                phased_calling[sample]['A'][-1].add(allele)
                if not any([allele in alls for alls in unphased_calling[sample]['hom']]): # check if is homozygous
                    continue
                phased_calling[sample]['B'][-1].add(allele)
        else: # No break
            for phase in "AB":
                # remove empty
                phased_calling[sample][phase] = [c for c in phased_calling[sample][phase] if len(c) > 0] 
                # Add default allele as last priority
                phased_calling[sample][phase].append({"CYP2D6*1",})
    return phased_calling

def valid_calling(calling, cont_graph, homozygous):
    """Check if a calling is valid based on homozygous contained alleles."""
    # Find ancestors for both phases
    # TODO make more efficient by not generating all ancestors
    ancestors = {}
    for phase in 'AB':
        ancestors[phase] = set()
        for alleles in calling[phase][:-1]: # All non-default matches for a phase
            for allele in alleles: 
                if allele not in cont_graph.nodes():
                    continue
                for ancestor in nx.ancestors(cont_graph, allele): # All alleles that are contained in this allele
                    if sort_types(ancestor) in (1, 2):
                        ancestors[phase].add(ancestor)
    # If this calling is valid the overlap between the ancestors of both phases should be equal to the homozygous alleles
    overlap = ancestors['A'] & ancestors['B']
    if overlap == homozygous: 
        return True
    return False

def generate_alternative_callings(calling, cont_graph, extended, depth=1):
    """Generate alternative callings for a sample.
    
    This is useful for unphased not homozygous matches.
    """
    # TODO reduce runtime by not generating all possible callings or not copying
    # Need to copy to prevent changing the original calling
    # would not be needed in case of one iteration over this generator
    # All in phase A already a valid answer
    yield copy.deepcopy(calling) 
    for i, alleles in enumerate(calling['A'][:-1]): # All non-default matches for this sample
        for allele in list(alleles): # Maintain relation strength rank
            # Find contained alleles of this allele
            if allele in cont_graph.nodes() and allele not in extended: 
                # TODO check overlap graph
                deeper_matches = set()
                for underlying, _ in cont_graph.in_edges(allele): # Find underlying alleles
                    if sort_types(underlying) in (1, 2): 
                        deeper_matches.add(underlying)
            # Replace allele with underlying alleles
            if len(deeper_matches) > 0:
                deeper_calling = copy.deepcopy(calling)
                deeper_calling['A'][i] |= deeper_matches
                # print(depth, "replace", allele, "with", deeper_matches, "in A:", deeper_calling)
                # Recurse with extended allelles
                # Track which alleles have been extended to prevent duplicates
                # TODO also add replaced
                for alternative in generate_alternative_callings(deeper_calling, cont_graph, extended | {allele,},  depth+1):
                    yield alternative
            # Move alleles to other phase
            calling['A'][i].remove(allele)
            if len(calling['A'][i]) == 0: continue # Don't move since this would result in a duplicate (X/Y = Y/X)
            while len(calling['B']) < len(calling['A']): calling['B'].insert(0, set()) # Make B the same length as A
            calling['B'][i].add(allele)
            # print(depth, "move", allele, "to B:", calling)
            # Also a valid answer
            yield copy.deepcopy(calling)


def star_allele_calling_all(samples, nodes, edges, functions, supremals, reference, phased=True, detail_level=0):
    """Iterate over samples and call star alleles for each."""
    eq_graph = nx.Graph([(left, right) for left, right, relation in edges if relation == va.Relation.EQUIVALENT])
    cont_graph = nx.DiGraph([(left, right) for left, right, relation in edges if relation == va.Relation.IS_CONTAINED])
    overlap_graph = nx.Graph([(left, right) for left, right, relation in edges if relation == va.Relation.OVERLAP])
    # Call each sample
    callings = {sample.split('_')[0]: {} for sample in sorted(samples)} 
    for sample in samples:
        calling = star_allele_calling(sample, eq_graph, cont_graph, overlap_graph, functions, supremals, reference)
        sample_source, phasing = sample.split('_')
        callings[sample_source][phasing] = calling
    # Create a textual representation of the calling based on the amount of detail needed
    if phased: # Calling is phased
        representations = {sample: calling_to_repr(callings[sample], cont_graph, functions, **detail_from_level(detail_level)) for sample in callings}
        return representations
    if not phased: # Unphased calling should be separated
        sep_callings = separate_callings(callings, cont_graph, functions)
        representations = {}
        # TODO move this to a filter alternatives function
        # TODO allow for keeping multiple alternative representations
        # TODO check NA07348
        for sample, calling in sep_callings.items():
            if sample != "NA12815":
                continue
            if calling['A'] == calling['B']: # Already phased (homozygous)
                # TODO is this a good check for homozygous?
                representations[sample] = calling_to_repr(calling, cont_graph, functions, **detail_from_level(detail_level))
                continue
            # Generate alternative callings
            alternatives = generate_alternative_callings(calling, cont_graph, set())
            # TODO ordering
            # All homozygous alleles for the current sample
            homozygous = set([allele for alleles in callings[sample]['hom'] for allele in alleles if allele != "CYP2D6*1"])
            preferred = None # Using split method, *1/*x+... or *x/*y
            for alternative in alternatives:
                # Filter out non valid callings
                # TODO filter earlier at generation (not allowing move in in case of non-homozygous overlap)
                if not valid_calling(alternative, cont_graph, homozygous):
                    pass
                    # continue
                representation = calling_to_repr(alternative, cont_graph, functions, **detail_from_level(1))
                print(f"{sample}: {'+'.join(representation['A'])}/{'+'.join(representation['B'])}")
                # TODO Select preferred alternative
                preferred = alternative
                # if preferred is None: preferred = alternative
                # if len(representation['A']) == len(representation['B']) == 1:
                #     preferred = alternative
                #     break
            print(sample, "Preferred:", preferred)
            representations[sample] = calling_to_repr(preferred, cont_graph, functions, **detail_from_level(detail_level))
        return representations


def find_core_string(match):
    """Get the core allele from a match string wise."""
    if sort_types(match) == 1:
        return match
    elif sort_types(match) == 2:
        return match[:-4]
    else:
        raise Exception(f"Unexpected match type: {match}")

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

def relevance(sample, nodes, edges, functions, supremals, reference): 
    """Determine if extra variants may be relevant for calling."""
    def overlap(a1, a2): max(a1.start, a2.start) <= min(a1.end, a2.end)
    # Get all called corealleles
    # TODO Only cores needed?
    alleles = star_allele_calling_all([sample], nodes, edges, functions, supremals, reference, detail_level=1)
    alleles = alleles[sample.split("_")[0]][sample.split("_")[1]]
    # Find extra variants: variants in the sample but not in any of the called alleles (pruned edges)
    # TODO do this on a graph?
    variants = []
    for edge in edges:
        if edge[0] == sample and sort_types(edge[1]) in (3,5):
            variants.append(edge[1])
        elif edge[1] == sample and sort_types(edge[0]) in (3,5):
            variants.append(edge[0])
    if len(variants) == 0: # No variants found
        return None
    # Check the relevance of each variant
    variants_relevance = {}
    for variant in variants:
        # Check if variant possibly interferes with any allele (overlaps with supremal)
        interferes = any([overlap(supremals[variant], supremals[allele]) for allele in alleles if allele != "CYP2D6*1"])
        # Find the impact of the variant
        impact = functions[variant]
        if sort_types(variant) == 5: # Personal variant
            severity = severity_GO(impact)
        elif sort_types(variant) == 3: # pharmvar variant
            severity = severity_pharmvar(impact)
        if severity == 2 and not interferes: # Change on protein level but doesn't interfere with any allele
            severity = 1
        variants_relevance[variant] = severity != 1 # Only not relevant if benign and doesn't interfere
        # TODO be less conservative with unknown impact?
    return variants_relevance

def detail_from_level(level):
    """Return a dictionary of options for a specific detail level
    
    The output corresponds to the options of calling_to_repr.

    Different detail levels are available.
    0: Only print best core match(es)
    1: Print all direct core matches
    2: Print all direct matches including suballeles
    3: Print all direct matches without core allele lookup
    4: Print all direct matches including the default allele
    TODO implement more levels?
    """
    kwargs = {}
    kwargs["find_cores"] = level <= 2
    kwargs["prioritize_function"] = level <= 0
    kwargs["prioritize_strength"] = level <= 0
    kwargs["suballeles"] = level >= 2
    kwargs["default"] = level == 4
    return kwargs

def calling_to_repr(calling, cont_graph, functions, find_cores, suballeles, default, prioritize_function, prioritize_strength, reorder=True):
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
    reorder: reorder alleles based on star allele number.
    # TODO add indirectly found alleles
    # TODO add alternative descriptions
    # TODO hide unparsable?
    """
    def star_num(match):
        num = match.split('*')[1]
        if num =='?':
            return 0
        return float(num)
    def remove_contained(matches, cont_graph):
        """Filter out cores that are contained in other cores"""
        # TODO do this by not adding some cores in the first place?
        matches = list(set(matches)) # Remove duplicates
        for match1 in matches[:]:  
            for match2 in matches[:]:
                if match1 == match2:
                    continue
                if sort_types(match2) == 2: # Cores can be contained in subs
                    continue
                if match2 not in cont_graph.nodes():
                    continue
                if match1 in nx.ancestors(cont_graph, match2): # match1 is contained in match2
                    matches.remove(match1)
                    break
        return matches
    representation = {}
    for phase in calling:
        representation[phase] = [] # Flatten into single list
        for alleles in calling[phase]:
            for allele in alleles:
                # Add detail
                if find_cores and sort_types(allele) == 2: # Find cores of suballeles
                    cores = find_core_traversal(allele, cont_graph) # Cores of suballele
                    representation[phase].extend(cores)
                # Remove detail
                if not suballeles and sort_types(allele) == 2: # Don't represent suballeles
                    continue
                if not default and allele == "CYP2D6*1": # Don't represent default allele if there are other matches
                    if len(representation[phase]) > 0: # Works since default is always last
                        continue
                representation[phase].append(allele) # Add match
            if prioritize_strength and len(representation[phase]) > 0: # Only keep strongest related alleles
                break
        # Remove cores contained in other cores 
        # These can occur because of get_core_traversal
        representation[phase] = remove_contained(representation[phase], cont_graph)
        # Prioritize alleles of equal rank (after filtering)
        if prioritize_function and len(representation[phase]) > 1:
            prioritized = prioritize_calling(representation[phase], functions)
            representation[phase] = list(prioritized[0]) # Only keep the first priority
        # Sort alleles based on star allele number
        representation[phase].sort(key=star_num)
    if reorder: # Sort phases based on start allele number
        sorted_alleles = sorted(representation.values(), key=lambda matches: min([star_num(match) for match in matches]))
        representation = {phase: matches for phase, matches in zip(representation.keys(), sorted_alleles)}
    return representation
        
def find_path(s, t, cont_graph, eq_graph, overlap_graph, path=None, visited=None):
    """Find a path from s to t in the graphs"""
    # TODO fix
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