import algebra as va
import networkx as nx
import warnings
from .data import api_get
import re
import copy
from itertools import combinations, combinations_with_replacement
import math
from enum import IntEnum
import functools

all_functions = ("function not assigned", 'unknown function', 'uncertain function', 'normal function', 'decreased function', 'no function')

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

class Type(IntEnum):
    """ Enum for variant types, sorted by relevance. """
    CORE = 1 # Core allele
    SUB = 2 # Suballele
    VAR = 3 # PharmVar variant
    P_VAR = 4 # Personal variant
    SAMPLE = 5 # Input sample

def find_type(v):
    """ Find type of variant. """
    # TODO use datatypes instead of strings
    if '*' in v:
        if '.' in v:
            return Type.SUB # Suballele
        else:
            return Type.CORE # Core allele
    elif v[:2] in ('HG', 'NA'):
        return Type.SAMPLE # Sample
    else:
        if ':' in v:
            return Type.VAR # Variant
        return Type.P_VAR # Personal variant
    
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
    if len(annotated["unknown function"]) + len(annotated["uncertain function"]) + len(annotated["function not assigned"]) > 0:
        # Put all in the same position of ordering
        prioritized_matching.append(
            annotated["decreased function"] |
            annotated["normal function"] |
            annotated["unknown function"] |
            annotated["uncertain function"] | 
            annotated["function not assigned"]
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
    TODO test if made redundant for calling method
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
        m_eq = set([m for m in eq_graph[sample] if find_type(m) in (Type.CORE, Type.SUB)])
        if len(m_eq) > 1: # More equivalents found, not possible for correct dataset
            raise Exception(f"{sample}: multiple equivalent matches found: {matches['equivalent']}.")
        if len(m_eq) == 1: 
            # print(sample, m_eq, "equivalent")
            matches.append(m_eq) # Strongest match
    # Find directly contained alleles
    if sample in cont_graph.nodes():
        m_cont = set([m for m, _ in cont_graph.in_edges(sample) if find_type(m) in (Type.CORE, Type.SUB)])
        if len(m_cont) > 0: 
            # print(sample, m_cont, "contained")
            matches.append(m_cont) # Less strong match than equivalent
    # Find directly overlapping alleles
    # not looking at overlaps that are a result of contained alleles
    if sample in overlap_graph.nodes():
        m_ov = set([m for m in overlap_graph[sample] if find_type(m) in (Type.CORE, Type.SUB)])
        if len(m_ov) > 0: 
            # Overlap treated with lower priority than equivalent and contained
            # print(sample, m_ov, "overlap")
            matches.append(m_ov)
    # Add default allele, will have the lowest priority
    if len(matches) == 0:
        # print(sample, "default", set([m for m, _ in cont_graph.in_edges(sample) if find_type(m) not in (Type.CORE, Type.SUB)]))
        pass
    matches.append({"CYP2D6*1",})
    # Return all matches to be filtered later (based on detail level)
    return matches

def separate_callings(unphased_calling, cont_graph, functions):
    """Unpack a calling of a unphased sample into two phased callings."""
    raise DeprecationWarning("This is now handled by alternative calling method")
    samples = list(unphased_calling.keys())
    phased_calling = {sample: {'A': [], 'B': []} for sample in sorted(samples)} 
    # Infer phasing for some samples
    for sample in samples:
        # Handle unparsable samples
        if unphased_calling[sample]['all'][-1] == {'CYP2D6*?',}:
            phased_calling[sample]['A'].append({"CYP2D6*?",})
            phased_calling[sample]['B'].append({"CYP2D6*?",})
            continue
        for i, alleles in enumerate(unphased_calling[sample]['all'][:-1]): # All non-default alleles
            phased_calling[sample]['A'].append(set()) # Maintain relation strength rank
            phased_calling[sample]['B'].append(set())
            # TODO simplify the two loops below?
            # Add largest homozygous matches to both callings
            # Don't guess the heterozygous/mixed matches but add them all to A
            # This will be handled by alternative callings later
            for allele in alleles: # is largest instead of overlap
                phased_calling[sample]['A'][-1].add(allele)
                if not any([allele in alls for alls in unphased_calling[sample]['hom']]): # check if is homozygous
                    continue
                phased_calling[sample]['B'][-1].add(allele)
            # TODO Also add homozygous matches that are not in all to both callings ?
        # Add defaults
        for phase in "AB":
            # remove empty
            phased_calling[sample][phase] = [c for c in phased_calling[sample][phase] if len(c) > 0] 
            # Add default allele as last priority
            phased_calling[sample][phase].append({"CYP2D6*1",})
    return phased_calling

temp_cache = {}
def find_ancestor_variants(allele, eq_graph, cont_graph, ov_graph):
    """Find all variants that are ancestors of allele."""
    # TODO use ov_graph to find overlapping?
    global temp_cache
    if allele in temp_cache:
        return temp_cache[allele]
    ancestors = set()
    queue = {allele,}
    while len(queue) > 0:
        a = queue.pop()
        # Check if variant 
        if find_type(a) == Type.VAR:
            ancestors.add(a)
            continue # Don't find variants in variant
        # check contained
        if a in cont_graph.nodes():
            for cont in cont_graph.predecessors(a):
                queue.add(cont)
        # check equivalents
        if a in eq_graph.nodes():
            for eq in list(nx.shortest_path(eq_graph, a))[1:]: # not self but eq
                if find_type(eq) == Type.VAR:
                    ancestors.add(eq)
                elif eq in cont_graph.nodes(): # Look at contained in eq (needed when variant exactly equals star allele)
                    for cont in cont_graph.predecessors(eq):
                        queue.add(cont)
    temp_cache[allele] = ancestors
    return ancestors

def valid_calling(sample, calling, homozygous, cont_graph, eq_graph, ov_graph, most_specific, CNV_possible=False):
    """Check if a calling is valid based on homozygous contained alleles.
    
    The homozygous alleles should be present in both phases.
    The homozygous alleles should be present only once in each phase.

    TODO don't use in generation but generate only correct solutions
    TODO test CNV_possible option
    """
    raise DeprecationWarning("This is now handled by alternative calling method which checks itself more efficiently")
    # Find ancestors for both phases
    alleles = {p: set([a for alls in calling[p] for a in alls]) for p in "AB"}
    ancestors = {p: set([v for a in alleles[p] for v in find_ancestor_variants(a, eq_graph, cont_graph, ov_graph)]) for p in "AB"}
    # Find homozygous ancestors
    hom_ancestors = set([v for h in homozygous for v in find_ancestor_variants(h, eq_graph, cont_graph, ov_graph)])
    # Homozygous variants should be present in both phases
    for ha in hom_ancestors:
        for p in "AB":
            if ha not in ancestors[p]:
                return False
    # Non-homozygous variants should be present in only one phase
    for a in ancestors['A']:
        if a in hom_ancestors:
            continue
        if a in ancestors['B']:
            return False
    # Alleles should be the most specific representation
    # TODO implement correct test for most specific
    # for p in "AB":
    #     anc = set(ancestors[p])
    #     while len(anc) > 0:
    #         largest = None, set()
    #         for specific, predecessors in most_specific.items():
    #             if predecessors <= anc: # Calling includes this allele
    #                 if len(predecessors) > len(largest[1]):
    #                     largest = specific, predecessors
    #         if largest[0] is None: # No allele found
    #             # TODO is this a valid state or not?
    #             break
    #         if largest[0] not in alleles[p]: # Most specific allele not in calling
    #             return False
    #         anc -= largest[1]
    # Only one core 
    for p in "AB":
        cores = set([a for a in alleles[p] if find_type(a) == Type.CORE])
        cores = [c for c in cores if c != "CYP2D6*1"]
        if len(cores) > 1:
            return False
    # Find all variants that should be present in a calling
    all_ancestors = find_ancestor_variants(sample + "_all", eq_graph, cont_graph, ov_graph)
    if sample + "_all" in cont_graph.nodes():
        for ancestor in cont_graph.predecessors(sample + "_all"):
            if find_type(ancestor) == Type.VAR:
                all_ancestors.remove(ancestor)
    # Calling should contain all variants
    if all_ancestors != ancestors["A"] | ancestors["B"]:
        return False
    return True

def generate_alternative_callings_top_down(sample, start_calling, homozygous, cont_graph, eq_graph, ov_graph, functions, detail_level):
    """Generate unique valid alternative callings for a sample.
    
    Function attempts to only generate valid callings.
    Valid callings contain all variants contained in the same in the proper homozygous state.
    Furthermore they are the most specific description: *39+*34 = *2

    Finds new callings by extending existing callings.
    Extending is replacing an allele with its underlying alleles or variants.
    Extending is possible with all underlying or with itself and a homozygous allele (in the different phase).
    At each level of extending the function tries all possible placements.

    This method doesn't look ahead to the variants.
    Another approach is the bottom up approach which starts from the variants.

    TODO fix not working in all cases
    TODO only return valid results (without looking ahead)
    TODO only return unique results efficiently
    """
    raise DeprecationWarning("Doesn't work in all cases")
    def placements(calling):
        """Generate all possible placements for a calling."""
        # TODO test if valid
        # TODO don't place if not valid without looking ahead
        alleles = set([a for alleles in calling['A'][:-1] for a in alleles])
        for r in range(1, math.floor(len(alleles) / 2) + 1): # All unique combinations (2C5 = 3C5)
            for movement in combinations(alleles, r):
                movement = set(movement)
                # Construct calling
                new_calling = {
                    "A": [set([a for a in alls if a not in movement]) for alls in calling['A']],
                    "B": [set([a for a in alls if a in movement]) for alls in calling['A']],
                }
                for i, alleles in enumerate(reversed(calling['B'])):
                    new_calling['B'][-(i+1)] |= alleles
                yield new_calling
                if r == 1 and len(alleles) == 2: break # only one calling possible
        # Add homozygous to other phase
        # TODO valid?
        new_calling = {
            "A": [set([a for a in alls]) for alls in calling['A']],
            "B": [set(homozygous), {"CYP2D6*1",}],
        }
        yield new_calling
        # Don't move, all on one phase
        # last in order since this notation is not preferred
        yield calling 

    # Find star allele definitions to check for most specific
    # TODO handle suballele detail level mismatch
    definitions = {}
    for node in cont_graph.nodes():
        # if find_type(node) != Type.CORE or detail_level > 1 and find_type(node) != Type.SUB:
        if find_type(node) not in (Type.CORE, Type.SUB):
            continue
        ancestors = find_ancestor_variants(node, eq_graph, cont_graph, ov_graph)
        if not ancestors:
            continue
        definitions[node] = ancestors
    # BFS over valid alternative callings
    queue = [start_calling]
    # Track unique representations 
    unique = set()
    while len(queue) > 0:
        calling = queue.pop()
        # Try different placements
        for placement in placements(calling):
            # print(placement)
            if not valid_calling(sample, placement, homozygous, cont_graph, eq_graph, ov_graph, definitions):
                continue
            representation = calling_to_repr(placement, cont_graph, functions, **detail_from_level(detail_level), reorder=True)
            representation = f"{'+'.join(representation['A'])}/{'+'.join(representation['B'])}"
            if representation in unique:
                continue
            unique.add(representation)
            print(representation)
            # TODO don't have to extend here? already valid with less specific
            yield placement
        # Extend the alleles at this level
        for i, alleles in enumerate(calling['A'][:-1]): # all non default matches
            for allele in alleles:
                # Find underlying alleles and variants
                # TODO overlap?
                underlying = set([u for u in cont_graph.predecessors(allele) if find_type(u) in (Type.CORE, Type.SUB, Type.VAR)])
                # Extend by replacing with underlying
                # TODO avoid duplicates (A -> a; B -> b and B -> b; A -> a in next iteration)
                if not underlying: continue # Don't replace equivalent alleles
                new_calling = {
                    "A": [set([a for a in alls if a != allele]) for alls in calling['A']],
                    "B": [{"CYP2D6*1",}]
                }
                new_calling['A'][i] |= underlying
                queue.append(new_calling)

def generate_alternative_callings_bottom_up(sample, homozygous, cont_graph, eq_graph, ov_graph, functions, detail_level):
    """Generate all valid alternative calling in a bottom up approach.
    
    Should yield the same results as generate_alternative_callings.
    But is more efficient as it considers only unique valid callings.

    Only returns valid answers as judged by valid_calling.

    TODO fix most specific finding
    TODO consider overlap
    TODO change how alternatives are generated based on detail level
        now returning wrong results for a detail level 1
    TODO prevent duplicates being generated
    TODO make more efficient by reverse hash of definitions
    """
    raise DeprecationWarning("Finding of most specific allele is not correct")
    # Find star allele definitions to check for most specific
    definitions = {}
    alleles = {}
    for node in cont_graph.nodes():
        # Only look at core allele definitions if suballeles are not considered
        if find_type(node) != Type.CORE and (detail_level <= 1 or find_type(node) != Type.SUB):
            continue
        if find_type(node) == Type.SUB:
            core = find_core_string(node)
            if core not in alleles:
                alleles[core] = set()
            alleles[core].add(node)
        ancestors = find_ancestor_variants(node, eq_graph, cont_graph, ov_graph)
        if not ancestors:
            pass # TODO handle (overlap)
        definitions[node] = ancestors
    # Find what variants the sample contains indirectly
    all_variants = find_ancestor_variants(sample + "_all", eq_graph, cont_graph, ov_graph)
    if sample + "_all" in cont_graph.nodes():
        for variant in cont_graph.predecessors(sample + "_all"):
            if find_type(variant) != Type.VAR:
                continue
            all_variants.remove(variant)
    # Split into homozygous and heterozygous
    hom_variants = set()
    for h in homozygous:
        for a in find_ancestor_variants(h, eq_graph, cont_graph, ov_graph):
            hom_variants.add(a)
    het_variants = all_variants - hom_variants
    # Find all possible combinations of heterozygous variants
    unique = set()
    for r in range(0, math.floor(len(het_variants) / 2) + 1): # All unique combinations (5C2 = 5c3)
        for combination in combinations(het_variants, r):
            # Combination forms two unordered groups of variants
            combination = set(combination)
            variants = {
                "A": combination | hom_variants,
                "B": all_variants - combination | hom_variants,
            }
            # print(variants)
            # Convert to calling
            calling = {}
            for phase in ("A", "B"):
                calling[phase] = []
                core_covers = set()
                core_alleles = set()
                # Find most specific allele descriptions for variants
                # TODO correct heuristic for most specific?
                # TODO use td structure for this
                for core, subs in alleles.items():
                    covers = set()
                    most_specific = set()
                    for allele in subs | {core,}:
                        if core == "CYP2D6*1":
                            continue
                        ancestors = definitions[allele]
                        if ancestors <= variants[phase]: # definition is subset of variants
                            if len(covers | ancestors) >= len(covers): # Covers more of the variants
                                if covers | ancestors == ancestors: # This allele is more specific than the current coverage
                                    most_specific = {allele,}
                                    covers = set(ancestors)
                                elif not (ancestors <= covers): # This allele should be present together with the current one
                                    most_specific.add(allele)
                                    covers |= ancestors                
                    # Add variants not covered by any allele definition
                    most_specific |= (variants[phase] - covers)
                    if len(covers) > len(core_covers):
                        core_covers = covers
                        core_alleles = most_specific
                calling[phase].append(core_alleles)
                calling[phase].append({"CYP2D6*1",})
            # Return valid calling
            # if not valid_calling(sample, calling, homozygous, cont_graph, eq_graph, ov_graph, definitions):
            #     raise Exception("Invalid calling generated")
            representation = calling_to_repr(calling, cont_graph, functions, **detail_from_level(detail_level), reorder=True)
            representation = f"{'+'.join(representation['A'])}/{'+'.join(representation['B'])}"
            # Only return unique representations (duplicates can be generated because of detail level)
            if representation in unique:
                continue
            print(calling)
            print(representation)
            unique.add(representation)
            # print(calling)
            # print(representation)
            # print()
            yield calling
            # TODO valid?
            if r == 1 and len(het_variants) == 2: # only one combination possible
                break

def generate_alternative_callings(sample, homozygous_alleles, hom_variants, cont_graph, eq_graph, ov_graph, functions, filter_default=False):
    """ 
    Checks which alleles the sample contains and sees if a valid calling (*x/*y) n_cores = 2 is possible.
    If so, test if the distribution of variants is valid by homozygosity.
    If not remove detail by extending the alleles with their ancestors until a valid calling is found.

    TODO stopping condition
    TODO fix runtime with suballeles
    TODO ordering
    TODO handle overlap?
    TODO do not create new state but mutate (optimisation)
    """
    def valid(_ancestors, hom_variants, heterozygous):
        """ Check if the distribution of variants is valid by homozygosity """
        # TODO do not allow absent homs (with new definition of homs not needed)
        # Hom must be present in both phases
        for hom in hom_variants:
            c = 0
            if hom in _ancestors[0]: c += 1
            if hom in _ancestors[1]: c += 1
            if c == 1: # Cannot be present in single phase but can be absent due to extending
                return False
        # Het must be present in only one phase
        for het in heterozygous:
            c = 0
            if het in _ancestors[0]: c += 1
            if het in _ancestors[1]: c += 1
            if c > 1: # Cannot be present in both phases but can be absent due to extending
                return False
        return True
    def generate_callings(state):
        """ Generate all possible callings for a given state """
        cores = list(set((find_core_string(a) for a in state if 
                            find_type(a) != Type.VAR and 
                            find_type(a) != Type.P_VAR and 
                            find_core_string(a) != "CYP2D6*1")))
        # Can only result in valid state when it consists of 2 cores or less (multiplied by count)
        # Here, the suballeles of the same core are grouped and suballeles of 1 are allowed
        if len(cores) <= 2 and not \
            (len(cores) == 2 and 
                any((state.count(a) > 1 for a in state if find_core_string(a) in cores))): 
            # Find base calling, alleles of different cores must be in different phases
            _calling_base = [set(), set()]
            free = set()
            for c in state:
                if find_core_string(c) == "CYP2D6*1": # Are free unless homozygous
                    if state.count(c) == 1:
                        free.add(c)
                    else:
                        _calling_base[1].add(c)
                        _calling_base[0].add(c)
                    continue
                # Distribute alleles over phases (group cores together)
                if c in _calling_base[0]: # Twice in state
                    _calling_base[1].add(c)
                i = cores.index(find_core_string(c))
                _calling_base[i].add(c)
            # Check possible distributions of free moving alleles
            for r in range(0, len(free)+1):
                for f in combinations(free, r):
                    f = set(f)
                    _calling = [_calling_base[0] | f, _calling_base[1] | free - f]
                    _pattern = [set(), set()]
                    for i in range(2):
                        for allele in _calling[i]: 
                            _pattern[i] |= allele_definitions[allele]
                    yield _calling, _pattern
    # Find star allele definitions TODO move
    allele_definitions = {}
    suballeles = {}
    for node in cont_graph.nodes():
        if find_type(node) != Type.CORE and find_type(node) != Type.SUB:
            continue
        if find_type(node) == Type.SUB:
            core = find_core_string(node)
            if core not in suballeles: suballeles[core] = set()
            suballeles[core].add(node)
        ancestors = find_ancestor_variants(node, eq_graph, cont_graph, ov_graph)
        allele_definitions[node] = ancestors
    # Find contained alleles of sample
    if sample + "_all" in cont_graph.nodes():
        alleles = list((a for a in cont_graph.predecessors(sample + "_all") if find_type(a) != Type.VAR and find_type(a) != Type.P_VAR))
    elif sample + "_all" in eq_graph.nodes(): 
        # Can occur for simple homozygous callings
        alleles = list((a for a in eq_graph[sample + "_all"] if find_type(a) != Type.VAR and find_type(a) != Type.P_VAR))
    else:
        # No alleles found, empty calling 
        # only alternative will be {}/{} which is returned as *1/*1
        alleles = list()
    # Contained alleles consisting of homozygous variants
    homozygous_alleles = list(homozygous_alleles) 
    # [Optional] Ignore suballeles of default allele as these will be filtered out later (optimisation)
    if filter_default:
        alleles = [a for a in alleles if find_core_string(a) != "CYP2D6*1"]
        homozygous_alleles = [a for a in homozygous_alleles if find_core_string(a) != "CYP2D6*1"]
    # Find fundamental variants of sample
    variants = find_ancestor_variants(sample + "_all", eq_graph, cont_graph, ov_graph)
    if sample + "_all" in cont_graph.nodes(): 
        # Filter out 'extra variants'
        variants -= set(cont_graph.predecessors(sample + "_all"))
    # Check which variants are heterozygous (homozygous known from phasing)
    hom_variants = set(hom_variants)
    het_variants = variants - hom_variants 
    # Only use homozygous variants that form a complete allele
    # QUESTION is this valid?
    hom_variants = set((a for a in hom_variants if any((a in allele_definitions[h] for h in homozygous_alleles))))
    # Add initial state
    queue = []
    queue.append((list(alleles), 0, True, 0))
    # Add some initial alleles twice
    # when these contain a homozygous variant that is not present in another allele 
    # TODO limit? (optimisation)
    # as these may be needed to arrive at a valid state
    # TODO can merge with above?
    # TODO simplify
    # for a in alleles:
    #     hom_anc = allele_definitions[a] & hom_variants
        # for o in alleles:
        #     if o == a:
        #         continue
        #     if find_core_string(o) == find_core_string(a):
        #         continue
        #     hom_anc -= allele_definitions[o]
        # if len(hom_anc) > 0 and a not in homozygous_alleles: # Contains hom but isn't hom itself
        #     if len(queue) == 1:
        #         new_state = list(alleles)
        #         new_state.insert(new_state.index(a), a)
        #         queue.insert(0, (new_state, 0, False)) # Don't call on first (not a valid state)
        #     else:
        #         queue[0][0].insert(queue[1][0].index(a), a)
    # If homozygous alleles are already present a valid state must include these
    for a in alleles:
        if a in homozygous_alleles:
            for q in queue:
                q[0].insert(q[0].index(a), a)
    count = 0
    # print(*queue, sep="\n")
    while len(queue) > 0:
        state, extended, call, parent = queue.pop(0)
        any_valid = False
        if call:
            count += 1
            # print(count, state)
            # Only try generating a calling of a valid number of cores
            # (65.1,2.2,10.1 will never form a valid calling of two real alleles)
            for _calling, _pattern in generate_callings(state):
                if valid(_pattern, hom_variants, het_variants):
                    calling = {"A": [_calling[0], {"CYP2D6*1",}], "B": [_calling[1], {"CYP2D6*1",}]}
                    yield calling
                    any_valid = True
        # Stop as can extend no further
        if extended >= len(state):
            continue
        # Stop as all further will be less specific
        # TODO fix
        # TODO refine ...
        if any_valid:
            continue
        # 1) Do not extend this one but continue
        # Ensures that all possible callings are generated
        queue.append((list(state), extended + 1, False, count))
        # Extend alleles with underlying alleles
        # Find underlying alleles not contained in other alleles
        # Compresses combinations that may arise due to duplicate containment (both 65 and 2.2 contain 2)
        extend = state[extended]
        removed = set()
        underlying = list()
        _underlying = set()
        # TODO use ancestors function limited to depth 1?
        if extend in cont_graph.nodes():
            for c in cont_graph.predecessors(extend): 
                _underlying.add(c)
        if extend in eq_graph.nodes():
            for c in eq_graph[extend]: 
                if find_type(c) != Type.VAR:
                    continue # TODO handle eq suballeles
                _underlying.add(c)
        for u in _underlying:
            # Lose detail by ignoring variants
            if find_type(u) == Type.VAR or find_type(u) == Type.P_VAR:
                removed.add(u)
                continue
            # Check if u is contained in or equal to other alleles (ignore extended)
            # Assume that the inverse does not happen TODO check
            if any((u in nx.ancestors(cont_graph, a) for a in state[extended+1:])):
                continue
            if u in state[extended+1:]:
                continue
            if filter_default and find_core_string(u) == "CYP2D6*1":
                continue
            underlying.append(u)
            # Present in state twice
            # if u in homozygous_alleles: 
            #     underlying.append(u)
        # TODO Don't extend if this does not remove any information (optimisation)
        # Prevents different paths finding the same conclusion (e.g. 2.2,10.4>2,10 and 65>2,10)
        # if len(underlying) == len(removed) == 0:
        #     continue
        # Do not extend homozygous alleles as this will never result in a more specific valid calling 
        # QUESTION is this valid / needed
        # if extend in homozygous_alleles:
        #     continue
        # 2) Replace allele with underlying alleles
        new_state = [state[i] for i in range(len(state)) if i != extended]
        for u in underlying:
            new_state.insert(extended, u)
        # print(" ", state)
        # print(" ", "extend", extend, "at", extended, 'with', underlying)
        # print(" ", new_state)
        queue.append((new_state, extended + 1 - 1, True, count))
    # print(count)

def order_callings(calling, functions, no_default=True, shortest=True, no_uncertain=True):
    """Order alternative callings by clinical relevance.
    
    no_default: prefer callings with not all on one side, (e.g. *X/*Y over *1+*X/*Y)
    shortest: Prefer shorter callings, equal to without '+' (e.g. *Q/*Z over *X+*Y/*Z)
    no_uncertain: Prefer certain function over uncertain ones.
    
    TODO instead use a generation
    TODO prefer worse function over better one?
    TODO make this a comparison function?
    """
    raise DeprecationWarning("Order will be fixed by generation")
    score = 0
    any_default = False
    for phase in ("A", "B"):
        matches = set()
        uncertain_function = False
        for alleles in calling[phase][:-1]:
            for allele in alleles:
                if find_type(allele) not in (Type.CORE, Type.SUB):
                    continue
                matches.add(allele)
                function = functions[allele]
                if all_functions.index(function) <= 2:
                    uncertain_function = True
        if no_uncertain and uncertain_function:
            score += 1
        if shortest:
            score += len(matches) * 10 # Assume no more than 99 matches
        if len(matches) == 0:
            any_default = True
    if no_default and any_default and score >= 20: 
        score += 1000 # Prefer callings with not all on one side
    return score
                    

def star_allele_calling_all(samples, nodes, edges, functions, supremals, reference, homozygous=None, phased=True, detail_level=1, reorder=True):
    """Iterate over samples and call star alleles for each."""
    eq_graph = nx.Graph([(left, right) for left, right, relation in edges if relation == va.Relation.EQUIVALENT])
    cont_graph = nx.DiGraph([(left, right) for left, right, relation in edges if relation == va.Relation.IS_CONTAINED])
    ov_graph = nx.Graph([(left, right) for left, right, relation in edges if relation == va.Relation.OVERLAP])

    # Call each sample
    callings = {sample.split('_')[0]: {} for sample in sorted(samples)} 
    for sample in samples:
        calling = star_allele_calling(sample, eq_graph, cont_graph, ov_graph, functions, supremals, reference)
        sample_source, phasing = sample.split('_')
        callings[sample_source][phasing] = calling
    # Create a textual representation of the calling based on the amount of detail needed
    if phased: # Calling is phased
        representations = {sample: calling_to_repr(callings[sample], cont_graph, functions, **detail_from_level(detail_level), reorder=reorder) for sample in callings}
        return representations
    if not phased: # Unphased calling should be separated
        representations = {}
        # TODO make this into function
        # TODO allow for keeping of multiple alternative representations?
        # test_i = 0
        for sample, calling in callings.items():
            # DEBUG
            # if sample != "NA10859": continue # Small tree
            # if sample != "HG00421": continue # Common basic difficult pattern
            # if sample != "HG00337": continue # Simple straightforward solution
            # if sample != "HG00423": continue # nearly fully homozygous
            # if sample != "NA19143": continue # Most complex bu
            # if sample != "HG00589": continue # Need for allowing individual variants
            # if sample != "NA12815": continue # Calling contains allele with contained allele
            # if sample != "HG00111": continue # Simple homozygous (eq)
            # if sample != "NA18861": continue # Homozygous
            # if sample != "HG00276": continue # Fully homozygous with multiple direct subs
            # if sample != "HG00373": continue # Unparsable
            # if sample != "NA19143": continue # Example of multiple homozygous needed in start state
            # if sample != "NA18518": continue # Example of restricting the hom variants to entire alleles
            # if sample != "NA12006": continue # Previously wrong answer 
            # if sample != "HG01190": continue # Most complex td TODO fix runtime?
            # if sample != "NA19122": continue # Example of multiple homozygous contained

            # test_i += 1
            # if test_i >= 5:
            #     exit()
            # Filter out unparsable
            if calling["all"][-1] == {"CYP2D6*?",}:
                representations[sample] = {"A": ["CYP2D6*?",], "B": ["CYP2D6*?",]}
                print(sample)
                print(f"{','.join(representations[sample]['A'])}/{','.join(representations[sample]['B'])}")
                print()
                continue
            # All homozygous alleles for the current sample
            homozygous_alleles = set([allele for alleles in calling['hom'] for allele in alleles if allele != "CYP2D6*1"])
            # Generate unique valid alternative callings
            alternatives = generate_alternative_callings(sample, homozygous_alleles, homozygous[sample], cont_graph, eq_graph, ov_graph, functions, filter_default=True)
            print(sample)
            preferred = None
            unique_repr = set()
            for alternative in alternatives:
                representation = calling_to_repr(alternative, cont_graph, functions, **detail_from_level(detail_level), reorder=True)
                representation = f"{','.join(representation['A'])}/{','.join(representation['B'])}"
                # filter out non-unique representations (arise due to detail level mismatch)   
                if representation in unique_repr:
                    continue
                unique_repr.add(representation) 
                if preferred is None:
                    preferred = alternative
                print(representation)
            # Select the most relevant alternative
            if not preferred:
                # TODO handle differently?
                raise Exception("No valid alternative callings found for ", sample)
            representations[sample] = calling_to_repr(preferred, cont_graph, functions, **detail_from_level(detail_level), reorder=reorder)
            print()
        return representations


def find_core_string(match):
    """Get the core allele from a match string wise."""
    if find_type(match) == Type.CORE:
        return match
    elif find_type(match) == Type.SUB:
        return match[:-4]
    else:
        raise Exception(f"Unexpected match type: {match}")

def find_core_traversal(match, cont_graph):
    """Get the coreallele(s) from a match based on traversal.
    
    Returns a set since a suballele can contain more corealleles.
    Will not return indirectly contained corealleles.
    This won't find the corealleles that are not contained in the suballele by design.
    """
    raise DeprecationWarning("Not used as find core string is used instead")
    cores = []
    if match in cont_graph.nodes(): # Only have to check contained since everything is connected to the core
        for node, _ in cont_graph.in_edges(match):
            if find_type(node) == Type.CORE: # Found core contained in match, don't traverse further
                cores.append(node)
            elif find_type(node) == Type.SUB: # Traverse further
                cores.extend(find_core_traversal(node, cont_graph))
    return cores

def detail_from_level(level):
    """Return a dictionary of options for a specific detail level
    
    The output corresponds to the options of calling_to_repr.

    Different detail levels are available.
    0: Only print best core match(es)
    1: Print all direct core matches
    2: Print all direct sub allele matches
    3: Print all direct matches without core allele lookup
    4: Print all direct matches including the default allele
    TODO implement more levels?
    """
    if level == 0:
        warnings.warn("Detail level 0 may lose some useful information")
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
                if find_cores and find_type(match2) == Type.SUB: # Cores can be contained in subs
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
                if find_cores and find_type(allele) == Type.SUB: # Find cores of suballeles
                    core = find_core_string(allele) # Cores of suballele
                    if core != "CYP2D6*1": # Wild type will already be present
                        representation[phase].append(core)
                # Remove detail
                if not suballeles and find_type(allele) == Type.SUB: # Don't represent suballeles
                    continue
                if find_type(allele) in (Type.VAR, Type.P_VAR): # No variants in calling TODO make parameter
                    continue
                if not default and allele == "CYP2D6*1": # Don't represent default allele if there are other matches
                    if len(representation[phase]) > 0: # Works since default is always last
                        continue
                representation[phase].append(allele) # Add match
            if prioritize_strength and len(representation[phase]) > 0: # Only keep strongest related alleles
                break
        # Remove suballeles of default allele if there are other core matches
        if not default and len(set([find_core_string(a) for a in representation[phase]])) > 1:
            representation[phase] = [a for a in representation[phase] if find_core_string(a) != "CYP2D6*1"]
        # Remove alleles contained in others 
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
