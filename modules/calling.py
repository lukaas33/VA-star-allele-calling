import algebra as va
import networkx as nx
import warnings
from .data import api_get
from .other_sources import find_id_hgvs, get_annotation_entrez, severity_GO, severity_pharmvar
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
            matches.append(m_eq) # Strongest match
    # Find directly contained alleles
    if sample in cont_graph.nodes():
        m_cont = set([m for m, _ in cont_graph.in_edges(sample) if find_type(m) in (Type.CORE, Type.SUB)])
        if len(m_cont) > 0: 
            matches.append(m_cont) # Less strong match than equivalent
    # Find directly overlapping alleles
    # not looking at overlaps that are a result of contained alleles
    if sample in overlap_graph.nodes():
        m_ov = set([m for m in overlap_graph[sample] if find_type(m) in (Type.CORE, Type.SUB)])
        if len(m_ov) > 0: 
            # Overlap treated with lower priority than equivalent and contained
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

def find_ancestor_variants(allele, eq_graph, cont_graph, ov_graph):
    """Find all variants that are ancestors of allele."""
    # TODO use ov_graph?
    # TODO return smallest variant?
    queue = {allele,}
    while len(queue) > 0:
        a = queue.pop()
        # Check if variant 
        if find_type(a) == Type.VAR:
            yield a
            continue # Don't find variants in variant
        # check contained
        if a in cont_graph.nodes():
            for cont in cont_graph.predecessors(a):
                queue.add(cont)
        # check equivalents
        if a in eq_graph.nodes():
            for eq in list(nx.shortest_path(eq_graph, a))[1:]: # not self but eq
                if find_type(eq) == Type.VAR:
                    yield eq
                elif eq in cont_graph.nodes(): # Look at contained in eq (needed when variant exactly equals star allele)
                    for cont in cont_graph.predecessors(eq):
                        queue.add(cont)

def valid_calling(sample, calling, homozygous, cont_graph, eq_graph, ov_graph, most_specific, CNV_possible=False):
    """Check if a calling is valid based on homozygous contained alleles.
    
    The homozygous alleles should be present in both phases.
    The homozygous alleles should be present only once in each phase.

    TODO don't use in generation but generate only correct solutions
    TODO test CNV_possible option
    """
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
    for p in "AB":
        anc = set(ancestors[p])
        while len(anc) > 0:
            largest = None, set()
            for specific, predecessors in most_specific.items():
                if predecessors <= anc: # Calling includes this allele
                    if len(predecessors) > len(largest[1]):
                        largest = specific, predecessors
            if largest[0] is None: # No allele found
                # TODO is this a valid state or not?
                break
            if largest[0] not in alleles[p]: # Most specific allele not in calling
                return False
            anc -= largest[1]
    
    # Find all variants that should be present in a calling
    all_ancestors = set(find_ancestor_variants(sample + "_all", eq_graph, cont_graph, ov_graph))
    if sample + "_all" in cont_graph.nodes():
        for ancestor in cont_graph.predecessors(sample + "_all"):
            if find_type(ancestor) == Type.VAR:
                all_ancestors.remove(ancestor)
    # Calling should contain all variants
    if all_ancestors != ancestors["A"] | ancestors["B"]:
        return False
    return True

def old_generate_alternative_callings(sample, calling, homozygous, cont_graph, eq_graph, ov_graph, functions, detail_level):
    """Generate valid alternative callings for a sample in a BFS manner.
    
    Generates most specific first.
    Option to stop after first valid calling is found (find_all).

    Only returns valid callings based on homozygous matches.

    Only returns calling with a unique representation according to the required detail level.

    Need to (deep)copy if generator not immediately consumed.

    # TODO improve runtime for suballeles
    # TODO prevent duplicates in generation instead of filtering (at move for instance)
    # TODO prevent invalid generation (move to invalid state or extend with invalid state)
    # TODO don't deepcopy everything?
    # TODO fix runtime NA19174
    # TODO validate on suballeles
    """
    raise DeprecationWarning("Use generate_alternative_callings instead")
    def underlying(allele, cont_graph):
        """Return the underlying alleles of some allele."""	
        # TODO look at overlapping
        for contained, _ in cont_graph.in_edges(allele): # Find underlying alleles
            if find_type(contained) in (Type.CORE, Type.SUB):
                yield contained
    def move_alleles(calling):
        """Move alleles to different phase"""
        # TODO don't move to invalid state (keep homozygous in both)
        # TODO don't move to have alleles within a phase be contained in each other
        # Move some alleles
        for i, alleles in enumerate(calling['A'][:-1]): # All non-default matches for this sample; maintain rank
            alleles = tuple(alleles)
            # TODO Prevent duplicates here (X/Y = Y/X)
            for r in range(1, len(alleles)): # Don't move all (r < n) since this would result in duplicates
                for move in combinations(alleles, r):
                    moved_alt = copy.deepcopy(calling)
                    # Make B the same length as A
                    while len(moved_alt['B']) < len(moved_alt['A']): 
                        moved_alt['B'].insert(0, set())
                    # Change phase 
                    moved_alt['A'][i] -= set(move)
                    moved_alt['B'][i] |= set(move)
                    yield moved_alt
        # Don't move anything (*1/*X+*Y)
        # Last yielded since this notation is not preferred
        yield copy.deepcopy(calling)
    # Find star allele definitions consisting of only star alleles
    extendable = {}
    for node in cont_graph.nodes():
        if find_type(node) not in (Type.CORE, Type.SUB):
            continue
        predecessors = set(cont_graph.predecessors(node))
        if not predecessors:
            continue
        if all([find_type(s) in (Type.CORE, Type.SUB) for s in predecessors]):
            extendable[node] = predecessors
    # print(valid_calling(sample, {"A": [{"CYP2D6*4", "CYP2D6*10", "CYP2D6*74"}], "B": [{"CYP2D6*41"}]}, homozygous, cont_graph, eq_graph, ov_graph, extendable))
    # Perform BFS on callings
    queue = [(calling, set(), 0)]
    unique = set()
    while len(queue) > 0:
        any_valid = False
        # First alternative calling in the queue
        alternative, extended, depth = queue.pop(0) 
        # Try different phasing at the current depth
        for moved_alt in move_alleles(alternative):
            # Check if the representation of this alternative is unique
            if valid_calling(sample, moved_alt, homozygous, cont_graph, eq_graph, ov_graph, extendable):
                representation = calling_to_repr(moved_alt, cont_graph, functions, **detail_from_level(detail_level), reorder=True)
                representation = f"{'+'.join(representation['A'])}/{'+'.join(representation['B'])}"
                if representation in unique: continue # Only return unique
                print(representation)
                unique.add(representation)
                any_valid = True
                yield moved_alt
        # Don't look deeper when calling valid (don't replace description with same information)
        # if any_valid: # TODO test
        #     continue
        # Look deeper in the tree
        for i, alleles in enumerate(alternative['A'][:-1]): # All non-default alleles
            alleles = set(alleles)
            for r in range(1, len(alleles)+1):
                # Prevent extending a previous level (will be allowed if a allele occurs again in a different level)
                for extend in combinations(alleles - extended, r=r): 
                    # print("extend", extend)
                    extend_alt = copy.deepcopy(alternative)
                    extend = set(extend)
                    # TODO can merge cases?
                    # 1. Replace allele with underlying alleles
                    # TODO is this approach valid?
                    if not any_valid: # Don't replace if already valid as this adds no information
                        for allele in extend:                        
                            underlying_alleles = set(underlying(allele, cont_graph))
                            extend_alt['A'][i] -= {allele,}
                            extend_alt['A'][i] |= underlying_alleles
                        queue.append((extend_alt, extended | alleles - underlying_alleles, depth+1)) 
                    # 2. Add underlying alleles to extended together with parent allele if homozygous
                    extend_alt = copy.deepcopy(alternative)
                    for allele in extend:
                        # TODO allow variable number of underlying alleles?
                        # TODO don't extend if number of a homozygous allele is > 2
                        underlying_alleles = set(underlying(allele, cont_graph))
                        for underlying_allele in underlying_alleles:
                            # Extend if homozygous or contains homozygous  
                            if not any(h in nx.ancestors(cont_graph, underlying_allele) | {underlying_allele,} for h in homozygous):
                                continue
                            # Add homozygous underlying and don't remove parent
                            # Must be in the other phase as the parent for this to be valid (will be handled by move_alleles)
                            underlying_alleles.add(underlying_allele)
                            extend_alt['A'][i] |= {underlying_allele,}
                    queue.append((extend_alt, extended | alleles - underlying_alleles, depth+1))

def generate_alternative_callings_recursive(calling, cont_graph, homozygous, extended, find_all=False, depth=1):
    """Generate valid alternative callings for a sample.
    
    This is useful for unphased not homozygous matches.
    Need to copy yielded callings to prevent changing the original calling
    would not be needed in case of one iteration over this generator

    Find_all can be used to find all alternative callings, if false valid callings are extended into less specific callings.

    # TODO breadthfirst instead of depth first (lower are more specific)
    # TODO don't extend if valid solution with more specific alleles is found (improve this to work with HG00421)
    # TODO fix runtime of HG01190 and of with suballeles
    """
    raise DeprecationWarning("This function is deprecated, use generate_alternative_callings instead")	
    def move_alleles(calling):
        """Move alleles to other phase"""
        # Move
        for i, alleles in enumerate(calling['A'][:-1]): # All non-default matches for this sample; maintain rank
            # Don't move all since this would result in a duplicate (X/Y = Y/X)
            for r in range(1, len(alleles)):
                # TODO how to avoid duplicates here?
                for move in combinations(alleles, r=r):
                    moved_calling = copy.deepcopy(calling)
                    # Make B the same length as A
                    while len(moved_calling['B']) < len(moved_calling['A']): 
                        moved_calling['B'].insert(0, set())
                    # Change phase 
                    moved_calling['A'][i] -= set(move)
                    moved_calling['B'][i] |= set(move)
                    yield moved_calling
        # Don't move
        # Last since this notation is not preferred
        yield copy.deepcopy(calling)
    # Alternatives from moving at this depth
    any_valid = False
    for alternative in move_alleles(copy.deepcopy(calling)):
        # Return only valid callings
        if valid_calling(alternative, cont_graph, homozygous):
            any_valid = True
            yield alternative
    # Go deeper by extending one allele
    # Don't extend if already valid TODO test
    if not any_valid or find_all:
        # Find allele to extend
        extend_allele, extend_underlying = None, None
        for i, alleles in enumerate(calling['A'][:-1]): 
            for allele in alleles:
                # Avoid infinite recursion
                # TODO test
                if allele in extended:
                    continue
                # Find underlying alleles of this allele
                if allele not in cont_graph.nodes(): 
                    continue
                # TODO check overlap graph
                underlying = set()
                for contained, _ in cont_graph.in_edges(allele): # Find underlying alleles
                    if find_type(contained) not in (Type.CORE, Type.SUB):
                        continue
                    underlying.add(contained)
                if len(underlying) == 0:
                    continue
                extend_allele, extend_underlying = allele, underlying
                break # Don't extend further at this depth
            else: # No break
                continue
            break
        if extend_allele is not None:
            # 1. Go deeper without extending or replacing the current
            for alternative in generate_alternative_callings(calling, cont_graph, homozygous, extended | {extend_allele,}, depth+1):
                yield alternative
            # 2. Go deeper with allele replaced by underlying
            replaced_calling = copy.deepcopy(calling)
            replaced_calling['A'][i].remove(extend_allele)
            replaced_calling['A'][i] |= extend_underlying
            # print(depth, "replace", allele, "with", underlying, "-", replaced_calling)
            for alternative in generate_alternative_callings(replaced_calling, cont_graph, homozygous, extended | {extend_allele,}, depth+1):
                yield alternative
            # 3. Go deeper with extended alleles
            # Don't remove but add underlying which may be contained in it (if they are homozygous)
            for r in range(1, len(extend_underlying)): 
                for extend in combinations(extend_underlying, r):
                    if not all([a in homozygous for a in extend]): # Can only extend an allele with itself and a homozygous allele
                        continue
                    extended_calling = copy.deepcopy(calling)
                    extended_calling['A'][i] |= set(extend)
                    for alternative in generate_alternative_callings(extended_calling, cont_graph, homozygous, extended | {extend_allele,}, depth+1):
                        yield alternative

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
        ancestors = set(find_ancestor_variants(node, eq_graph, cont_graph, ov_graph))
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

    TODO consider overlap
    TODO change how alternatives are generated based on detail level
        now returning wrong results for a detail level 1
    TODO check for duplicates being generated
    TODO make more efficient by reverse hash of definitions
    """
    # Find star allele definitions to check for most specific
    definitions = {}
    for node in cont_graph.nodes():
        # if find_type(node) != Type.CORE or detail_level > 1 and find_type(node) != Type.SUB:
        if find_type(node) not in (Type.CORE, Type.SUB):
            continue
        ancestors = set(find_ancestor_variants(node, eq_graph, cont_graph, ov_graph))
        if not ancestors:
            continue
        definitions[node] = ancestors
    # Find what variants the sample contains indirectly
    all_variants = set(find_ancestor_variants(sample + "_all", eq_graph, cont_graph, ov_graph))
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
            # Convert to calling
            calling = {}
            for phase in ("A", "B"):
                calling[phase] = []
                while len(variants[phase]) > 0:
                    # Find most specific allele descriptions for variants
                    most_specific = {}
                    n = 0
                    for allele, ancestors in definitions.items():
                        if ancestors <= variants[phase]: # definition is subset of variants
                            if len(ancestors) > n: # More specific
                                most_specific = allele
                                n = len(ancestors)
                    if len(calling[phase]) == 0: calling[phase].append(set())
                    if n == 0: # No allele definition found
                        # TODO allow this?
                        calling[phase][-1] |= variants[phase]
                        break
                    # Found, replace variants with allele
                    variants[phase] -= definitions[most_specific]
                    calling[phase][-1].add(most_specific)
                calling[phase].append({"CYP2D6*1",})
            # Return valid calling
            # if not valid_calling(sample, calling, homozygous, cont_graph, eq_graph, ov_graph, definitions):
            #     raise Exception("Invalid calling generated")
            representation = calling_to_repr(calling, cont_graph, functions, **detail_from_level(detail_level), reorder=True)
            representation = f"{'+'.join(representation['A'])}/{'+'.join(representation['B'])}"
            # Only return unique representations (duplicates can be generated because of detail level)
            if representation in unique:
                continue
            unique.add(representation)
            yield calling
            # TODO valid?
            if r == 1 and len(het_variants) == 2: # only one combination possible
                break

def order_callings(calling, functions, no_default=True, shortest=True, no_uncertain=True):
    """Order alternative callings by clinical relevance.
    
    no_default: prefer callings with not all on one side, (e.g. *X/*Y over *1+*X/*Y)
    shortest: Prefer shorter callings, equal to without '+' (e.g. *Q/*Z over *X+*Y/*Z)
    no_uncertain: Prefer certain function over uncertain ones.
    
    TODO prefer worse function over better one?

    TODO make this a comparison function?
    """
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
                    

def star_allele_calling_all(samples, nodes, edges, functions, supremals, reference, phased=True, detail_level=0, reorder=True):
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
        sep_callings = separate_callings(callings, cont_graph, functions)
        representations = {}
        # TODO make this into function
        # TODO allow for keeping of multiple alternative representations?
        # test_i = 0
        for sample, calling in sep_callings.items():
            # if sample == "NA19174": continue
            # if sample != "NA19174": continue
            # test_i += 1
            # if test_i > 14:
            #     exit()
            if calling['A'] == calling['B']: # Already phased (homozygous)
                representations[sample] = calling_to_repr(calling, cont_graph, functions, **detail_from_level(detail_level), reorder=reorder)
                print(sample, "(hom)")
                print(f"{'+'.join(representations[sample]['A'])}/{'+'.join(representations[sample]['B'])}")
                print()
                continue
            # All homozygous alleles for the current sample
            homozygous = set([allele for alleles in callings[sample]['hom'] for allele in alleles if allele != "CYP2D6*1"])
            # Generate unique valid alternative callings
            alternatives = generate_alternative_callings_bottom_up(sample, homozygous, cont_graph, eq_graph, ov_graph, functions, detail_level)
            print(sample, "(alt)")
            alternatives = list(alternatives)
            alternatives.sort(key=lambda a: order_callings(a, functions))
            for alternative in alternatives:
                representation = calling_to_repr(alternative, cont_graph, functions, **detail_from_level(detail_level), reorder=True)
                representation = f"{'+'.join(representation['A'])}/{'+'.join(representation['B'])}"# ({order_callings(alternative, functions)})"            
                print(representation)
            # Select the most relevant alternative
            if len(alternatives) > 0:
                preferred = alternatives[0] 
            else:
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
    cores = []
    if match in cont_graph.nodes(): # Only have to check contained since everything is connected to the core
        for node, _ in cont_graph.in_edges(match):
            if find_type(node) == Type.CORE: # Found core contained in match, don't traverse further
                cores.append(node)
            elif find_type(node) == Type.SUB: # Traverse further
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
        if edge[0] == sample and find_type(edge[1]) in (Type.VAR, Type.P_VAR):
            variants.append(edge[1])
        elif edge[1] == sample and find_type(edge[0]) in (Type.VAR, Type.P_VAR):
            variants.append(edge[0])
    if len(variants) == 0: # No variants found
        return {}
    # Check the relevance of each variant
    variants_relevance = {}
    for variant in variants:
        # Check if variant possibly interferes with any allele (overlaps with supremal)
        interferes = any([overlap(supremals[variant], supremals[allele]) for allele in alleles if allele != "CYP2D6*1"])
        # Find the impact of the variant
        impact = functions[variant]
        if find_type(variant) == Type.P_VAR: 
            severity = severity_GO(impact)
        elif find_type(variant) == Type.VAR:
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
                if find_type(match2) == Type.SUB: # Cores can be contained in subs
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
                    cores = find_core_traversal(allele, cont_graph) # Cores of suballele
                    representation[phase].extend(cores)
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