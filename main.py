import requests
import pickle
import os.path
import os
import algebra as va

def api_get(url, params={}):
    """ Wrapper for making API calls. 
    """
    result = requests.get(url, params=params)
    if result.status_code != 200:
        raise Exception(f"API query to {url} was invalid")
    try:
        response = result.json()
    except ValueError:
        response = result.text
    return response

def reference_get():
    """ Get reference sequence

    TODO find proper API, is now hardcoded on webqrequest
    """
    url = "https://www.ncbi.nlm.nih.gov/sviewer/viewer.cgi?tool=portal&save=file&log$=seqview&db=nuccore&report=fasta&id=568815576&&ncbi_phid=null"
    foldername = "data"
    filename = f"{foldername}/NC000022.11.fasta"
    if not os.path.isfile(filename):
        response = api_get(url)
        if not os.path.exists(foldername):
            os.mkdir(foldername)
        with open(filename, 'w') as file:
            file.write(response)
    with open(filename, 'r') as file:
        file.readline() # Skip header
        response = file.read()
    # TODO needed to remove newlines?
    return response


def pharmvar_api_get(target):
    """ Call the Pharmvar API and return the data as a Python object. 

    More information can be found here:
    https://www.pharmvar.org/documentation
    """
    baseurl = "https://www.pharmvar.org/api-service/"
    url = baseurl + target
    default_params = {
        "exclude-sub-alleles": "false", # Get suballeles
        "include-reference-variants": "false", # Don't include reference variants
        "include-retired-alleles": "false", # Don't include retired definitions
        "include-retired-reference-sequences": "false", # Don't include retired references
        "reference-collection": "GRCh38", # Use newest HC reference
        #  Dont'filter by "position" or "reference-sequence"
    }    
    return api_get(url, default_params)

def pharmvar_get(target):
    """ Get data from api and save it unless it is already stored.

    Reduces waiting time when experimenting.
    Stores it by target and the output shouldn't change that often.
    TODO implement some retention method.
    """
    foldername = "temp"
    filename = f"{foldername}/{target.replace('/', '_').replace('*', '')}.pkl"
    if os.path.isfile(filename):
        with open(filename, 'rb') as file:
            return pickle.load(file)
    else:
        data = pharmvar_api_get(target)
        if not os.path.exists(foldername):
            os.mkdir(foldername)
        with open(filename, 'wb') as file:
            pickle.dump(data, file)
        return data
        

# Get the reference sequence relevant for the (current) gene of interest
reference_sequence = reference_get()


# List genes as symbols in Pharmvar
genes = pharmvar_get("genes/list") 
# All information associated with the (current) gene of interest
gene = pharmvar_get("genes/CYP2D6") 
# Group suballeles by core alleles and index by the star-allele notation of the core allele
corealleles = {allele["alleleName"]: allele for allele in gene["alleles"] if allele["alleleType"] == "Core"} 
suballeles = {coreallele: [sub_allele for sub_allele in gene["alleles"] if sub_allele["coreAllele"] == coreallele] for coreallele in corealleles.keys()}
# All information associated with the (current) allele of interest
coreallele = corealleles["CYP2D6*3"]
suballele = suballeles["CYP2D6*3"]
# Get the variants contained in the sub- and core allele of interest as HGVS notation
coreallele_variants = [variant["hgvs"] for variant in coreallele["variants"]]
for suballele in suballele:
    suballele_variants = [variant["hgvs"] for variant in suballele["variants"]]
    # Convert HGVS notation to sequence notation
    s_coreallele_variants = []
    for variant in coreallele_variants:
        s_coreallele_variants += va.variants.parse_hgvs(variant) 
    s_suballele_variants = []
    for variant in suballele_variants:
        s_suballele_variants += va.variants.parse_hgvs(variant) 
    # Find relation between core and suballeles
    # TODO should you find the relation between the individual variants or the alleles? (right now the 2nd)
    relation = va.compare(reference_sequence, s_coreallele_variants, s_suballele_variants)
    # Expect containment or equivalence
    if relation not in (va.Relation.EQUIVALENT, va.Relation.IS_CONTAINED):
        print(coreallele["alleleName"], coreallele_variants)
        print(suballele["alleleName"], suballele_variants)
        raise ValueError(f"Unexpected relationship between core allele {coreallele['alleleName']} and suballele {suballele['alleleName']}: {relation}")

# TODO what is a variant group
# TODO Check if names follow a logical format
# TODO Check if core allele equivalent or contained in each suballele