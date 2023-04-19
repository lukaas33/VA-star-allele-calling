import requests
import pickle
import os
import os.path
import vcf
import algebra as va
import warnings

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
    """
    # TODO find proper API, is now hardcoded on web request
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
        response = response.replace('\n', '').replace('\r', '') # Remove newlines
    return response

def pharmvar_api_get(target):
    """ Call the Pharmvar API and return the data as a Python object. 

    More information can be found here:
    https://www.pharmvar.org/documentation
    """
    # TODO take version number into account
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

def cache_set(data, name):
    """Store data to avoid generating/retrieving it too often."""
    foldername = "temp"
    filename = f"{foldername}/{name}.pkl"
    if not os.path.exists(foldername):
        os.mkdir(foldername)
    with open(filename, 'wb') as file:
        pickle.dump(data, file)

def cache_get(name):
    """Get stored data"""
    # TODO implement some retention method.
    foldername = "temp"
    filename = f"{foldername}/{name}.pkl"
    if os.path.isfile(filename):
        with open(filename, 'rb') as file:
            return pickle.load(file)
    raise ValueError("Cached file not found")

def pharmvar_get(target):
    """ Get data from api and save it unless it is already stored.

    Reduces waiting time when experimenting.
    Stores it by target and the output shouldn't change that often.
    """
    name = target.replace('/', '_').replace('*', '')
    try: 
        return cache_get(name)
    except:
        data = pharmvar_api_get(target)
        cache_set(data, name)
        return data

def parse_samples(directory, reference, phased=False):
    """ Parse sample VCF files as variant objects."""
    # TODO detect if phased
    # TODO handle mixed phasing
    # TODO handle unphased samples
    samples = {}
    for filename in os.listdir(directory):
        sample_name = filename.split('.')[0]
        if phased:
            allele = {"A": {}, "B": {}} # First and second allele
        else:
            allele = {"hom": {}, "het": {}, "all": {}} # Homozygous, heterozygous and all variants
        # TODO check that no records are double
        with open(os.path.join(directory, filename), 'r') as file:
            reader = vcf.Reader(file)
            # QUESTION: what are the filter, quality, format, info fields?
            # QUESTION: what is the format of alt
            for record in reader:
                # Validate if reference of vcf files is on the reference sequence
                if record.REF != reference[record.start:record.end]:
                    raise ValueError("Reference sequence does not match")
                # Check ALT field. Can have multiple values for different phases, etc.
                if len(record.ALT) > 1: 
                    raise ValueError("Multiple ALT alleles not supported") # TODO handle different alt values?
                # Check multiple samples
                if len(record.samples) > 1:
                    raise ValueError("Multiple samples not supported") # TODO handle
                # Create variant with Zero based half-open positions
                variant = va.Variant(record.start, record.end, record.ALT[0].sequence) 
                hgvs = va.variants.to_hgvs([variant]) # TODO add prefix and reference (but allow differentiation from pharmvar variants)
                phasing = record.samples[0].data[0]
                # Store variant in correct alleles
                if phased:
                    if '|' not in phasing: 
                        raise ValueError("Sample is not phased")
                    # Add to correct phased allele
                    # 0|1 or 1|0 shows that the variant is one of the two alleles, and can be used to group variants into two alleles
                    # 1|1 
                    for p, phase in zip("AB", phasing.split('|')):
                        if phase == '1':
                            allele[p][hgvs] = variant
                else:
                    if '|' in phasing:
                        raise ValueError("Sample is phased")
                    # Add to unphased allele
                    # 1/1 is homozygous, here perfect phasing information is available
                    # 1/0 or 0/1 is heterozygous and no phasing information is available
                    allele["all"][hgvs] = variant
                    if phasing == "1/1":
                        allele["hom"][hgvs] = variant
                    elif phasing == "1/0" or phasing == "0/1": # TODO handle 1/0 differently?
                        allele["het"][hgvs] = variant
                    else:
                        raise ValueError("Unknown phasing", phasing)    
            # Add alleles to sample
            for p in allele:
                samples[f"{sample_name}_{p}"] = allele[p]
    return samples