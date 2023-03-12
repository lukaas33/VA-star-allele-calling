import requests
import pickle
import os
import os.path
import vcf
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

def parse_samples():
    """ Parse sample VCF files as variant objects."""
    directory = "data/samples"
    samples = {}
    for filename in os.listdir(directory):
        variants = []
        with open(os.path.join(directory, filename), 'r') as file:
            reader = vcf.Reader(file)
            # TODO validate input, on reference sequence?
            # QUESTION: what are the filter, quality, format, info fields?
            # QUESTION: what is the format of alt
            for record in reader:
                if len(record.ALT) > 1: # TODO handle different alt values
                    raise ValueError("Multiple ALT alleles not supported")
                # TODO check if positions are correct
                variant = va.Variant(record.start, record.end, record.ALT[0].sequence) 
                variants.append(variant)
        name = filename.split('.')[0]
        samples[name] = variants 
        break # TODO remove
    return samples