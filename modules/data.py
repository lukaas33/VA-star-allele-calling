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

def parse_samples(reference):
    """ Parse sample VCF files as variant objects."""
    directory = "data/samples"
    samples = {}
    for filename in os.listdir(directory):
        name = filename.split('.')[0]
        phased_allele = [[], []]
        with open(os.path.join(directory, filename), 'r') as file:
            reader = vcf.Reader(file)
            # QUESTION: what are the filter, quality, format, info fields?
            # QUESTION: what is the format of alt
            for record in reader:
                # Validate if reference of vcf files is on the reference sequence
                if record.REF != reference[record.start:record.end]:
                    raise ValueError("Reference sequence does not match")
                if len(record.ALT) > 1: # TODO handle different alt values
                    # Can have multiple values for different phases, etc.
                    raise ValueError("Multiple ALT alleles not supported")
                # Zero based half-open
                variant = va.Variant(record.start, record.end, record.ALT[0].sequence) 
                if len(record.samples) > 1:
                    raise ValueError("Multiple samples not supported (how to interpret this?))")
                phasing = record.samples[0].data[0]
                # Add to correct phased allele
                # 0|1 is the same as 1|0 between samples (not within sample) 
                for i, phase in enumerate(phasing.split('|')):
                    if phase == '1':
                        hgvs = va.variants.to_hgvs([variant]) # TODO add prefix and reference (but allow differentiation from pharmvar variants)
                        phased_allele[i].append((hgvs, variant))
        # Add phased alleles to samples
        for i in range(2):
            phased_name = name + "AB"[i]
            if len(phased_allele[i]) == 0: # Empty alleles will be treated as *1
                pass 
            samples[phased_name] = phased_allele[i]
        # Add unphased 'allele' to samples
        # Also store double variants since this is relevant for calling
        all = []
        double = []
        for v in phased_allele[0] + phased_allele[1]:
            if v in all: # present twice
                double.append(v)
                continue
            for v2 in all: # overlaps so cannot be in same allele QUESTION is this correct?
                s, e = v[1].start, v[1].end
                s2, e2 = v2[1].start, v2[1].end
                if max(s, s2) <= min(e, e2): # overlap
                    double.append(v)
                    break
            else: # no overlap and not already present
                all.append(v)
        samples[name] = all
        if len(double) > 0:
            samples[name + '+'] = double 
    return samples