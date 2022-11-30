#Input functions used by rule all

def get_fastqc_output(species_list, timing):
    # Get the sample id of all samples of a given species
    paths = []
    for species in species_list:
        type_ = str(type(config[f'dataset_{species}']))
        if type_ == "<class 'dict'>":
            ids = list(config[f'dataset_{species}'].keys())
        elif type_ == "<class 'list'>":
            ids = config[f'dataset_{species}']
        path = ['logs/fastqc/{}/{}/{}.log'.format(timing, species, id) for id in ids] # timing: before of after trimming
        paths.append(path)
    # Return the fastqc log file paths of all species samples
    return [item for sublist in paths for item in sublist]


def get_coco_output(species_list, coco_mode):
    # Get the sample id of all samples of a given species
    extension = {'coco_cc': 'tsv', 'coco_cb': 'bedgraph'}
    ext = extension[coco_mode]
    paths = []
    for species in species_list:
        type_ = str(type(config[f'dataset_{species}']))
        if type_ == "<class 'dict'>":
            ids = list(config[f'dataset_{species}'].keys())
        elif type_ == "<class 'list'>":
            ids = config[f'dataset_{species}']
        path = ['results/{}/{}/{}.{}'.format(coco_mode, species, id, ext) for id in ids] # coco_mode: coco_cc or coco_cb
        paths.append(path)
    # Return the coco output file paths of all species samples
    return [item for sublist in paths for item in sublist]

    

# Function for wildcard constraints

def join_list(l, subl, remove=False):
    # From a list l, return a string of all items in subl joined by '|'
    small_list = [a for a in l if a in subl]
    # If we want to remove (instead of only keeping) items of subl from l
    if remove==True:
        small_list = [a for a in l if a not in subl]
    return "{}".format("|".join(small_list))



# Input functions for general rules

def get_species_genome(species):
    # Get the fasta of the genome of a given species
    species = str(species)
    if species == 'tetrahymena_thermophila':
        path = rules.download_tetrahymena_genome.output.genome
    elif 'saccharomyces' in species:
        path = rules.download_yeast_genome.output.genome
    elif species in ['homo_sapiens', 'mus_musculus']:
        path = rules.download_mammal_genome.output.genome
    return path


def get_species_gtf(species):
    # Get the gtf of the genome of a given species
    species = str(species)
    if species == 'tetrahymena_thermophila':
        path = rules.download_tetrahymena_gtf.output.gtf
    elif 'saccharomyces' in species:
        path = rules.download_yeast_gtf.output.gtf
    elif species == 'homo_sapiens':
        path = rules.download_human_gtf.output.gtf
    elif species == 'mus_musculus':
        path = rules.download_mouse_gtf.output.gtf
    return path


def get_coco_cc_output_per_species(species):
    # Get the sample id of all samples of a given species
    type_ = str(type(config[f'dataset_{species}']))
    if type_ == "<class 'dict'>":
        ids = list(config[f'dataset_{species}'].keys())
    elif type_ == "<class 'list'>":
        ids = config[f'dataset_{species}']
    path = ['results/coco_cc/{}/{}.tsv'.format(species, id) for id in ids] 
    # Return the coco output file paths of all the samples of a given species
    return path
