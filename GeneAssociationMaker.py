import json


def make_associations(ecoli_model_path):
    with open(ecoli_model_path, 'r') as f:
        ecoli_model = json.load(f)
    """ 
    model keys: ['metabolites', 'reactions', 'genes', 'id', 'compartments', 'version']
    ecoli_model['metabolites'] = list of metabolites containing keys:
            ['id', 'name', 'compartment', 'charge', 'formula', 'notes', 'annotation']
    ecoli_model['reactions'] = list of reactions containing keys:
            ['id', 'name', 'metabolites', 'lower_bound', 'upper_bound', 'gene_reaction_rule', 'notes', 'annotation']
            reaction['notes'] = {'original_bigg_ids': ['2AGPGAT180']}
    ecoli_model['genes'] = list of genes containing keys:
            ['id', 'name', 'notes', 'annotation']
    ecoli_model['id'] = iAF1260b
    ecoli_model['compartments'] = {'c': 'cytosol', 'e': 'extracellular space', 'p': 'periplasm'}
    ecoli_model['version'] = 1
    """
    ecoli_reactions = ecoli_model['reactions']
    print(ecoli_reactions[120])
    reactions_genes_list = dict()
    for (reaction_name, reaction_dict) in reactions:
        reaction_genes_assoc = reaction_dict.gene_product_association
        all_reaction_genes = get_genes(reaction_genes_assoc)
        reactions_genes_list[reaction_name] = all_reaction_genes
