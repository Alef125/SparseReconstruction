# import xml.etree.ElementTree as ETree
import json
import numpy as np


def get_fields(fields_str):
    repaired_str = repr(fields_str)[1:-1]
    print(repaired_str)
    fields = repaired_str.split("\\t")
    return fields


# def parse_xml(filename):
#     with open(filename, 'r') as f:
#         data = json.load(f)
#         print(data)
    # tree = ETree.parse(filename)
    # sbml_root = tree.getroot()  # <sbml ...>  {model}  </sbml>
    # model_root = sbml_root[0]  # <mode ...>  {listOfUnitDefinitions, listOfObjectives, listOfParameters,
    # #                                         listOfCompartments, listOfSpecies, listOfGeneProducts,
    # #                                         listOfReactions}  </model>
    # child = model_root[0]
    # print(child.tag)
    # print(child.attrib)
    # for child2 in child:
    #     print(child2.tag)
    #     print(child2.attrib)

# file_path = "../Data/Knock-Out Experiments - EcoCyc/LB enriched/LB_enriched_-_NG.txt"
# with open(file_path, 'r') as f:
#     lines = f.readlines()
#     print(lines[0])
#     data_fields = get_fields(lines[0][:-1])
#
# print(data_fields)


def is_reaction_knocked_out(out_gene, rec_associated_genes):
    if not rec_associated_genes:  # empty dict
        return False
    assoc_kind = list(rec_associated_genes.keys())[0]
    if assoc_kind == "GPARef":
        if rec_associated_genes["GPARef"] == out_gene:
            return True
        return False
    elif assoc_kind == "GPAAnd":
        for and_part in rec_associated_genes["GPAAnd"]:
            if is_reaction_knocked_out(out_gene, and_part):
                return True
        return False
    elif assoc_kind == "GPAOr":
        for or_part in rec_associated_genes["GPAOr"]:
            if not is_reaction_knocked_out(out_gene, or_part):
                return False
        return True
    else:
        raise Exception


def get_knocked_out_reactions(knocked_out_gene, gene_associations):
    knocked_out_reactions = []
    for reaction_id, rec_associated_genes in gene_associations.items():
        if is_reaction_knocked_out(knocked_out_gene, rec_associated_genes):
            knocked_out_reactions.append(reaction_id)
    return knocked_out_reactions


def get_fields_name(fields_str):
    repr_str = repr(fields_str)[1:-3]
    fields = repr_str.split("\\t")
    return fields


def get_knock_out_genes_names(line_str, fields):
    repr_str = repr(line_str)[1:-3]
    line_parts = repr_str.split("\\t")
    line_dict = dict(zip(fields, line_parts))
    return line_dict


def get_gene_main_name(names_str):
    gene_names = names_str.split(" // ")
    return gene_names[0]


def knock_outs_gene_names(knock_out_lines):
    data_fields = get_fields_name(knock_out_lines[0])
    num_data_lines = len(knock_out_lines) - 1
    growth_genes_list = []
    for i in range(1, num_data_lines):
        growth_gene_dict = get_knock_out_genes_names(knock_out_lines[i], data_fields)
        growth_gene_names_str = growth_gene_dict["Names"]
        growth_gene_main_name = get_gene_main_name(growth_gene_names_str)
        growth_genes_list.append(growth_gene_main_name)
    return growth_genes_list


def get_gene_mask(knocked_out_reactions, all_reactions_ids_list):
    num_reactions = len(all_reactions_ids_list)
    gene_mask = np.ones(num_reactions)
    print(knocked_out_reactions)
    print(all_reactions_ids_list[:30])
    for i in range(num_reactions):
        reaction_id = all_reactions_ids_list[i]
        if reaction_id in knocked_out_reactions:
            gene_mask[i] = 0
    return gene_mask


def make_knock_out_masks(all_reactions_ids_list, knock_outs_info, associations):
    for knock_out_gene in knock_outs_info:
        knocked_out_reactions = get_knocked_out_reactions(knock_out_gene, associations)
        the_gene_mask = get_gene_mask(knocked_out_reactions, all_reactions_ids_list)


def save_knock_out_bounds(all_reactions_ids_list, gene_associations_path, knock_outs_path, save_path):
    with open(gene_associations_path, 'r') as f:
        associations = json.load(f)
    with open(knock_outs_path, 'r') as f:
        lines = f.readlines()
        knock_outs_info = knock_outs_gene_names(lines)

    make_knock_out_masks(all_reactions_ids_list, knock_outs_info, associations)
