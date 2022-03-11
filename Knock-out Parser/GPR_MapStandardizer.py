"""
GPR_MapStandardizer.py

This script, converts common formats of GPR map or Genes Association Data into the standard form used in the
main preprocessing code.

The standard output is a dictionary (.json file) in the format of {gene_id : reactions_list},
in the which the reactions_list is the list of reactions which would be shut down if we knock-out the gene
defined with the id gene_id.
"""

import json
import pandas as pd


def get_all_associated_genes(gpr_expression: dict) -> list:
    """
    :param gpr_expression: A dict in a format like:
                    {"GPAOr":[{"GPARef":"ulaD"},{"GPARef":"sgbH"}]}
    :return: All genes existing in the gpr_expression
    """
    all_associated_genes = []
    for id_key, expr in gpr_expression.items():
        if id_key == "GPARef":
            all_associated_genes.append(expr)
        else:
            for gpr_item in expr:
                all_associated_genes = all_associated_genes + get_all_associated_genes(gpr_item)
    return list(set(all_associated_genes))


def is_shutter_gene(gpr_expression: dict, gene_id: str) -> bool:
    """

    :param gpr_expression: A dict in a format like:
                    {"GPAOr":[{"GPARef":"ulaD"},{"GPARef":"sgbH"}]}
    :param gene_id: A gene existing in the gpr_expression
    :return: True if the gene's knock-out will shut the expression, and False otherwise
    """
    for id_key, expr in gpr_expression.items():
        if id_key == "GPARef":
            return expr == gene_id
        elif id_key == "GPAAnd":
            for item_expression in expr:  # Shutting at least one item
                if is_shutter_gene(gpr_expression=item_expression, gene_id=gene_id):
                    return True
            return False
        elif id_key == "GPAOr":
            for item_expression in expr:  # Shutting all item
                if not is_shutter_gene(gpr_expression=item_expression, gene_id=gene_id):
                    return False
            return True
        else:
            print("Your GPR format with key " + id_key + " is not standard")
            raise Exception


def find_knocker_out_genes_from_gpa(gpr_expression: dict) -> list:
    """
    :param gpr_expression: A dict in a format like:
                    {"GPAOr":[{"GPARef":"ulaD"},{"GPARef":"sgbH"}]}
    :return: A list of gene which their knocking-out will knock-out the entire expression
    """
    all_associated_genes = get_all_associated_genes(gpr_expression)
    knocker_out_genes = []
    for gene_id in all_associated_genes:
        if is_shutter_gene(gpr_expression=gpr_expression, gene_id=gene_id):
            knocker_out_genes.append(gene_id)
    return knocker_out_genes


def convert_rule_to_gpa(gpr_rule: str) -> dict:
    """
    :param gpr_rule: A str in a format like:
            ( BSU29690 and BSU29700 and BSU29710 ) or ( BSU08060 and BSU08070 and BSU08090 ) or BSU26640
    :return: A dict in a format like:
                    {'GPAOr': [{'GPARef': 'BSU28560'},
                               {'GPAOr': [{'GPARef': 'BSU18250'},
                                          {'GPAOr': [{'GPARef': 'BSU10360'},
                                                     {'GPARef': 'BSU10270'}]
                                          }]
                                }]
                    }
    """

    # print(gpr_expression)
    if gpr_rule[0:2] == '( ':
        close_parentheses_idx = gpr_rule.index(")")

        if len(gpr_rule) == close_parentheses_idx + 1:  # gpr_expression = "( ... )"
            return convert_rule_to_gpa(gpr_rule[2:-2])
        elif gpr_rule[close_parentheses_idx+1:close_parentheses_idx+5] == ' or ':
            operator = "GPAOr"
            second_term_start_idx = close_parentheses_idx + 5
        elif gpr_rule[close_parentheses_idx+1:close_parentheses_idx+6] == ' and ':
            operator = "GPAAnd"
            second_term_start_idx = close_parentheses_idx + 6
        else:
            print("Undefined Operator")
            raise Exception

        first_parentheses_term = gpr_rule[2:close_parentheses_idx - 1]
        first_member = convert_rule_to_gpa(first_parentheses_term)
        second_term = gpr_rule[second_term_start_idx:]
        second_member = convert_rule_to_gpa(second_term)
        return {operator: [first_member, second_member]}

    if ' ' in gpr_rule:
        space_idx = gpr_rule.index(' ')
        if gpr_rule[space_idx:space_idx + 4] == ' or ':
            operator = "GPAOr"
            second_term_start_idx = space_idx + 4
        elif gpr_rule[space_idx:space_idx + 5] == ' and ':
            operator = "GPAAnd"
            second_term_start_idx = space_idx + 5
        else:
            print("Undefined Operator")
            raise Exception
        first_parentheses_term = gpr_rule[:space_idx]
        first_member = convert_rule_to_gpa(first_parentheses_term)
        second_term = gpr_rule[second_term_start_idx:]
        second_member = convert_rule_to_gpa(second_term)
        return {operator: [first_member, second_member]}
    else:
        return {"GPARef": gpr_rule}

# ######################################################################################################


class GPRMapConverter:
    def __init__(self,
                 gene_assoc_data_filepath: str,
                 input_reactions_nomenclature: str = None,
                 output_reactions_nomenclature: str = None,
                 translation_filepath: str = None):
        """
        :param gene_assoc_data_filepath: The file path for gene_assoc_data.json
        :param input_reactions_nomenclature: The column name in the translation_file
                                             corresponding to the reaction names in the gene_assoc_data.json
        :param output_reactions_nomenclature: The column name in the translation_file
                                              corresponding to the desired reaction names

                                  ** Note **: This argument *should* be the nomenclature used in your
                                              Internal_Rxns_Bounds.csv and EX_Rxns_Bounds.csv files.

        :param translation_filepath: The file path for translation_file, which contains all different names
                                     for reactions
        """
        self.gene_assoc_data_filepath = gene_assoc_data_filepath
        self.gene_assoc_data = None
        self.read_genes_assoc_data()
        # #############################################################
        self.input_reactions_nomenclature = input_reactions_nomenclature
        self.output_reactions_nomenclature = output_reactions_nomenclature
        # ##############################################
        self.translation_filepath = translation_filepath
        self.translation_dict = None
        self.make_translation_dict()

    def read_genes_assoc_data(self):
        """
        This method, reads gene_assoc_data from self.gene_assoc_data_filepath
        """
        with open(self.gene_assoc_data_filepath, 'r') as json_file:
            self.gene_assoc_data = json.load(json_file)

    def make_translation_dict(self):
        """
        This method, reads translation_file.csv from self.translation_filepath,
            then makes the self.translation_dict dictionary with keys as reactions names in
            self.input_reactions_nomenclature and values as reactions names in self.output_reactions_nomenclature
        """
        if self.translation_filepath:
            self.translation_dict = {}
            translation_file = pd.read_csv(self.translation_filepath)
            if self.input_reactions_nomenclature:
                if self.input_reactions_nomenclature not in translation_file.columns:
                    print("input_reactions_nomenclature does not exist in translation_file")  # ToDo: to Exception
                    raise Exception
            if self.output_reactions_nomenclature:
                if self.output_reactions_nomenclature not in translation_file.columns:
                    print("output_reactions_nomenclature does not exist in translation_file")
                    raise Exception
            for index, row in translation_file.iterrows():
                key = row[self.input_reactions_nomenclature]
                value = row[self.output_reactions_nomenclature]
                self.translation_dict[key] = value

    def make_genes_to_reactions_ko_dict(self, filepath_to_save: str, gpr_type: str = "GPA"):
        """
        :param filepath_to_save: The path to save the final dictionary as a .json file.
        :param gpr_type: Defines the type of self.gene_assoc_data:
                         "GPA": A dict in the format of e.g.:
                         {"R_KG6PDC": {"GPAOr":[{"GPARef":"ulaD"},{"GPARef":"sgbH"}]}}
                         "Rule": A dict in the format of e.g.:
                         {"ACTD2": "( BSU29690 and BSU29700 and BSU29710 ) or ( BSU08060 and BSU08070 and BSU08090 )"}
        :return: Saves the standard .json GPR file.
                 If self.translation_dict is None, the reactions ids will be the same with the ids in the
                 input self.gene_assoc_data data. Otherwise, they will be translated from
                 self.input_reactions_nomenclature into self.output_reactions_nomenclature.
        """
        genes_to_reactions_ko_dict = {}
        for rxn_name, gpr_expression in self.gene_assoc_data.items():
            # ############### Translating the reaction name ###################
            output_rxn_name = rxn_name
            if self.translation_dict:
                try:
                    output_rxn_name = self.translation_dict[rxn_name]
                except KeyError:
                    print("The reaction " + rxn_name + " does not exist in the translation file")
                    raise Exception
            # ############### Find knocker-out genes for this reaction ###################
            if gpr_type == "GPA":
                gpr_gpa_format = gpr_expression
            elif gpr_type == "Rule":
                gpr_gpa_format = convert_rule_to_gpa(gpr_rule=gpr_expression)
            else:
                print("gpr_type is not defined correctly")
                raise Exception
            knocker_out_genes = find_knocker_out_genes_from_gpa(gpr_expression=gpr_gpa_format)
            if "" in knocker_out_genes:
                knocker_out_genes.remove("")
            # ############ Add this reaction to the corresponding knocker-out genes ################
            for knocker_out_gene in knocker_out_genes:
                if knocker_out_gene in genes_to_reactions_ko_dict.keys():
                    genes_to_reactions_ko_dict[knocker_out_gene].append(output_rxn_name)
                else:
                    genes_to_reactions_ko_dict[knocker_out_gene] = [output_rxn_name]
            # #################################################################
        with open(filepath_to_save, 'w', encoding='utf-8') as f:
            json.dump(genes_to_reactions_ko_dict, f, ensure_ascii=False, indent=4)


# obj = GPRMapConverter(gene_assoc_data_filepath="../Data/Palsson B.Subtilis Reconstruction/Genes Associations.json",
#                       input_reactions_nomenclature='R_Rxn',
#                       output_reactions_nomenclature='Rxn',
#                       translation_filepath="../Data/Palsson B.Subtilis Reconstruction/Reactions Translation File.csv")
obj = GPRMapConverter(gene_assoc_data_filepath="../../Data/Palsson B.Subtilis Reconstruction/B_Subtilis GPR rules.json",
                      input_reactions_nomenclature='Base id',
                      output_reactions_nomenclature='BiGG id',
                      translation_filepath="../../Data/Palsson B.Subtilis Reconstruction/B_Subtilis Rxn Translation.csv")
obj.make_genes_to_reactions_ko_dict(filepath_to_save="../../Data/Palsson B.Subtilis Reconstruction/Organism GPR.json",
                                    gpr_type="Rule")
