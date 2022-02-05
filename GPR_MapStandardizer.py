"""
GPR_MapStandardizer.py

This script, converts common formats of GPR map or Genes Association Data into the standard form used in the
main preprocessing code.
"""

import json
import pandas as pd


def find_knocker_out_genes_from_gpr(gpr_expression: dict) -> list:
    """
    :param gpr_expression:
    :return: A list of gene which their knocking-out will knock-out the entire expression
    """
    # ToDo: Break down this huge function
    if gpr_expression.keys():
        if 'GPARef' in gpr_expression.keys():
            if gpr_expression['GPARef'] != 'SPONTANEOUS':
                return [gpr_expression['GPARef']]
        elif 'GPAOr' in gpr_expression.keys():
            return []  # A single gene knock-out in OR expression doesn't shut the reaction
        elif 'GPAAnd' in gpr_expression.keys():
            knocker_out_genes = []
            for peptide_complex in gpr_expression['GPAAnd']:
                for id_key, expr in peptide_complex.items():
                    if id_key == 'GPARef':
                        knocker_out_genes.append(expr)
                    elif id_key == 'GPAAnd':
                        [knocker_out_genes.append(ref_expr['GPARef']) for ref_expr in expr]
                    else:
                        raise Exception
            return knocker_out_genes
        else:
            raise Exception  # ToDo: Format checking is not complete here
    return []


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
        :param translation_filepath: The file path for translation_file, which contains all different names
                                     for reactions
        """
        self.gene_assoc_data_filepath = gene_assoc_data_filepath
        self.gene_assoc_data = None
        self.read_genes_assoc_data()
        self.input_reactions_nomenclature = input_reactions_nomenclature
        self.output_reactions_nomenclature = output_reactions_nomenclature
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
                    print("output_reactions_nomenclature does not exist in translation_file")  # ToDo: to Exception
                    raise Exception
            for index, row in translation_file.iterrows():
                key = row[self.input_reactions_nomenclature]
                value = row[self.output_reactions_nomenclature]
                self.translation_dict[key] = value

    def make_genes_to_reactions_ko_dict(self, filepath_to_save):
        genes_to_reactions_ko_dict = {}
        for rxn_name, gpr_expression in self.gene_assoc_data.items():
            # ############### Translating the reaction name ###################
            output_rxn_name = rxn_name
            if self.translation_dict:
                output_rxn_name = self.translation_dict[rxn_name]  # ToDo: Exception if not in translation_dict.keys()
            # ############ Add this reaction to the corresponding knocker-out genes ################
            knocker_out_genes = find_knocker_out_genes_from_gpr(gpr_expression=gpr_expression)
            for knocker_out_gene in knocker_out_genes:
                if knocker_out_gene in genes_to_reactions_ko_dict.keys():
                    genes_to_reactions_ko_dict[knocker_out_gene].append(output_rxn_name)
                else:
                    genes_to_reactions_ko_dict[knocker_out_gene] = [output_rxn_name]
            # #################################################################
        with open(filepath_to_save, 'w', encoding='utf-8') as f:
            json.dump(genes_to_reactions_ko_dict, f, ensure_ascii=False, indent=4)


obj = GPRMapConverter(gene_assoc_data_filepath="../Data/Palsson B.Subtilis Reconstruction/Genes Associations.json",
                      input_reactions_nomenclature='R_Rxn',
                      output_reactions_nomenclature='Rxn',
                      translation_filepath="../Data/Palsson B.Subtilis Reconstruction/Reactions Translation File.csv")
obj.make_genes_to_reactions_ko_dict(filepath_to_save="../Data/Palsson B.Subtilis Reconstruction/Organism GPR.json")
