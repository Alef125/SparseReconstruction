"""
GenesKO_Standardizer.py
This code, converts common formats for Knock-Out Essentiality Data, such as .csv, into the standard .json file
used in further processes
"""

import json
import pandas as pd
import warnings


class GenesKOStandardizer:
    def __init__(self,
                 genes_ko_growth_filepath: str,
                 medium_name: str,
                 input_genes_nomenclature: str = None,
                 output_genes_nomenclature: str = None,
                 translation_filepath: str = None):
        """
        :param genes_ko_growth_filepath: The file path for knock-out experiment
        :param medium_name: Name of the media in which the KO experiment has been done
        :param input_genes_nomenclature: The column name in the translation_file
                                         corresponding to the genes names in the genes_ko_growth_filepath.csv
        :param output_genes_nomenclature: The column name in the translation_file
                                          corresponding to the desired genes names.

                              ** Note **: This argument *should* be the nomenclature used in you GPR rules.

        :param translation_filepath: The file path for translation_file, which contains all different names for genes
        """
        self.genes_ko_growth_filepath = genes_ko_growth_filepath
        self.genes_ko_growth_file = None
        self.read_genes_ko_growth_file()
        # ##############################
        self.medium = medium_name
        # #######################################################
        self.input_genes_nomenclature = input_genes_nomenclature
        self.output_genes_nomenclature = output_genes_nomenclature
        # ###############################################
        self.translation_filepath = translation_filepath
        self.translation_dict = None
        self.make_translation_dict()

    def read_genes_ko_growth_file(self):
        """
        This method reads the genes_ko_growth_file.csv from self.genes_ko_growth_filepath.
        That .csv file should have two columns named 'Gene' and 'Growth'.
        :return:
        """
        self.genes_ko_growth_file = pd.read_csv(self.genes_ko_growth_filepath)
        if 'Gene' not in self.genes_ko_growth_file.columns:
            print("No column named Gene is in this file")
            raise Exception  # ToDo: text in Exception
        if 'Growth' not in self.genes_ko_growth_file.columns:
            print("No column named Growth is in this file")
            raise Exception

    def make_translation_dict(self):
        """
        This method, reads translation_file.csv from self.translation_filepath,
            then makes the self.translation_dict dictionary with keys as genes names in
            self.input_genes_nomenclature and values as genes names in self.output_genes_nomenclature
        """
        if self.translation_filepath:
            self.translation_dict = {}
            translation_file = pd.read_csv(self.translation_filepath)
            if self.input_genes_nomenclature:
                if self.input_genes_nomenclature not in translation_file.columns:
                    print("input_reactions_nomenclature does not exist in translation_file")  # ToDo: to Exception
                    raise Exception
            if self.output_genes_nomenclature:
                if self.output_genes_nomenclature not in translation_file.columns:
                    print("output_reactions_nomenclature does not exist in translation_file")
                    raise Exception
            for index, row in translation_file.iterrows():
                key = row[self.input_genes_nomenclature]
                value = row[self.output_genes_nomenclature]
                self.translation_dict[key] = value

    def make_genes_ko_growth_dict(self, filepath_to_save):
        """
        :param filepath_to_save: The path to save the .json file
        :return: Nothing, Saves the .json file
        """
        genes_ko_dicts = []
        for index, row in self.genes_ko_growth_file.iterrows():
            gene_name = row['Gene']
            # ############### Translating the reaction name ###################
            output_gene_name = gene_name
            if self.translation_dict:
                try:
                    output_gene_name = self.translation_dict[gene_name]
                except KeyError:
                    warn_text = "The gene " + gene_name + " does not exist in the translation file"
                    warnings.warn(warn_text)
                    continue
                    # raise Exception
            # ################# Making the KO experiment dict ###################
            genes_ko_dicts.append(
                {'ko_gene_id': output_gene_name,
                 'medium': self.medium,
                 'growth': row['Growth']}
            )
        with open(filepath_to_save, 'w', encoding='utf-8') as f:
            json.dump(genes_ko_dicts, f, ensure_ascii=False, indent=4)


obj = GenesKOStandardizer(
    genes_ko_growth_filepath="../../Data/Palsson B.Subtilis Reconstruction/Genes KO Growth Data.csv",
    medium_name="LB_Rich_Medium",
    input_genes_nomenclature="name",
    output_genes_nomenclature="Base id",
    translation_filepath="../../Data/Palsson B.Subtilis Reconstruction/B_Subtilis Gene Translation.csv")
obj.make_genes_ko_growth_dict(filepath_to_save="../../Data/Palsson B.Subtilis Reconstruction/Genes KO Growth.json")
