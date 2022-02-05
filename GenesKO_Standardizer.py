"""
GenesKO_Standardizer.py
This code, converts common formats for Knock-Out Essentiality Data, such as .csv, into the standard .json file
used in further processes
"""

import json
import pandas as pd


class GenesKOStandardizer:
    def __init__(self, genes_ko_growth_filepath, medium_name):
        self.genes_ko_growth_filepath = genes_ko_growth_filepath
        self.genes_ko_growth_file = None
        self.read_genes_ko_growth_file()
        self.medium = medium_name

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

    def make_genes_ko_growth_dict(self, filepath_to_save):
        """
        :param filepath_to_save: The path to save the .json file
        :return: Nothing, Saves the .json file
        """
        genes_ko_dicts = []
        for index, row in self.genes_ko_growth_file.iterrows():
            genes_ko_dicts.append(
                {'ko_gene_id': row['Gene'],
                 'medium': self.medium,
                 'growth': row['Growth']}
            )
        with open(filepath_to_save, 'w', encoding='utf-8') as f:
            json.dump(genes_ko_dicts, f, ensure_ascii=False, indent=4)


obj = GenesKOStandardizer(genes_ko_growth_filepath="../Data/Palsson B.Subtilis Reconstruction/Genes KO Growth Data.csv",
                          medium_name="LB_Rich_Medium")
obj.make_genes_ko_growth_dict(filepath_to_save="../Data/Palsson B.Subtilis Reconstruction/Genes KO Growth.json")
