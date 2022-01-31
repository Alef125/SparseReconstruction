"""
This script, provides as class to convert some prevalent format for the source utilization growth data, such as .csv,
into the standard .json format to be used in the further reconstruction steps.
"""

import pandas as pd
import json


def get_list_from_comma_separate(comma_separate_str):
    return comma_separate_str.split(", ")


class SourceUtilGrowthData:
    def __init__(self, csv_filepath, medium_name):
        self.csv_filepath = csv_filepath
        self.medium = medium_name

    def make_json_file(self, filepath_to_save):
        csv_file = pd.read_csv(self.csv_filepath)
        util_dicts = []
        for index, row in csv_file.iterrows():
            util_dicts.append(
                {'sources_id': get_list_from_comma_separate(row['IDs']),
                 'medium': self.medium,
                 'growth': row['Growth'],
                 'confidence_sc': row['Confidence Score']}
            )
        with open(filepath_to_save, 'w', encoding='utf-8') as f:
            json.dump(util_dicts, f, ensure_ascii=False, indent=4)


obj = SourceUtilGrowthData(csv_filepath="../Data/Palsson B.Subtilis Reconstruction/Growth_Biolog.csv",
                           medium_name="minimal_media")
obj.make_json_file(filepath_to_save="../Data/Palsson B.Subtilis Reconstruction/Growth_Biolog.json")
