"""
SourceUtilBoundsMaker
This code, gets two files:
    1. Sources Utilization Data: a .json file characterizing growth or non-growth in a medium with specific sources
    2. Medium Exchange Bounds: a .csv file characterizing the basic bounds for exchange reactions in a medium
and produces four .csv files as Growth Utilization Exchange Bounds:
    1. growth lower bounds
    2. growth upper bounds
    3. non-growth lower bounds
    4. non-growth upper bounds.

For more information, see the document.
"""

import json
import warnings
import pandas as pd


def find_exchange_by_name(metabolite_name, all_exchange_names):
    exchange_possible_name = 'EX_' + metabolite_name + '(e)'
    exceptions_dict = {'citr-L': 'EX_citr', 'abt-L': 'EX_abt(e)', 'f6p': 'Sink_f6p', 'glc-D': 'EX_glc(e)'}
    if metabolite_name in exceptions_dict.keys():
        exchange_name = exceptions_dict[metabolite_name]
    else:
        exchange_name = exchange_possible_name
    if exchange_name not in all_exchange_names:
        print(metabolite_name + ' not found in the list')
        raise Exception
    return exchange_name


def modify_uptakes(medium_bounds, sources_ids, exchanges_ids_list):
    """
    :param medium_bounds: A DataFrame denoting the base bounds for the medium
    :param sources_ids: List of are sources available (limited) in the medium
    :param exchanges_ids_list: List of all exchange reactions
    :return: medium_bounds: modified medium bounds as a DataFrame
    """
    for source_id in sources_ids:
        exchange_id = find_exchange_by_name(source_id, exchanges_ids_list)
        medium_bounds.loc[medium_bounds['ID'] == exchange_id, 'Lower Bound'] = -5
    return medium_bounds


class SourceUtilBoundsMaker:
    def __init__(self, sources_util_filepath, media_filepath_dict: dict):
        """
        :param sources_util_filepath: A string denoting the filepath for sources_util.json file.
                                      e.g. list items: {'sources_id': ['leu-L', 'nh4', 'pi', 'so4'],
                                                        'medium': 'minimal_media',
                                                        'growth': False,
                                                        'confidence_sc': 2.0}
        :param media_filepath_dict: A dictionary in the format of  {'media_name': "medium_bounds.csv"},
                                    denoting the filepath for each medium_bounds.csv
                                    Note: the first column (ID) for all media bounds file *should be identical*.
        """
        self.sources_util_filepath = sources_util_filepath
        self.sources_util = None
        self.media_filepath_dict = media_filepath_dict
        self.media_dict = {}
        self.exchanges_ids_list = None
        self.growth_lower_bounds = None
        self.growth_upper_bounds = None
        self.non_growth_lower_bounds = None
        self.non_growth_upper_bounds = None

    def load_sources_util(self):
        if self.sources_util is None:
            with open(self.sources_util_filepath, 'r') as json_file:
                self.sources_util = json.load(json_file)

    def add_medium(self, medium_name, medium_filepath):
        self.media_filepath_dict[medium_name] = medium_filepath

    def load_media_bounds(self):
        # if not self.media_dict:
        for medium_name, medium_filepath in self.media_filepath_dict.items():
            self.media_dict[medium_name] = pd.read_csv(medium_filepath)

    def initiate_growth_bounds(self, exchanges_ids_list):
        self.exchanges_ids_list = exchanges_ids_list
        self.growth_lower_bounds = pd.DataFrame(self.exchanges_ids_list, columns=['ID'])
        self.growth_upper_bounds = pd.DataFrame(self.exchanges_ids_list, columns=['ID'])

    def initiate_non_growth_bounds(self, exchanges_ids_list):
        self.exchanges_ids_list = exchanges_ids_list
        self.non_growth_lower_bounds = pd.DataFrame(self.exchanges_ids_list, columns=['ID'])
        self.non_growth_upper_bounds = pd.DataFrame(self.exchanges_ids_list, columns=['ID'])

    def make_growth_bounds(self):
        self.load_sources_util()
        self.load_media_bounds()
        counter = 0
        for source_data in self.sources_util:
            if source_data['growth']:  # == True
                if source_data['medium'] not in self.media_dict.keys():
                    warnings.warn("medium " + source_data['medium'] + " not specified")
                    continue
                medium_bounds = self.media_dict[source_data['medium']].copy()
                if counter == 0:
                    self.initiate_growth_bounds(exchanges_ids_list=medium_bounds['ID'].tolist())
                else:
                    if medium_bounds['ID'].tolist() != self.exchanges_ids_list:
                        print("Exchange reactions of media are not Identical")
                        raise Exception
                # Setting the lower bounds for the sources in this dict to -5:
                medium_bounds = modify_uptakes(medium_bounds=medium_bounds,
                                               sources_ids=source_data['sources_id'],
                                               exchanges_ids_list=self.exchanges_ids_list)
                confidence_str = '1'
                if 'confidence_sc' in source_data.keys():
                    confidence_str = str(source_data['confidence_sc'])
                new_column_name = str(counter + 1) + ', confidence: ' + confidence_str
                self.growth_lower_bounds['l' + new_column_name] = medium_bounds['Lower Bound'].tolist()
                self.growth_upper_bounds['u' + new_column_name] = medium_bounds['Upper Bound'].tolist()
                counter += 1

    def make_non_growth_bounds(self):
        self.load_sources_util()
        self.load_media_bounds()
        counter = 0
        for source_data in self.sources_util:
            if not source_data['growth']:  # == False
                if source_data['medium'] not in self.media_dict.keys():
                    warnings.warn("medium " + source_data['medium'] + " not specified")
                    continue
                medium_bounds = self.media_dict[source_data['medium']].copy()
                if counter == 0:
                    self.initiate_non_growth_bounds(exchanges_ids_list=medium_bounds['ID'].tolist())
                else:
                    if medium_bounds['ID'].tolist() != self.exchanges_ids_list:
                        print("Exchange reactions of media are not Identical")
                        raise Exception
                # Setting the lower bounds for the sources in this dict to -5:
                medium_bounds = modify_uptakes(medium_bounds=medium_bounds,
                                               sources_ids=source_data['sources_id'],
                                               exchanges_ids_list=self.exchanges_ids_list)
                confidence_str = '1'
                if 'confidence_sc' in source_data.keys():
                    confidence_str = str(source_data['confidence_sc'])
                new_column_name = str(counter + 1) + ', confidence: ' + confidence_str
                self.non_growth_lower_bounds['l' + new_column_name] = medium_bounds['Lower Bound'].tolist()
                self.non_growth_upper_bounds['u' + new_column_name] = medium_bounds['Upper Bound'].tolist()
                counter += 1

    def save_all_bounds(self, folder_to_save):
        self.growth_lower_bounds.to_csv(folder_to_save + 'g_lower_bounds.csv', index=False)
        self.growth_upper_bounds.to_csv(folder_to_save + 'g_upper_bounds.csv', index=False)
        self.non_growth_lower_bounds.to_csv(folder_to_save + 'ng_lower_bounds.csv', index=False)
        self.non_growth_upper_bounds.to_csv(folder_to_save + 'ng_upper_bounds.csv', index=False)


media_filepath_dict1 = {'minimal_media': "../Data/Palsson B.Subtilis Reconstruction/Biolog_Medium_Bounds.csv"}
obj = SourceUtilBoundsMaker(sources_util_filepath="../Data/Palsson B.Subtilis Reconstruction/Growth_Biolog.json",
                            media_filepath_dict=media_filepath_dict1)
obj.make_growth_bounds()
obj.make_non_growth_bounds()
obj.save_all_bounds(folder_to_save="../Data/Palsson B.Subtilis Reconstruction/Util Bounds/")
