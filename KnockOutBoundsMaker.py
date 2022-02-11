"""
KnockOutBoundsMaker
This code, gets three files:
    1. Reactions KO Essentiality Data: a .json file characterizing growth or non-growth corresponding to each
     knocking-out experiment in a medium
    2. Medium Exchange Bounds: a .csv file characterizing the basic bounds for exchange reactions in a medium
    3. Internal Reactions Bounds: a .csv file characterizing the default bounds for internal reactions
and produces four .csv files as Growth Utilization Exchange Bounds:
    1. growth lower bounds
    2. growth upper bounds
    3. non-growth lower bounds
    4. non-growth upper bounds.

Note: If not "all" the reactions from an KO experiment are found in the reactions list,
      we ignore that experiment as partial shutting downs can be problematic.
For more information, see the document.
"""

import json
import warnings
import pandas as pd
import os


def modify_ko_bounds(total_bounds, ko_rxns_ids):
    """
    :param total_bounds: A DataFrame denoting the base bounds for all the reactions in the medium
    :param ko_rxns_ids: List of reactions which are knocked-out in this specific experiment
    :return: all_reactions_found_flag: True if all knocked-out reactions were completely shut down,
             to prevent misleading KO data generation
    """
    all_reactions_found_flag = True
    all_reactions_ids = total_bounds['ID'].tolist()
    for ko_rxn_id in ko_rxns_ids:
        if ko_rxn_id in all_reactions_ids:
            total_bounds.loc[total_bounds['ID'] == ko_rxn_id, 'Lower Bound'] = 0
            total_bounds.loc[total_bounds['ID'] == ko_rxn_id, 'Upper Bound'] = 0
        else:
            warning_text = "The reaction with ID " + ko_rxn_id + " does not exist in the list"
            warnings.warn(warning_text)
            all_reactions_found_flag = False
    return all_reactions_found_flag, total_bounds


class KnockOutBoundsMaker:
    def __init__(self, reactions_ko_filepath: str, media_filepath_dict: dict, internal_rxns_filepath: str):
        """
        :param reactions_ko_filepath: A string denoting the filepath for reactions_ko.json file.
                                      e.g. list items: {'ko_rxns_ids': ["TRPS1", "TRPS3", "TRPS2"],
                                                        'medium': 'LB_Rich_Medium',
                                                        'growth': true}
        :param media_filepath_dict: A dictionary in the format of  {'media_name': "medium_bounds.csv"},
                                    denoting the filepath for each medium_bounds.csv
                                    Note: the first column (ID) for all media bounds file *should be identical*.
        :param internal_rxns_filepath: A string denoting the filepath for internal_rxns_bounds.csv file.
        """
        self.reactions_ko_filepath = reactions_ko_filepath
        self.reactions_ko_list = None
        self.load_reactions_ko_list()
        # ############################################
        self.media_filepath_dict = media_filepath_dict
        self.media_dict = {}
        # ##################################################
        self.internal_rxns_filepath = internal_rxns_filepath
        self.internal_rxns_df = None
        self.load_internal_rxns_bounds()
        # ##############################
        self.all_exchange_ids = None
        self.growth_lower_bounds = None
        self.growth_upper_bounds = None
        self.non_growth_lower_bounds = None
        self.non_growth_upper_bounds = None

    def load_reactions_ko_list(self):
        """
        This method, loads the reactions_ko_list.json from reactions_ko_filepath
        :return: -
        """
        with open(self.reactions_ko_filepath, 'r') as json_file:
            self.reactions_ko_list = json.load(json_file)

    def add_medium(self, medium_name: str, medium_filepath: str):
        """
        :param medium_name: The name of the medium to be added
        :param medium_filepath: The path to the medium bounds .csv file
        :return: -. Adds this new medium to the self.media_filepath_dict
        """
        self.media_filepath_dict[medium_name] = medium_filepath

    def load_media_bounds(self):
        """
        This method, loads the media_bounds .json files based on the paths in self.media_filepath_dict
        :return:
        """
        # if not self.media_dict:
        for medium_name, medium_filepath in self.media_filepath_dict.items():
            self.media_dict[medium_name] = pd.read_csv(medium_filepath)

    def load_internal_rxns_bounds(self):
        """
        This method, loads the internal reactions bounds .csv file from self.internal_rxns_filepath
        :return: -
        """
        self.internal_rxns_df = pd.read_csv(self.internal_rxns_filepath)

    def initiate_growth_bounds(self, internal_ids_list: list, exchanges_ids_list: list):
        """
        :param internal_ids_list: The list of all internal reactions ids
        :param exchanges_ids_list: The list of all exchange reactions ids (same in all the media)
        :return: -. Initializes the growth bounds as a DataFrame
        """
        self.all_exchange_ids = exchanges_ids_list
        all_reactions_ids = internal_ids_list + exchanges_ids_list
        self.growth_lower_bounds = pd.DataFrame(all_reactions_ids, columns=['ID'])
        self.growth_upper_bounds = pd.DataFrame(all_reactions_ids, columns=['ID'])

    def initiate_non_growth_bounds(self, internal_ids_list: list, exchanges_ids_list: list):
        """
        :param internal_ids_list: The list of all internal reactions ids
        :param exchanges_ids_list: The list of all exchange reactions ids (same in all the media)
        :return: -. Initializes the non-growth bounds as a DataFrame
        """
        self.all_exchange_ids = exchanges_ids_list
        all_reactions_ids = internal_ids_list + exchanges_ids_list
        self.non_growth_lower_bounds = pd.DataFrame(all_reactions_ids, columns=['ID'])
        self.non_growth_upper_bounds = pd.DataFrame(all_reactions_ids, columns=['ID'])

    def make_growth_bounds(self):
        """
        This method, fills self.growth_lower_bounds and self.growth_upper_bounds based on
        self.reactions_ko_list
        :return: -
        """
        self.load_media_bounds()
        counter = 0
        for ko_data in self.reactions_ko_list:
            if ko_data['growth']:  # == true
                if ko_data['medium'] not in self.media_dict.keys():
                    warnings.warn("medium " + ko_data['medium'] + " not specified")
                    continue
                medium_bounds = self.media_dict[ko_data['medium']].copy()
                if counter == 0:
                    self.initiate_growth_bounds(internal_ids_list=self.internal_rxns_df['ID'].tolist(),
                                                exchanges_ids_list=medium_bounds['ID'].tolist())
                else:
                    if medium_bounds['ID'].tolist() != self.all_exchange_ids:
                        print("Exchange reactions of media are not Identical")
                        raise Exception
                # Setting the bounds for the knocked-out reaction in this dict to 0:
                total_bounds_df = pd.concat([self.internal_rxns_df, medium_bounds],
                                            ignore_index=True)
                is_valid, modified_bounds = modify_ko_bounds(total_bounds=total_bounds_df,
                                                             ko_rxns_ids=ko_data['ko_rxns_ids'])
                if is_valid:
                    new_column_name = str(counter + 1)
                    self.growth_lower_bounds['l' + new_column_name] = modified_bounds['Lower Bound'].tolist()
                    self.growth_upper_bounds['u' + new_column_name] = modified_bounds['Upper Bound'].tolist()
                    counter += 1

    def make_non_growth_bounds(self):
        """
        This method, fills self.non_growth_lower_bounds and self.non_growth_upper_bounds based on
        self.reactions_ko_list
        :return: -
        """
        self.load_media_bounds()
        counter = 0
        for ko_data in self.reactions_ko_list:
            if not ko_data['growth']:  # == false
                if ko_data['medium'] not in self.media_dict.keys():
                    warnings.warn("medium " + ko_data['medium'] + " not specified")
                    continue
                medium_bounds = self.media_dict[ko_data['medium']].copy()
                if counter == 0:
                    self.initiate_non_growth_bounds(internal_ids_list=self.internal_rxns_df['ID'].tolist(),
                                                    exchanges_ids_list=medium_bounds['ID'].tolist())
                else:
                    if medium_bounds['ID'].tolist() != self.all_exchange_ids:
                        print("Exchange reactions of media are not Identical")
                        raise Exception
                # Setting the bounds for the knocked-out reaction in this dict to 0:
                total_bounds_df = pd.concat([self.internal_rxns_df, medium_bounds],
                                            ignore_index=True)
                is_valid, modified_bounds = modify_ko_bounds(total_bounds=total_bounds_df,
                                                             ko_rxns_ids=ko_data['ko_rxns_ids'])
                if is_valid:
                    new_column_name = str(counter + 1)
                    self.non_growth_lower_bounds['l' + new_column_name] = modified_bounds['Lower Bound'].tolist()
                    self.non_growth_upper_bounds['u' + new_column_name] = modified_bounds['Upper Bound'].tolist()
                    counter += 1

    def save_all_bounds(self, folder_to_save: str):
        """
        :param folder_to_save: Folder path to save all four growth and non-growth bound
        :return: -
        """
        if not os.path.exists(folder_to_save):
            os.makedirs(folder_to_save)
        self.growth_lower_bounds.to_csv(folder_to_save + 'g_lower_bounds.csv', index=False)
        self.growth_upper_bounds.to_csv(folder_to_save + 'g_upper_bounds.csv', index=False)
        self.non_growth_lower_bounds.to_csv(folder_to_save + 'ng_lower_bounds.csv', index=False)
        self.non_growth_upper_bounds.to_csv(folder_to_save + 'ng_upper_bounds.csv', index=False)


media_filepath_dict_ = {'LB_Rich_Medium': "../Data/Palsson B.Subtilis Reconstruction/LB_Medium_Bounds.csv"}
obj = KnockOutBoundsMaker(reactions_ko_filepath="../Data/Palsson B.Subtilis Reconstruction/Reactions KO Growth.json",
                          media_filepath_dict=media_filepath_dict_,
                          internal_rxns_filepath="../Data/Palsson B.Subtilis Reconstruction/Internal_Rxns_Bounds.csv")
obj.make_growth_bounds()
obj.make_non_growth_bounds()
obj.save_all_bounds(folder_to_save="../Data/Palsson B.Subtilis Reconstruction/KO Bounds/")
