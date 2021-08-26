import json


class Metabolite:
    def __init__(self, met_id, name, notes, local_index):
        self._id = met_id
        self._name = name
        self._notes = notes
        self._local_index = local_index

    def get_id(self):
        return self._id

    def get_name(self):
        return self._name

    def get_notes(self):
        return self._notes

    def get_local_index(self):
        return self._local_index


class Reaction:
    def __init__(self, rec_id, name, metabolites, notes, local_index):
        self._id = rec_id
        self._name = name
        self._metabolites = metabolites
        self._notes = notes
        self._local_index = local_index

    def get_id(self):
        return self._id

    def get_name(self):
        return self._name

    def get_metabolites(self):
        return self._metabolites

    def get_notes(self):
        return self._notes

    def get_local_index(self):
        return self._local_index


def get_fields(fields_str):
    repaired_str = repr(fields_str)[1:-1]
    print(repaired_str)
    fields = repaired_str.split("\\t")
    return fields


def make_metabolites_dict(metabolites_dicts_list):
    metabolites_dict = {}
    num_metabolites = len(metabolites_dicts_list)
    for index in range(num_metabolites):
        a_metabolite_dict = metabolites_dicts_list[index]
        the_metabolite_id = a_metabolite_dict['id']
        the_metabolite = Metabolite(met_id=the_metabolite_id,
                                    name=a_metabolite_dict['name'],
                                    notes=a_metabolite_dict['notes'],
                                    local_index=index)
        metabolites_dict[the_metabolite_id] = the_metabolite
    return metabolites_dict


def make_reactions_dict(reactions_dicts_list):
    reactions_dict = {}
    num_reactions = len(reactions_dicts_list)
    for index in range(num_reactions):
        a_reaction_dict = reactions_dicts_list[index]
        the_reaction_id = a_reaction_dict['id']
        the_reaction = Reaction(rec_id=the_reaction_id,
                                name=a_reaction_dict['name'],
                                metabolites=a_reaction_dict['metabolites'],
                                notes=a_reaction_dict['notes'],
                                local_index=index)
        reactions_dict[the_reaction_id] = the_reaction
    return reactions_dict


def save_lines(lines, save_path):
    with open(save_path, 'w') as f:
        for line in lines:
            f.write("%s\n" % line)


def save_sparse_stoichiometry_matrix(all_metabolites, all_reactions, save_path):
    num_metabolites = len(all_metabolites.keys())
    num_reactions = len(all_reactions.keys())
    first_row_str = str(num_metabolites) + ' ' + str(num_reactions)
    lines = [first_row_str]
    for rec_id, the_reaction in all_reactions.items():
        rec_index = the_reaction.get_local_index()
        the_rec_metabolites = the_reaction.get_metabolites()
        for met_id, the_coeff in the_rec_metabolites.items():
            met_index = all_metabolites[met_id].get_local_index()
            line_str = str(met_index) + ' ' + str(rec_index) + ' ' + str(the_coeff)
            lines.append(line_str)
    save_lines(lines, save_path)


def save_bigg_model(model_filename, save_path):
    # with open(metabolites_filename, 'r') as f:
    #     lines = f.readlines()
    #     headers = lines[0]
    #     """
    #     a list of 15721 metabolites
    #     """
    #     print(lines[-1])

    with open(model_filename, 'r') as f:
        model_data = json.load(f)
        """
        keys = ['metabolites', 'reactions', 'genes', 'id', 'compartments', 'version']
        model_data['metabolites'] = a list consists of 15638 dicts, each metabolite dict has keys:
            ['id', 'name', 'compartment', 'notes', 'annotation']
        model_data['reactions'] = a list consists of 28301 dicts, each reaction dict has keys:
            ['id', 'name', 'metabolites', 'lower_bound', 'upper_bound', 'gene_reaction_rule', 'notes', 'annotation']
            Note: all lower_bounds are 0.0 and all upper_bounds are 1000.0
        model_data['genes'] =  []  (Blank!)
        model_data['id'] = bigg_universal
        model_data['compartments'] = {}  (Blank!)
        model_data['version'] = 1
        """
        all_metabolites = make_metabolites_dict(model_data['metabolites'])
        all_reactions = make_reactions_dict(model_data['reactions'])
        save_sparse_stoichiometry_matrix(all_metabolites, all_reactions, save_path)


def main():
    # metabolites_filename = "../Data/BiGG Universal Model/bigg_models_metabolites.txt"
    # reactions_filename = "../Data/BiGG Universal Model/bigg_models_reactions.txt"
    model_filename = "../Data/BiGG Universal Model/universal_model.json"
    save_path = "../Data/BiGG Universal Model/S.txt"
    save_bigg_model(model_filename, save_path)

    # file_path = "../Data/Knock-Out Experiments - EcoCyc/LB enriched/LB_enriched_-_NG.txt"
    # with open(file_path) as f:
    #     lines = f.readlines()
    #     print(lines[0])
    #     data_fields = get_fields(lines[0][:-1])
    #
    # print(data_fields)


if __name__ == '__main__':
    main()
