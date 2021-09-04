import json
import GeneKnockOutParser
# import GeneAssociationMaker


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
    def __init__(self, rec_id, name, metabolites, notes, lower_bound, upper_bound, bigg_id, local_index):
        self._id = rec_id
        self._name = name
        self._metabolites = metabolites
        self._notes = notes
        self._lower_bound = lower_bound
        self._upper_bound = upper_bound
        self._local_index = local_index
        self._bigg_id = bigg_id

    def get_id(self):
        return self._id

    def get_name(self):
        return self._name

    def get_metabolites(self):
        return self._metabolites

    def get_notes(self):
        return self._notes

    def get_lower_bound(self):
        return self._lower_bound

    def get_upper_bound(self):
        return self._upper_bound

    def get_bigg_id(self):
        return self._bigg_id

    def get_local_index(self):
        return self._local_index


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


def get_bigg_ids(all_reactions):
    bigg_ids = []
    for rec_id, the_reaction in all_reactions.items():
        bigg_ids.append(the_reaction.get_bigg_id())
    return bigg_ids


def make_reactions_dict(reactions_dicts_list):
    reactions_dict = {}
    num_reactions = len(reactions_dicts_list)
    for index in range(num_reactions):
        a_reaction_dict = reactions_dicts_list[index]
        the_reaction_id = a_reaction_dict['id']
        original_bigg_id = a_reaction_dict['notes']['original_bigg_ids'][0]
        the_reaction = Reaction(rec_id=the_reaction_id,
                                name=a_reaction_dict['name'],
                                metabolites=a_reaction_dict['metabolites'],
                                notes=a_reaction_dict['notes'],
                                lower_bound=a_reaction_dict['lower_bound'],
                                upper_bound=a_reaction_dict['upper_bound'],
                                bigg_id=original_bigg_id,
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
            line_str = str(met_index+1) + ' ' + str(rec_index+1) + ' ' + str(the_coeff)  # +1 for python indexing from 0
            lines.append(line_str)
    save_lines(lines, save_path + "/S.txt")


def save_lower_and_upper_bounds(all_reactions, bounds_path):
    num_reactions = len(all_reactions.keys())
    lbs = [str(num_reactions)]
    ubs = [str(num_reactions)]
    for rec_id, the_reaction in all_reactions.items():
        rec_lb = the_reaction.get_lower_bound()
        rec_ub = the_reaction.get_upper_bound()
        lbs.append(str(rec_lb))
        ubs.append(str(rec_ub))
    save_lines(lbs, bounds_path + "/l.txt")
    save_lines(ubs, bounds_path + "/u.txt")


def save_bigg_model(model_filename, save_path, do_save):
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
            Note: It seems that all lower_bounds are 0.0 and all upper_bounds are 1000.0
            'notes': {'original_bigg_ids': ['HMR_3915']}
        model_data['genes'] =  []  (Blank!)
        model_data['id'] = bigg_universal
        model_data['compartments'] = {}  (Blank!)
        model_data['version'] = 1
        """
        all_metabolites = make_metabolites_dict(model_data['metabolites'])
        all_reactions = make_reactions_dict(model_data['reactions'])
        if do_save:
            save_sparse_stoichiometry_matrix(all_metabolites, all_reactions, save_path)
            save_lower_and_upper_bounds(all_reactions, save_path)
        return get_bigg_ids(all_reactions)


def main():
    # metabolites_filename = "../Data/BiGG Universal Model/bigg_models_metabolites.txt"
    # reactions_filename = "../Data/BiGG Universal Model/bigg_models_reactions.txt"
    model_filename = "../Data/BiGG Universal Model/universal_model.json"
    universal_model_path = "../Data/BiGG Universal Model"
    all_reactions_bigg_ids = save_bigg_model(model_filename, universal_model_path, do_save=False)
    print(all_reactions_bigg_ids)

    # ecoli_model_path = "../Data/Escherichia coli str. K-12 substr. MG1655/iAF1260b.json"
    # GeneAssociationMaker.make_associations(ecoli_model_path)

    associations_path = "../Data/Escherichia coli str. K-12 substr. MG1655/ecoli_genes_associations.json"
    knock_outs_path = "../Data/Knock-Out Experiments - EcoCyc/LB enriched/LB_enriched_-_G.txt"
    bounds_save_path = "../Data/Knock-Out Experiments - EcoCyc/LB enriched"
    GeneKnockOutParser.save_knock_out_bounds(all_reactions_ids_list=all_reactions_bigg_ids,
                                             gene_associations_path=associations_path,
                                             knock_outs_path=knock_outs_path,
                                             universal_model_path=universal_model_path,
                                             save_path=bounds_save_path)


if __name__ == '__main__':
    main()
