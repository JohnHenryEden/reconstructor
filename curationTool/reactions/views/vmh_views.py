
import json
from django.views.decorators.csrf import csrf_exempt
from reactions.models import User, Reaction, MetabolitesAddedVMH, ReactionsAddedVMH, Subsystem
from reactions.reaction_info import construct_vmh_formula
from reactions.utils.search_vmh import search_metabolites_vmh, is_name_in_vmh
from reactions.utils.to_smiles import any_to_smiles
from reactions.utils.utils import capitalize_first_letter
from reactions.utils.gen_vmh_abbrs import gen_metabolite_abbr
from reactions.utils.add_to_vmh_utils import check_met_names_abbrs_vmh, check_names_abbrs_vmh, gather_reaction_details, rxn_prepare_json_paths_and_variables, met_prepare_json_paths_and_variables, add_reaction_matlab, add_metabolites_matlab, smiles_to_inchikeys, smiles_to_charged_formula, get_nonfound_metabolites
from reactions.utils.MatlabSessionManager import MatlabSessionManager
from reactions.utils.search_vmh import check_reaction_vmh, get_from_vmh
import requests
import os
from django.http import JsonResponse
import logging


def get_vmh_subsystems():
    BASE_URL = 'https://www.vmh.life/'
    endpoint = f"{BASE_URL}_api/subsystems/"
    subsystems = []

    # Fetch subsystems from VMH
    while True:
        response = requests.get(endpoint, verify=False)
        data = response.json()['results']
        subsystems.extend([subsystem['name'] for subsystem in data])
        endpoint = response.json().get('next')
        if not endpoint:
            break

    return subsystems


def get_subsystems(request):
    try:
        subsystems = get_vmh_subsystems()

        # Fetch subsystems from the database
        stored_subsystems = Subsystem.objects.values_list('name', flat=True)

        # Merge the two lists
        combined_subsystems = set(subsystems).union(set(stored_subsystems))
        combined_subsystems = list(combined_subsystems)

        return JsonResponse({'subsystem_list': combined_subsystems})

    except Exception as e:
        return JsonResponse({'error': True, 'message': str(e)}, status=500)


@csrf_exempt
def update_subsystems(request):
    if request.method == 'POST':
        try:
            data = json.loads(request.body)
            subsystems = data.get('subsystems', [])
            # Add new subsystems to the database
            for subsystem in subsystems:
                Subsystem.objects.get_or_create(name=subsystem)

            return JsonResponse({'message': 'Subsystems updated successfully'})
        except Exception as e:
            return JsonResponse({'error': True, 'message': str(e)}, status=500)
    return JsonResponse(
        {'error': True, 'message': 'Invalid request method'}, status=400)


@csrf_exempt
def prepare_add_to_vmh(request):
    if request.method != 'POST':
        return JsonResponse({'status': 'error',
                             'message': 'Invalid request method. Use POST instead.'},
                            status=400)

    try:
        request_data = json.loads(request.body)
        reactionIds = request_data['reactionIds']
    except json.JSONDecodeError:
        return JsonResponse(
            {'status': 'error', 'message': 'Invalid JSON data.'}, status=400)
    except KeyError:
        return JsonResponse({'status': 'error',
                             'message': 'Missing reactionIds in the request data.'},
                            status=400)
    except Exception as e:
        return JsonResponse(
            {'status': 'error', 'message': 'An unexpected error occurred.'}, status=500)

    try:
        reaction_objs = [
            Reaction.objects.get(
                pk=int(reactionId)) for reactionId in reactionIds]

        in_vmh = [
            reaction.vmh_found and not reaction.vmh_found_similar for reaction in reaction_objs]
        if True in in_vmh:
            reaction_objs_in_vmh = [
                reaction for reaction, found in zip(
                    reaction_objs, in_vmh) if found]
            names = [reaction.short_name for reaction in reaction_objs_in_vmh]
            return JsonResponse(
                {'status': 'error', 'message': f'The following reactions are already in VMH: {", ".join(names)}'})

        subs_in_vmh = [json.loads(reaction.subs_found) if reaction.subs_found else [
        ] for reaction in reaction_objs]
        prods_in_vmh = [json.loads(reaction.prod_found) if reaction.prod_found else [
        ] for reaction in reaction_objs]

        matlab_session = MatlabSessionManager()
        subs_abbr = [
            [
                gen_metabolite_abbr(
                    sub,
                    sub_type,
                    sub_name,
                    search_metabolites_vmh,
                    matlab_session) for sub,
                sub_type,
                sub_name in zip(
                    json.loads(
                        reaction.substrates),
                    json.loads(
                        reaction.substrates_types),
                    json.loads(
                        reaction.substrates_names))] for reaction in reaction_objs]
        prods_abbr = [
            [
                gen_metabolite_abbr(
                    prod,
                    prod_type,
                    prod_name,
                    search_metabolites_vmh,
                    matlab_session) for prod,
                prod_type,
                prod_name in zip(
                    json.loads(
                        reaction.products),
                    json.loads(
                        reaction.products_types),
                    json.loads(
                        reaction.products_names))] for reaction in reaction_objs]

        subs_need_new_names = [[] for _ in reactionIds]
        prods_need_new_names = [[] for _ in reactionIds]

        for idx, in_vmh_list in enumerate(subs_in_vmh):
            for j, sub_in_vmh in enumerate(in_vmh_list):
                if not sub_in_vmh:
                    need_new_name = is_name_in_vmh(json.loads(
                        reaction_objs[idx].substrates_names)[j])
                    subs_need_new_names[idx].append(need_new_name)
                else:
                    subs_need_new_names[idx].append(False)
        for idx, in_vmh_list in enumerate(prods_in_vmh):
            for j, prod_in_vmh in enumerate(in_vmh_list):
                if prod_in_vmh:
                    prods_need_new_names[idx].append(False)
                else:
                    need_new_name = is_name_in_vmh(json.loads(
                        reaction_objs[idx].products_names)[j])
                    prods_need_new_names[idx].append(need_new_name)

        reaction_abbrs = ['' for _ in reactionIds]

        return JsonResponse({
            'status': 'success',
            'subs_in_vmh': subs_in_vmh,
            'prods_in_vmh': prods_in_vmh,
            'subs_abbr': subs_abbr,
            'prods_abbr': prods_abbr,
            'subs_need_new_names': subs_need_new_names,
            'prods_need_new_names': prods_need_new_names,
            'reaction_abbrs': reaction_abbrs,
        })
    except Exception as e:
        return JsonResponse({'status': 'error',
                             'message': 'An error occurred while processing reactions.'},
                            status=500)


def bypass_search_func(metabolites, types, *args, **kwargs):
    # Always return False for found and None for abbreviation
    return [False], [None]


@csrf_exempt  # Remove this if you want to deal with CSRF tokens later
def create_formula_abbr(request):
    if request.method == 'POST':
        # Parse JSON data from the request body
        try:
            data = json.loads(request.body)
            metabolite = data.get('metabolite')
            mtype = data.get('mtype')
            metabolite_name = data.get('metabolite_name')
        except json.JSONDecodeError:
            return JsonResponse({'error': 'Invalid JSON'}, status=400)

        # Ensure all required fields are present
        if not all([metabolite, mtype, metabolite_name]):
            return JsonResponse({'error': 'Missing data'}, status=400)

        # Initialize Matlab session manager
        matlab_session = MatlabSessionManager()

        # Use the bypass function to force the else clause
        abbr = gen_metabolite_abbr(
            metabolite,
            mtype,
            metabolite_name,
            bypass_search_func,
            matlab_session)

        # Return the abbreviation as JSON
        return JsonResponse({'abbr': abbr})

    return JsonResponse({'error': 'Invalid request method'}, status=400)


def add_to_vmh(request):
    """
    Main function to handle the request for adding reactions to VMH.
    """
    req_body = json.loads(request.body)

    userID = req_body.get('userID')
    user = User.objects.get(pk=userID)
    if not user:
        return JsonResponse(
            {'status': 'error', 'message': 'Invalid user key.'}, status=404)
    if not user.cred_add_to_vmh:
        return JsonResponse({'status': 'error',
                             'message': 'User does not have permission to add to VMH.'},
                            status=403)
    user_name = user.name
    user_full_name = user.name
    reactions = req_body.get('reactions')
    reaction_ids = []
    not_enough_info, no_comments, not_balanced = [], [], []
    met_added_info = {}
    names_list = [reaction['short_name'] for reaction in reactions]
    abbr_list = [reaction['abbreviation'] for reaction in reactions]

    # Check for abbreviations, names, and confidence scores
    missing_names = [reaction['short_name'] == '' for reaction in reactions]
    missing_abbrs = [reaction['abbreviation'] == '' for reaction in reactions]
    missing_conf_scores = [reaction['confidence_score']
                           == '" "' for reaction in reactions]

    if True in missing_names:
        return JsonResponse(
            {'status': 'error', 'message': 'Please enter a description for reaction'})

    if True in missing_abbrs:
        return JsonResponse({'status': 'error',
                             'message': 'Please enter an abbreviation'})

    if True in missing_conf_scores:
        return JsonResponse(
            {'status': 'error', 'message': 'Please enter a confidence score for all reactions'})

    for name in names_list:
        if names_list.count(name) > 1:
            return JsonResponse(
                {'status': 'error', 'message': f'Reaction with name `{name}` is repeated in the list.'})
    for abbr in abbr_list:
        if abbr_list.count(abbr) > 1:
            return JsonResponse(
                {'status': 'error', 'message': f'Reaction with abbreviation `{abbr}` is repeated in the list.'})

    name_in_vmh, abbr_in_vmh = check_names_abbrs_vmh(
        [(reaction['short_name'], reaction['abbreviation']) for reaction in reactions])
    if True in list(name_in_vmh.values()):
        name_in_vmh_reactions = [
            reaction for reaction, in_vmh in zip(
                reactions, name_in_vmh.values()) if in_vmh]
        return JsonResponse(
            {
                'status': 'error',
                'message': f'The following reaction descriptions are already in VMH: {", ".join([reaction["short_name"] for reaction in name_in_vmh_reactions])}'})
    if True in list(abbr_in_vmh.values()):
        abbr_in_vmh_reactions = [
            reaction for reaction, in_vmh in zip(
                reactions, abbr_in_vmh.values()) if in_vmh]
        return JsonResponse(
            {
                'status': 'error',
                'message': f'The following reaction abbreviations are already in VMH: {", ".join([reaction["abbreviation"] for reaction in abbr_in_vmh_reactions])}'})

    reactions_new_subsInfo = [
        json.loads(
            reaction['substrates_info']) for reaction in reactions]
    reactions_new_prodsInfo = [json.loads(
        reaction['products_info']) for reaction in reactions]
    reactions_subs_found = [
        json.loads(
            Reaction.objects.get(
                pk=reaction['pk']).subs_found) for reaction in reactions]
    reactions_prods_found = [
        json.loads(
            Reaction.objects.get(
                pk=reaction['pk']).prod_found) for reaction in reactions]
    subs_names_vmh, subs_abbr_vmh, prods_names_vmh, prods_abbr_vmh = check_met_names_abbrs_vmh(
        reactions_new_subsInfo, reactions_new_prodsInfo, reactions_subs_found, reactions_prods_found)
    if True in list(subs_names_vmh.values()):
        subs_names_in_vmh = [
            sub for sub in subs_names_vmh.keys() if subs_names_vmh[sub]]
        return JsonResponse(
            {
                'status': 'error',
                'message': f'The following substrates have metabolite names that are already in VMH: {", ".join(subs_names_in_vmh)}'})
    if True in list(subs_abbr_vmh.values()):
        subs_abbrs_in_vmh = [
            sub for sub in subs_abbr_vmh.keys() if subs_abbr_vmh[sub]]
        return JsonResponse(
            {
                'status': 'error',
                'message': f'The following substrates have metabolite abbreviations that are already in VMH: {", ".join(subs_abbrs_in_vmh)}'})
    if True in list(prods_names_vmh.values()):
        prods_names_in_vmh = [
            prod for prod in prods_names_vmh.keys() if prods_names_vmh[prod]]
        return JsonResponse(
            {
                'status': 'error',
                'message': f'The following products have metabolite names that are already in VMH: {", ".join(prods_names_in_vmh)}'})
    if True in list(prods_abbr_vmh.values()):
        prods_abbrs_in_vmh = [
            prod for prod in prods_abbr_vmh.keys() if prods_abbr_vmh[prod]]
        return JsonResponse(
            {
                'status': 'error',
                'message': f'The following products have metabolite abbreviations that are already in VMH: {", ".join(prods_abbrs_in_vmh)}'})
    subs_abbr = []
    prods_abbr = []
    for reaction in reactions:

        obj = Reaction.objects.get(pk=reaction['pk'])
        obj.short_name = reaction['short_name']
        reaction_ids.append(obj.id)
        # Update substrate names and abbreviations
        subs_info = json.loads(reaction['substrates_info'])
        new_subs_names = [
            capitalize_first_letter(
                sub['name']) for sub in subs_info]
        new_subs_abbrs = [sub['abbreviation'] for sub in subs_info]
        subs_abbr.append(new_subs_abbrs)
        obj.substrates_names = json.dumps(new_subs_names)
        # Update product names and abbreviations
        prods_info = json.loads(reaction['products_info'])
        new_prods_names = [
            capitalize_first_letter(
                prod['name']) for prod in prods_info]
        new_prods_abbrs = [prod['abbreviation'] for prod in prods_info]
        prods_abbr.append(new_prods_abbrs)
        obj.products_names = json.dumps(new_prods_names)
        # Update references, external links, and comments
        references, ext_links, comments = json.loads(
            reaction['references']), json.loads(
            reaction['ext_links']), json.loads(
            reaction['comments'])
        new_references, new_ext_links, new_comments = [], [], []
        for ref in references:
            ref['user_name'] = user_name
            if not ('PMID' in ref['info'] or 'DOI' in ref['info']):
                ref['info'] = f"{ref['ref_type']}:{ref['info']}"
            new_references.append(ref)
        for link in ext_links:
            link['user_name'] = user_name
            new_ext_links.append(link)
        for comment in comments:
            comment['user_name'] = user_name
            new_comments.append(comment)
        if f"Created and Added to VMH by: {user_full_name}" not in list(
                map(lambda x: x['info'], new_comments)):
            new_comments.append(
                {'info': f"Created and Added to VMH by: {user_full_name}", 'user_name': user_name})
        not_enough_info.append(
            len(new_references) < 1 and len(new_ext_links) < 1)
        no_comments.append(len(new_comments) < 2)
        obj.references = new_references if new_references else None
        obj.ext_links = new_ext_links if new_ext_links else None
        obj.comments = new_comments if new_comments else None
        obj.confidence_score = reaction.get(
            'confidence_score', 0)  # Add confidence score
        if not (
            json.loads(
                obj.balanced_charge)[0] and json.loads(
                obj.balanced_count)[0]):
            not_balanced.append(True)
        else:
            not_balanced.append(False)
        # Save the updated reaction object
        obj.save()
    reaction_objs = [
        Reaction.objects.get(
            pk=reaction_id) for reaction_id in reaction_ids]

    if True in not_enough_info:
        not_enough_info_reactions = [
            reaction for reaction, not_enough in zip(
                reaction_objs, not_enough_info) if not_enough]
        return JsonResponse(
            {
                'status': 'error',
                'message': f'The following reactions do not have at least one reference or external link: {", ".join([reaction.short_name for reaction in not_enough_info_reactions])}'})
    if True in no_comments:
        no_comments_reactions = [
            reaction for reaction,
            no_comment in zip(
                reaction_objs,
                no_comments) if no_comment]
        return JsonResponse(
            {
                'status': 'error',
                'message': f'The following reactions do not have at least one comment: {", ".join([reaction.short_name for reaction in no_comments_reactions])}'})
    if True in not_balanced:
        not_balanced_reactions = [
            reaction for reaction,
            not_balanced in zip(
                reaction_objs,
                not_balanced) if not_balanced]
        return JsonResponse(
            {
                'status': 'error',
                'message': f'The following reactions are not balanced: {", ".join([reaction.short_name for reaction in not_balanced_reactions])}'})
    all_vmh = False if any([any(element != 'VMH' for element in json.loads(reaction.substrates_types)) or any(
        element != 'VMH' for element in json.loads(reaction.products_types)) for reaction in reaction_objs]) else True
    reaction_identifiers, reaction_names = [
        reaction['abbreviation'] for reaction in reactions], [
        reaction.short_name for reaction in reaction_objs]
    matlab_session = MatlabSessionManager()
    if not all_vmh:
        unique_abbrs, unique_mols, unique_types, unique_names = get_nonfound_metabolites(
            reaction_objs, subs_abbr, prods_abbr, search_func=search_metabolites_vmh)
        if unique_abbrs and unique_mols and unique_types:
            abbrs = unique_abbrs
            smiles, errors = any_to_smiles(
                unique_mols, unique_types, None, side=None)
            smiles = ['' if v is not None else smiles[i]
                      for i, v in enumerate(errors)]
            inchikeys = smiles_to_inchikeys(smiles)
            names = unique_names
            formulas, charges = smiles_to_charged_formula(smiles)
            json_paths = met_prepare_json_paths_and_variables(
                abbrs, names, formulas, charges, inchikeys, smiles)
            matlab_result = add_metabolites_matlab(json_paths, matlab_session)
            for path in json_paths:
                os.remove(path)
            if matlab_result['status'] == 'success':
                met_ids = matlab_result['met_ids']
                met_added_info = {
                    abbr: [
                        met_id,
                        formula,
                        inchikey] for abbr,
                    met_id,
                    formula,
                    inchikey in zip(
                        abbrs,
                        met_ids,
                        formulas,
                        inchikeys)}
                for abbr, info in met_added_info.items():
                    MetabolitesAddedVMH.objects.create(
                        user=user,
                        user_name=user_name,
                        metabolite_id=info[0],
                        metabolite_formula=info[1],
                        metabolite_abbr=abbr,
                    )
            else:
                matlab_session.quit()
                return JsonResponse(
                    {'status': 'error', 'message': matlab_result['message']})
    reaction_formulas = [
        construct_vmh_formula(
            reaction_objs[idx],
            subs_abbr[idx],
            prods_abbr[idx]) for idx in range(
            len(reaction_objs))]
    reaction_directions, reaction_subsystems, reaction_references, reaction_external_links, reaction_gene_info, reaction_comments, reaction_confidence_scores = gather_reaction_details(
        reaction_objs)
    json_paths = rxn_prepare_json_paths_and_variables(
        reaction_identifiers,
        reaction_names,
        reaction_formulas,
        reaction_directions,
        reaction_subsystems,
        reaction_references,
        reaction_external_links,
        reaction_gene_info,
        reaction_comments,
        reaction_confidence_scores)
    matlab_result = add_reaction_matlab(json_paths, matlab_session)
    # Cleanup temporary JSON files
    for path in json_paths:
        os.remove(path)

    if matlab_result['status'] == 'success':
        rxn_added_info = {
            abbr: [
                rxn_id,
                reaction_formula] for abbr,
            rxn_id,
            reaction_formula in zip(
                reaction_identifiers,
                matlab_result['rxn_ids'],
                reaction_formulas)}
        for abbr, info in rxn_added_info.items():
            ReactionsAddedVMH.objects.create(
                user=user,
                user_name=user_name,
                reaction_id=info[0],
                reaction_formula=info[1],
                reaction_abbr=abbr,
            )
        matlab_session.quit()
        return JsonResponse({'status': 'success',
                             'rxn_added_info': rxn_added_info,
                             'met_added_info': met_added_info})
    else:
        matlab_session.quit()
        return JsonResponse(
            {'status': 'error', 'message': matlab_result['message']})
