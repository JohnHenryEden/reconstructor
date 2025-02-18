import json
import pandas as pd
from reactions.utils.search_vmh import search_vmh
from reactions.utils.to_mol import any_to_mol
from reactions.utils.utils import capitalize_first_letter
from reactions.utils.get_from_rhea import get_from_rhea
from reactions.models import SavedMetabolite, User, Reaction
import requests
from collections import defaultdict
from django.views.decorators.csrf import csrf_exempt
from django.forms.models import model_to_dict
from django.http import JsonResponse
from django.db.models import Q
import rdkit.Chem as Chem

def verify_metabolite(request):
    """
    Verifies if a metabolite input from the user is in VMH.
    """
    main_input = request.POST.get('metabolite')
    input_type = request.POST.get('type')

    if main_input.strip() == '' and input_type in [
            'VMH', 'SwissLipids', 'ChEBI ID', 'ChEBI Name', 'PubChem ID', 'Saved']:
        return JsonResponse({'error': True, 'message': 'No input provided'})
    if input_type == 'VMH':
        BASE_URL = 'https://www.vmh.life/'
        endpoint = f"{BASE_URL}_api/metabolites/?abbreviation={main_input}"
        response = requests.get(endpoint, verify=False)
        if response.status_code == 200:
            data = response.json()
            res = data.get('results', [])
            if res:
                inchi_string, smile = res[0].get(
                    'inchiString', ''), res[0].get(
                    'smile', '')
                if not smile and not inchi_string:
                    return JsonResponse(
                        {
                            'error': True,
                            'message': f'Metabolite {main_input} does not have SMILES or inchi String on VMH'},
                        status=404)
                else:
                    return JsonResponse({'found': True,
                                         'abbr': data['results'][0]['abbreviation'],
                                         'name': data['results'][0]['fullName'],
                                         'miriam': data['results'][0]['miriam'],
                                         'input_type': input_type})
            else:
                return JsonResponse(
                    {'error': True, 'message': f'Metabolite `{main_input}` not found in VMH'}, status=404)
        else:
            return JsonResponse(
                {
                    'error': True,
                    'message': f'VMH API returned error {response.status_code} for metabolite `{main_input}`'},
                status=500)
        
    elif input_type == 'Saved':
        user_id = request.POST.get('userID')
        try:
            user = User.objects.get(pk=user_id)
        except (User.DoesNotExist, ValueError):
            return JsonResponse({'error': True, 'message': 'Cannot use saved metabolite without signing in'})

        # Look up the metabolite by primary key, checking if the user is the owner or shared with
        try:
            saved_met = SavedMetabolite.objects.get(
                Q(owner=user) | Q(shared_with=user),
                pk=main_input
            )
        except SavedMetabolite.DoesNotExist:
            return JsonResponse({'error': True, 'message': f'Saved metabolite with id `{main_input}` not found'}, status=404)

        return JsonResponse({
            'found': False,
            'abbr': saved_met.vmh_abbr,
            'name': saved_met.name,
            'input_type': input_type
        })
        
    else:
        mols, errors, names = any_to_mol(
            [main_input], [input_type], request, side=None)
        mol, error, name = mols[0], errors[0], names[0]
        if error:
            return JsonResponse({'error': True, 'message': error})
        found, miriam, abbr, nameVMH = search_vmh(
            mol, return_abbr=True, return_name=True)
        name = nameVMH if found else name
        name = capitalize_first_letter(name)
        inchi_key = Chem.MolToInchiKey(mol)
        user_id = request.POST.get('userID')
        if not user_id:
            user = None
        else:
            user = User.objects.get(pk=user_id)
        saved_exists, name_in_db = False, ''

        if user:
            # Check if the metabolite exists in the database for the current user (owner or shared)
            saved_met_obj = SavedMetabolite.objects.filter(
                Q(owner=user) | Q(shared_with=user),
                inchi_key=inchi_key
            ).distinct()
            
            if saved_met_obj.exists():
                saved_exists = True
                name_in_db = saved_met_obj.first().name

        return JsonResponse({
            'found': found,
            'abbr': abbr,
            'name': name,
            'miriam': miriam,
            'input_type': input_type,
            'saved_exists': saved_exists,
            'name_in_db': name_in_db
        })

def get_saved_metabolites(request):
    user_id = request.GET.get('user_id')
    if not user_id:
        return JsonResponse({'metabolites': []})

    user = User.objects.get(pk=user_id)
    
    # Get both owned and shared metabolites
    saved_metabolites = SavedMetabolite.objects.filter(
        Q(owner=user) | Q(shared_with=user)
    ).distinct()
    user_reactions = user.saved_reactions.all()
    metabolite_reactions = defaultdict(list)
    for reaction in user_reactions:
        try:
            subs = json.loads(reaction.substrates)
            prods = json.loads(reaction.products)
            subs_types = json.loads(reaction.substrates_types)
            prods_types = json.loads(reaction.products_types)
        except json.JSONDecodeError:
            subs, prods = [], []
        
        reaction_info = {
            'id': reaction.id,
            'name': reaction.short_name or f'Reaction {reaction.id}'
        }
        # Check all substrate/product IDs (stored as strings)
        for met_id, met_type in zip(subs + prods, subs_types + prods_types):
            if met_type == 'Saved':
                metabolite_reactions[str(met_id)].append(reaction_info)
    metabolites_list = []
    for met in saved_metabolites:
        external_links = {
            'keggId': met.keggId,
            'pubChemId': met.pubChemId,
            'cheBlId': met.cheBlId,
            'hmdb': met.hmdb,
            'metanetx': met.metanetx,
        }

        metabolites_list.append({
            'id': met.id,
            'name': met.name,
            'vmh_abbr': met.vmh_abbr or '',
            'reactions': metabolite_reactions[str(met.id)],
            'inchi_key': met.inchi_key, 
            'inchi': met.inchi,
            'smiles': met.smiles,
            'mol_w': met.mol_w,
            'mol_formula': met.mol_formula,
            'mol_file': met.mol_file,
            'external_links': external_links  # Include external links here
        })
    
    return JsonResponse({'metabolites': metabolites_list})

def fetch_rhea_rxn(request):
    if request.method == 'POST':
        try:
            data = json.loads(request.body)
            reaction_abbreviation = data.get('reactionAbbreviation')
            if reaction_abbreviation:
                result = get_from_rhea(reaction_abbreviation)
                return JsonResponse(result)
            else:
                return JsonResponse(
                    {'error': 'Missing reaction abbreviation'}, status=400)
        except json.JSONDecodeError:
            return JsonResponse({'error': 'Invalid JSON'}, status=400)
    else:
        return JsonResponse({'error': 'Invalid request'}, status=400)

# views.py
@csrf_exempt
def update_metabolite(request, metabolite_id):
    if request.method == 'PUT':
        try:
            metabolite = SavedMetabolite.objects.get(id=metabolite_id)
            data = json.loads(request.body)
            for field, value in data.items():
                if value.strip() == '' and field == 'name':
                    return JsonResponse({'status': 'error', 'message': f'{field} field cannot be empty'}, status=400)
                setattr(metabolite, field, value)
            metabolite.save()
            # Return updated data for the front end to patch the row without resetting the modal
            updated_metabolite = model_to_dict(metabolite, exclude=['shared_with'])
            return JsonResponse({'status': 'success', 'metabolite': updated_metabolite})
        except SavedMetabolite.DoesNotExist:
            return JsonResponse({'status': 'error', 'message': 'Metabolite not found'}, status=404)

@csrf_exempt
def delete_metabolite(request, metabolite_id):
    if request.method == 'DELETE':
        try:
            metabolite = SavedMetabolite.objects.get(id=metabolite_id)
            user_id = json.loads(request.body).get('user_id')
            user = User.objects.get(id=user_id) if user_id else None

            # Check if user is not the owner
            if metabolite.owner != user:
                if user in metabolite.shared_with.all():
                    metabolite.shared_with.remove(user)
                    metabolite.save()
                    return JsonResponse({'status': 'success', 'message': 'Metabolite removed from your list.'})
                else:
                    return JsonResponse({'status': 'error', 'message': 'Permission denied.'}, status=403)
            # User is the owner
            # Check if metabolite is shared and transfer ownership if possible
            shared_users = list(metabolite.shared_with.all())
            if shared_users:
                new_owner = shared_users[0]
                metabolite.owner = new_owner
                metabolite.shared_with.remove(new_owner)
                if user in metabolite.shared_with.all():
                    metabolite.shared_with.remove(user)
                metabolite.save()
                # Delete owner's reactions using this metabolite
                user_reactions = user.saved_reactions.all()
                affected_reactions = []
                for reaction in user_reactions:
                    subs = json.loads(reaction.substrates)
                    prods = json.loads(reaction.products)
                    if str(metabolite_id) in map(str, subs) or str(metabolite_id) in map(str, prods):
                        affected_reactions.append(reaction)
                # Remove reactions from user's saved_reactions
                user.saved_reactions.remove(*affected_reactions)
                return JsonResponse({'status': 'success', 'message': 'Ownership transferred and reactions removed.'})

            # Find all reactions containing this metabolite
            affected_reactions = []
            all_reactions = user.saved_reactions.all()
            for reaction in all_reactions:
                subs = json.loads(reaction.substrates)
                prods = json.loads(reaction.products)
                if str(metabolite_id) in map(str, subs) or str(metabolite_id) in map(str, prods):
                    affected_reactions.append(reaction)

            confirm = request.GET.get('confirm') == 'true'
            if not confirm and affected_reactions:
                reactions_list = [{'id': r.id, 'short_name': r.short_name} for r in affected_reactions]
                return JsonResponse({
                    'status': 'needs_confirmation',
                    'message': 'IMPORTANT: DELETING THIS METABOLITE WILL DELETE ALL REACTIONS THAT CONTAIN IT. ARE YOU SURE YOU WANT TO PROCEED?',
                    'reactions': reactions_list
                }, status=200)

            # Delete affected reactions and metabolite
            for reaction in affected_reactions:
                user.saved_reactions.remove(reaction)  # Use the ManyToManyField manager
            metabolite.delete()
            return JsonResponse({'status': 'success'})

        except SavedMetabolite.DoesNotExist:
            return JsonResponse({'status': 'error', 'message': 'Metabolite not found'}, status=404)
        except User.DoesNotExist:
            return JsonResponse({'status': 'error', 'message': 'User not found.'}, status=404)
        except Exception as e:
            return JsonResponse({'status': 'error', 'message': str(e)}, status=500)
        
def share_metabolites(request):
    try:
        data = json.loads(request.body)
        metabolite_ids = data.get('metabolite_ids', [])
        target_username = data.get('target_username')
        sharing_user_id = data.get('sharing_user_id')
        sharing_user = User.objects.get(pk=sharing_user_id)
        if not sharing_user:
            return JsonResponse({'status': 'error', 'message': 'Not logged in'}, status=401)

        if not metabolite_ids or not target_username:
            return JsonResponse({'status': 'error', 'message': 'Metabolite IDs and target username are required.'}, status=400)
        
        try:
            target_user = User.objects.get(name=target_username)

        except User.DoesNotExist:
            return JsonResponse({'status': 'error', 'message': 'Target user not found.'}, status=404)
        
        if target_user == sharing_user:
            return JsonResponse({'status': 'error', 'message': 'Cannot share metabolites with yourself.'}, status=400)
        #Check if the user is shared_with not owner for any of the metabolites
        metabolites_shared_with_user = SavedMetabolite.objects.filter(id__in=metabolite_ids, shared_with=sharing_user)
        if metabolites_shared_with_user.exists():
            names = [met.name for met in metabolites_shared_with_user]
            return JsonResponse({'status': 'error', 'message': 'Must be the owner to share metabolites. \nYou are not the owner of the following metabolites: ' + ', '.join(names)}, status=400)
        
        # Only share metabolites where the current user is the owner.
        metabolites = SavedMetabolite.objects.filter(id__in=metabolite_ids, owner=sharing_user)
        if not metabolites.exists():
            return JsonResponse({'status': 'error', 'message': 'No valid metabolites found to share.'}, status=400)
        shared_count = 0
        for metabolite in metabolites:
            if target_user not in metabolite.shared_with.all():
                metabolite.shared_with.add(target_user)
                shared_count += 1
        already_shared = len(metabolite_ids) - shared_count
        if already_shared > 0:
            msg = f'Shared {shared_count} metabolites with {target_username}. {already_shared} metabolites were already shared.'
        else:
            msg = f'Shared {shared_count} metabolites with {target_username}.'

        return JsonResponse({'status': 'success', 'message': msg})
    
    except Exception as e:
        return JsonResponse({'status': 'error', 'message': str(e)}, status=500)

def generate_abbreviation(request):
    if request.method == 'POST':
        # Implement your abbreviation generation logic here
        name = json.loads(request.body).get('name')
        generated_abbr = name[:3].upper()  # Example simple generation
        return JsonResponse({'abbr': generated_abbr})