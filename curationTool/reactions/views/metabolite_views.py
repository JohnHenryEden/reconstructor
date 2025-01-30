import json
import pandas as pd
from django.http import JsonResponse
from reactions.utils.search_vmh import search_vmh
from reactions.utils.to_mol import any_to_mol
from reactions.utils.utils import capitalize_first_letter
from reactions.utils.get_from_rhea import get_from_rhea
import requests


def verify_metabolite(request):
    """
    Verifies if a metabolite input from the user is in VMH.
    """
    main_input = request.POST.get('metabolite')
    input_type = request.POST.get('type')
    if main_input.strip() == '' and input_type in [
            'VMH', 'SwissLipids', 'ChEBI ID', 'ChEBI Name', 'PubChem ID']:
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
        return JsonResponse({'found': found,
                             'abbr': abbr,
                             'name': name,
                             'miriam': miriam,
                             'input_type': input_type})


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
