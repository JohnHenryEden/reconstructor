from django.http import JsonResponse
import json
from django.views.decorators.csrf import csrf_exempt
from reactions.models import Reaction
import re
import requests
from reactions.utils.utils import parse_xml
from django.shortcuts import render


def check_session_data(request):
    # Check if 'gene_info' is in session
    gene_info = request.session.get('gene_info', None)

    if gene_info:
        return JsonResponse({'status': 'success', 'gene_info': gene_info})
    else:
        return JsonResponse(
            {'status': 'error', 'message': 'No gene info in session.'})


# Use csrf_exempt if you don't want to deal with CSRF tokens in
# development. However, be cautious with this in production.
@csrf_exempt
def delete_gene_info_from_session(request):
    if request.method == 'POST':
        try:
            data = json.loads(request.body)
            info_to_delete = data.get('info_to_delete')

            # Get the gene info list from the session
            gene_info = request.session.get('gene_info', [])

            # Filter out the item to delete
            updated_gene_info = [
                info for info in gene_info if info['info'] != info_to_delete]

            # Update the session with the filtered list
            request.session['gene_info'] = updated_gene_info
            request.session.modified = True

            return JsonResponse({'status': 'success',
                                 'message': 'Gene info deleted from session.'})
        except Exception as e:
            return JsonResponse(
                {'status': 'error', 'message': str(e)}, status=400)
    else:
        return JsonResponse(
            {'status': 'error', 'message': 'Invalid request method.'}, status=405)


# Use csrf_exempt if you don't want to deal with CSRF tokens in
# development. However, be cautious with this in production.
@csrf_exempt
def clear_session(request):
    if request.method == 'POST':
        try:
            # Clear the entire session
            request.session.flush()  # This will clear all session data

            return JsonResponse({'status': 'success',
                                 'message': 'Session cleared successfully.'})
        except Exception as e:
            return JsonResponse(
                {'status': 'error', 'message': str(e)}, status=400)
    else:
        return JsonResponse(
            {'status': 'error', 'message': 'Invalid request method.'}, status=405)


@csrf_exempt
def parse_formula_with_compartments(request):
    if request.method == 'POST':
        try:
            body = json.loads(request.body)
        except json.JSONDecodeError as e:
            return JsonResponse(
                {'error': f'JSON decode error: {str(e)}'}, status=400)

        formulas = body.get('formulas', [])
        subs_comps = body.get('subs_comps', [])
        prods_comps = body.get('prods_comps', [])

        # Log received data for debugging

        detailed_formulas = []

        def extract_components(comp_str):
            """Extracts quantity and name from a component string."""
            match = re.match(r'(\d*)\s*(\w+)', comp_str.strip())
            if match:
                quantity = match.group(1) if match.group(1) else '1'
                name = match.group(2)
                return quantity, name
            return None, None

        for formula in formulas:
            # Normalize formula by ensuring spaces around the arrows
            formula = formula.replace('->', ' -> ').replace('<=>', ' <=> ')

            if ' -> ' in formula:
                direction = 'forward'
                substrates, products = formula.split(' -> ')
            elif ' <=> ' in formula:
                direction = 'reversible'
                substrates, products = formula.split(' <=> ')
            else:
                return JsonResponse(
                    {'error': f'Invalid formula format: {formula}'}, status=400)

            subs_list = [comp.strip() for comp in substrates.split(
                ' + ')] if substrates else []
            prods_list = [comp.strip()
                          for comp in products.split(' + ')] if products else []

            # Remove empty components
            subs_list = [comp for comp in subs_list if comp]
            prods_list = [comp for comp in prods_list if comp]

            # Check if the lengths of subs_comps and prods_comps match the
            # lengths of subs_list and prods_list
            if len(subs_list) != len(subs_comps):
                return JsonResponse(
                    {'error': 'Number of substrate compartments does not match number of substrates'}, status=400)

            if len(prods_list) != len(prods_comps):
                return JsonResponse(
                    {'error': 'Number of product compartments does not match number of products'}, status=400)

            subs_detailed = []
            for j, comp in enumerate(subs_list):
                quantity, name = extract_components(comp)
                subs_detailed.append(f"{quantity} {name}[{subs_comps[j]}]")

            prods_detailed = []
            for j, comp in enumerate(prods_list):
                quantity, name = extract_components(comp)
                prods_detailed.append(f"{quantity} {name}[{prods_comps[j]}]")

            if direction == 'forward':
                detailed_formula = ' + '.join(subs_detailed) + \
                    ' -> ' + ' + '.join(prods_detailed)
            else:  # reversible
                detailed_formula = ' + '.join(subs_detailed) + \
                    ' <=> ' + ' + '.join(prods_detailed)

            detailed_formulas.append(detailed_formula)

        return JsonResponse({'detailed_formulas': detailed_formulas})

    return JsonResponse({'error': 'Invalid request method'}, status=400)


@csrf_exempt  # Use csrf_exempt if CSRF token isn't being managed
def save_formula(request):
    if request.method == 'POST':
        try:
            data = json.loads(request.body)
            formula = data.get('formula', '')
            reaction_id = data.get('reaction_id', '')

            # Ensure the formula is saved with quotes
            quoted_formula = f'"{formula}"'

            # Find the reaction by the provided reaction ID
            try:
                reaction = Reaction.objects.get(id=reaction_id)
            except Reaction.DoesNotExist:
                return JsonResponse(
                    {'status': 'failed', 'error': 'Reaction not found'}, status=404)

            # Update the formulas field with the provided formula string with
            # quotes
            reaction.rxn_formula = quoted_formula
            reaction.save()

            return JsonResponse({'status': 'success'})
        except Exception as e:
            return JsonResponse(
                {'status': 'failed', 'error': str(e)}, status=500)
    return JsonResponse(
        {'status': 'failed', 'error': 'Invalid request method'}, status=400)


def get_pubmed_info(request, pmid):
    # Base URL for PubMed API
    base_url = 'https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi'
    params = {
        'db': 'pubmed',
        'id': pmid,
        'retmode': 'xml'
    }
    response = requests.get(base_url, params=params)

    if response.status_code == 200:
        # Parsing XML response is necessary here to extract needed information
        # This is a placeholder to indicate where parsing should occur
        # You might use libraries like xml.etree.ElementTree or lxml to parse
        # the XML

        pubmed_info = parse_xml(response.content)  # Implement this function
        if pubmed_info['title'] is None:
            return JsonResponse(
                {'status': 'error', 'message': 'PubMed API returned an empty response'}, status=500)
        return JsonResponse({
            'status': 'success',
            'author': pubmed_info['authors'],
            'title': pubmed_info['title'],
            'abstract': pubmed_info['abstract']
        })
    else:
        return JsonResponse({'status': 'error',
                             'message': f"PubMed API returned error {response.status_code}"},
                            status=500)


def get_doi_info(request, doi):
    doi = doi.replace('DOI:', '')
    base_url = f'https://api.crossref.org/works/{doi}'
    try:
        response = requests.get(base_url)
        if response.status_code == 200:
            data = response.json()['message']

            # Extract needed information
            authors = [
                f"{author['given']} {author['family']}" for author in data.get(
                    'author', [])]
            title = data.get('title', [None])[0]
            abstract = data.get('abstract', None)
            return JsonResponse({
                'status': 'success',
                'author': authors,
                'title': title,
                'abstract': abstract
            })
        else:
            return JsonResponse({'status': 'error',
                                 'message': 'CrossRef API returned an error'},
                                status=response.status_code)
    except Exception as e:
        return JsonResponse(
            {'status': 'error', 'message': 'Failed to fetch DOI data'}, status=500)


def chemdoodle_sketcher(request):
    return render(request, 'reactions/chemdoodle_sketcher.html')
