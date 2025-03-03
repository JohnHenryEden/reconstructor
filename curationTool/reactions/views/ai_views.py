"""
This module provides AI-based functionalities for predicting metabolic reactions using GPT models.

It includes:
- `get_ai_response`: Handles user requests for AI-generated predictions.
- `get_gpt_predictions`: Interacts with GPT to generate metabolic reaction predictions.
"""

import json
from django.shortcuts import render
from django.http import JsonResponse
from reactions.utils.GPT_functions import (
    get_vmh_met_from_inchi,
    metanetx_to_inchi,
    parse_metabolic_reactions_gpt,
    askGPT4
)
from reactions.views.user_views import validate_user_ID


def get_ai_response(request):
    """
    Handle POST requests to generate AI-based metabolic reaction predictions.

    Process:
        - Validates the user session and permissions.
        - Parses the input JSON containing a gene-related query.
        - Calls GPT-based functions to generate metabolic reaction predictions.
        - Returns the predictions in JSON format.

    Expects:
        JSON body with:
        - `key` (str): Input text containing gene names.
        - `temperature` (float, optional): Temperature parameter for GPT (default: 0.5).

    Returns:
        JsonResponse:
            - Success: Contains `predictions` (GPT-generated metabolic reactions).
            - Error: Returns an error message if the request is invalid or unauthorized.
    """
    if request.method == 'POST':
        user_id = request.session.get('userID')
        user = validate_user_ID(user_id)
        if not user.cred_add_to_vmh:
            return JsonResponse(
                {
                    'status': 'error',
                    'reason': 'permission_denied',
                    'error_message': 'User does not have permission to access this feature'},
                status=403)
        if user:
            try:
                # Parse the JSON data from the request body
                data = json.loads(request.body)
                # This should match the key you send in the JSON
                input_text = data.get('key')
                # Get temperature value, default to 0.5 if not provided
                temperature = data.get('temperature', 0.5)
                genes = input_text.split(' ')
                genes = [
                    gene for gene in genes if gene not in [
                        'and', 'or', 'AND', 'OR']]
                predictions = {}
                for gene in genes:
                    llm_response_html = get_gpt_predictions(gene, temperature)
                    predictions[gene] = llm_response_html
                return JsonResponse(
                    {'status': 'success', 'predictions': predictions})
            except json.JSONDecodeError:
                return JsonResponse(
                    {'status': 'error', 'error_message': 'Invalid JSON'}, status=400)
        else:
            return render(request, 'reactions/error.html',
                          {'error_message': 'Invalid key'})
    else:
        return JsonResponse(
            {'status': 'error', 'error_message': 'Unsupported method'}, status=405)


def get_gpt_predictions(gene, temperature):
    """
    Generate GPT-based metabolic reaction predictions for a given gene.

    Process:
        - Calls `askGPT4` to generate reaction predictions.
        - Parses the GPT output into structured metabolic reactions.
        - Converts metabolite identifiers to VMH abbreviations.

    Parameters:
        gene (str): The gene for which metabolic reactions are predicted.
        temperature (float): GPT temperature setting for response variability.

    Returns:
        str: A JSON string containing the list of VMH abbreviations for predicted reactions.
    """
    # Parsed GPT output into separate metabolite and gene reactions
    gpt_predictions_names = []
    # Function that asks ChatGPT to predict the reactions associated with the
    # gene
    reactions = askGPT4(gene, temperature)
    parsed_reactions = parse_metabolic_reactions_gpt(
        reactions)  # Function to parse ChatGPT's answer
    gpt_predictions_names.append(parsed_reactions)
    vmh_abbreviations = [
        [get_vmh_met_from_inchi(metanetx_to_inchi(met), met) for met in reaction]
        for reaction in gpt_predictions_names[0]
    ]
    return json.dumps(vmh_abbreviations)
   # return JsonResponse({'status': 'success', 'data': vmh_abbreviations})
