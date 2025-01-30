import json
from django.shortcuts import render
from django.http import JsonResponse
from reactions.utils.GPT_functions import get_vmh_met_from_inchi, metanetx_to_inchi, parse_metabolic_reactions_gpt, askGPT4
from reactions.views.user_views import validate_user_ID


def get_ai_response(request):
    if request.method == 'POST':
        userID = request.session.get('userID')
        user = validate_user_ID(userID)
        if not user.cred_add_to_vmh:
            return JsonResponse(
                {
                    'status': 'error',
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
    # Parsed GPT output into separate metabolite and gene reactions
    gpt_predictions_names = []
    # Function that asks ChatGPT to predict the reactions associated with the
    # gene
    reactions = askGPT4(gene, temperature)
    parsed_reactions = parse_metabolic_reactions_gpt(
        reactions)  # Function to parse ChatGPT's answer
    predicted_reaction_names = [
        reaction for reaction in parsed_reactions]  # Collect parsed reactions
    gpt_predictions_names.append(predicted_reaction_names)
    vmh_abbreviations = [
        [get_vmh_met_from_inchi(metanetx_to_inchi(met), met) for met in reaction]
        for reaction in gpt_predictions_names[0]
    ]
    return json.dumps(vmh_abbreviations)
   # return JsonResponse({'status': 'success', 'data': vmh_abbreviations})
