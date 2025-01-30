import json
import re
from django.shortcuts import render, redirect
from django.urls import reverse
from django.http import JsonResponse
from django.core import serializers
from django.views.decorators.csrf import csrf_exempt
from reactions.forms import ReactionForm
from reactions.models import User, CreatedReaction, Reaction, ReactionsAddedVMH
from django.views.decorators.http import require_POST
from reactions.utils.utils import safe_json_loads
from reactions.reaction_info import get_reaction_info
from reactions.utils.process_strings import construct_reaction_rxnfile
from reactions.utils.get_mol_info import get_mol_info
from reactions.utils.search_vmh import search_metabolites_vmh, check_reaction_vmh
from reactions.utils.to_mol import any_to_mol
from reactions.utils.utils import get_fields
from reactions.utils.RDT import RDT
from django.conf import settings
from django.views.decorators.http import require_GET
from .user_views import validate_user_ID
from django.shortcuts import get_object_or_404


def input_reaction(request):
    """
    Handles the POST request for a reaction input form.
    Processes the reaction data, RDT, Checks for the reaction in VMH,
    and returns the processed data.

    Input:
    - request: The HTTP request object from Django.

    Output:
    - HttpResponse: Renders a template with the reaction form or returns a JsonResponse with reaction data.
    """
    if request.method == 'POST':
        # Action to perform (either 'create' or 'edit')
        action = request.POST.get('action')
        form = ReactionForm(request.POST, request.FILES)
        user_id = request.POST.get('userID')
        if action == 'edit':
            reaction_id = request.POST.get('reaction_id')
            try:
                reaction = Reaction.objects.get(id=reaction_id)
            except Reaction.DoesNotExist:
                return JsonResponse(
                    {'message': 'Reaction not found', 'status': 'error'})
        else:
            # Create a new Reaction object if action is not 'edit'
            reaction = Reaction()

        # Get multiple substrates, products, and their stoichiometry as lists
        substrates_list = request.POST.getlist('substrates')
        products_list = request.POST.getlist('products')
        names_dict = request.POST.get('nameData')
        organs = request.POST.get('organs')

        names_dict = json.loads(names_dict)
        substrates_names = []
        products_names = []

        for key, value in names_dict.items():
            if 'substrate' in key:
                substrates_names.append(value)
            elif 'product' in key:
                products_names.append(value)
            else:
                raise ValueError(f"Invalid key: {key}")

        # Stoichiometry for substrates
        subs_sch = request.POST.getlist('subs_sch')
        prod_sch = request.POST.getlist(
            'prod_sch')  # Stoichiometry for products
        subs_comp = request.POST.getlist(
            'subs_comps')  # Compartments for substrates
        prod_comp = request.POST.getlist('prod_comps')
        substrates_types = request.POST.getlist('substrates_type')
        products_types = request.POST.getlist('products_type')
        direction = request.POST.get('direction')
        subs_sch = [int(s) for s in subs_sch]
        prod_sch = [int(s) for s in prod_sch]

        subs_mols, subs_errors, _ = any_to_mol(
            substrates_list, substrates_types, request, side='substrates')
        prod_mols, prod_errors, _ = any_to_mol(
            products_list, products_types, request, side='products')
        subsystem = request.POST.get('subsystem')

        all_errors = subs_errors + prod_errors
        if any(elem is not None for elem in all_errors):
            print(all_errors)

            error_message = "\n".join(
                [error for error in all_errors if error is not None])
            if request.headers.get('X-Requested-With') == 'XMLHttpRequest':
                return JsonResponse(
                    {'status': 'error', 'message': error_message})
            # Return error message in context for non-AJAX requests
            context = {'form': form, 'error_message': error_message}
            return render(request, 'reactions/Home_page.html', context)

        metabolite_formulas, metabolite_charges, metabolite_mol_file_strings, stereo_counts, stereo_locations_list = get_mol_info(
            subs_mols + prod_mols)
        metabolite_names = substrates_names + products_names
        reaction_rxn_file = construct_reaction_rxnfile(
            subs_mols, subs_sch, prod_mols, prod_sch, substrates_names, products_names)
        reaction.save()

        # Skip atom mapping if any product fields are empty or only one
        # substrate is provided
        skip_atom_mapping = request.POST.get('skipAtomMapping') == 'true' or (
            len(substrates_list) == 1 and len(products_list) == 0)
        if skip_atom_mapping:
            response_data = {'visualizations': [
                '/images/atom_mapping_skip.png']}
            balanced_count, (subs_atoms, prods_atoms), balanced_charge, (subs_charge,
                                                                         prods_charge), molc_formula, symb_to_name = get_reaction_info(reaction_rxn_file, direction)
        else:
            response_data = RDT(
                reaction_rxn_file,
                destination_path_png=f'media/images/visual{reaction.id}.png',
                destination_path_rxn=f'media/rxn_files/rxn{reaction.id}.rxn')
            balanced_count, (subs_atoms, prods_atoms), balanced_charge, (subs_charge, prods_charge), molc_formula, symb_to_name = get_reaction_info(
                f'media/rxn_files/rxn{reaction.id}.rxn', direction)

        vmh_found = check_reaction_vmh(
            substrates_list,
            products_list,
            subs_sch,
            prod_sch,
            substrates_types,
            products_types,
            subs_mols,
            prod_mols,
            direction,
            subsystem,
            subs_comp,
            prod_comp)

        if 'error' in response_data:
            context = {'form': form, 'error_message': response_data['error']}
            return render(request, 'reactions/Home_page.html', context)

        subs_found, subs_miriams = search_metabolites_vmh(
            substrates_list, substrates_types, request, side='substrates')
        prod_found, prod_miriams = search_metabolites_vmh(
            products_list, products_types, request, side='products')
        substrates_list = get_fields(
            request,
            substrates_list,
            substrates_types,
            settings.MEDIA_ROOT,
            settings.MEDIA_URL,
            side='substrates')
        products_list = get_fields(
            request,
            products_list,
            products_types,
            settings.MEDIA_ROOT,
            settings.MEDIA_URL,
            side='products')

        # Assign the values directly to the reaction instance
        reaction.Organs = json.dumps(organs)
        reaction.subs_sch = json.dumps(subs_sch)
        reaction.prods_sch = json.dumps(prod_sch)
        reaction.substrates_types = json.dumps(substrates_types)
        reaction.products_types = json.dumps(products_types)
        reaction.subs_comps = json.dumps(subs_comp)
        reaction.prods_comps = json.dumps(prod_comp)
        reaction.substrates_names = json.dumps(substrates_names)
        reaction.products_names = json.dumps(products_names)
        reaction.substrates = json.dumps(substrates_list)
        reaction.products = json.dumps(products_list)
        reaction.direction = direction
        reaction.subsystem = subsystem
        reaction.visualization = json.dumps(response_data['visualizations'])
        reaction.molc_formula = json.dumps([molc_formula])
        reaction.balanced_count = json.dumps([balanced_count])
        reaction.balanced_charge = json.dumps([balanced_charge])
        reaction.subs_atoms = json.dumps([subs_atoms])
        reaction.prods_atoms = json.dumps([prods_atoms])
        reaction.subs_charge = json.dumps([subs_charge])
        reaction.prods_charge = json.dumps([prods_charge])
        reaction.symb_to_name = json.dumps([symb_to_name])
        reaction.subs_found = json.dumps(subs_found)
        reaction.subs_miriams = json.dumps(subs_miriams)
        reaction.prod_found = json.dumps(prod_found)
        reaction.prod_miriams = json.dumps(prod_miriams)
        reaction.vmh_found = vmh_found['found']
        reaction.vmh_found_similar = vmh_found['similar']
        reaction.vmh_url = json.dumps(
            vmh_found['url']) if vmh_found['found'] else None
        reaction.vmh_formula = json.dumps(
            vmh_found['formula']) if vmh_found['found'] else None
        reaction.metabolite_names = json.dumps(metabolite_names)
        reaction.metabolite_formulas = json.dumps(metabolite_formulas)
        reaction.metabolite_charges = json.dumps(metabolite_charges)
        reaction.metabolite_mol_file_strings = json.dumps(
            metabolite_mol_file_strings)
        reaction.stereo_counts = json.dumps(stereo_counts)
        reaction.stereo_locations_list = json.dumps(stereo_locations_list)
        reaction.save()
        if request.headers.get('X-Requested-With') == 'XMLHttpRequest':
            data = {
                'visualization': json.loads(reaction.visualization),
                'molc_formula': json.loads(reaction.molc_formula),
                'balanced_count': json.loads(reaction.balanced_count),
                'balanced_charge': json.loads(reaction.balanced_charge),
                'subs_atoms': json.loads(reaction.subs_atoms),
                'prods_atoms': json.loads(reaction.prods_atoms),
                'subs_charge': json.loads(reaction.subs_charge),
                'prods_charge': json.loads(reaction.prods_charge),
                'symb_to_name': json.loads(reaction.symb_to_name),
                'subs_found': json.loads(reaction.subs_found),
                'subs_miriams': json.loads(reaction.subs_miriams),
                'prod_found': json.loads(reaction.prod_found),
                'prod_miriams': json.loads(reaction.prod_miriams),
                'reaction_id': reaction.id,
                'metabolite_names': json.loads(reaction.metabolite_names),
                'metabolite_formulas': json.loads(reaction.metabolite_formulas),
                'metabolite_charges': json.loads(reaction.metabolite_charges),
                'metabolite_mol_file_strings': json.loads(reaction.metabolite_mol_file_strings),
                'sterio_counts': json.loads(reaction.stereo_counts),
                'sterio_locations_list': json.loads(reaction.stereo_locations_list),
            }
            if vmh_found['found']:
                data['vmh_found'] = vmh_found['found']
                data['vmh_found_similar'] = vmh_found['similar']
                data['vmh_url'] = vmh_found['url']
                data['vmh_formula'] = vmh_found['formula']
            data['status'] = 'success'
            return JsonResponse(data)
    else:
        form = ReactionForm()
    return render(request, 'reactions/Home_page.html', {'form': form})


def get_reaction(request, reaction_id):
    """
    Fetches and returns details of a reaction by its ID.

    :param request: The HTTP request object.
    :param reaction_id: The ID of the reaction to retrieve.
    :return: JsonResponse containing the reaction details or an error message.
    """
    try:
        reaction = Reaction.objects.get(pk=reaction_id)
        reaction_data = {
            'Organs': reaction.Organs,
            'reaction_id': reaction.id,
            'short_name': reaction.short_name,
            'substrates': safe_json_loads(reaction.substrates),
            'products': safe_json_loads(reaction.products),
            'substrates_names': safe_json_loads(reaction.substrates_names),
            'products_names': safe_json_loads(reaction.products_names),
            'direction': reaction.direction,
            'subsystem': reaction.subsystem,
            'subs_comps': safe_json_loads(reaction.subs_comps),
            'prods_comps': safe_json_loads(reaction.prods_comps),
            'visualization': safe_json_loads(reaction.visualization),
            'rxn_formula': safe_json_loads(reaction.rxn_formula),
            'molc_formula': safe_json_loads(reaction.molc_formula),
            'balanced_count': safe_json_loads(reaction.balanced_count),
            'balanced_charge': safe_json_loads(reaction.balanced_charge),
            'subs_sch': safe_json_loads(reaction.subs_sch),
            'prod_sch': safe_json_loads(reaction.prods_sch),
            'subs_types': safe_json_loads(reaction.substrates_types),
            'prods_types': safe_json_loads(reaction.products_types),
            'subs_atoms': safe_json_loads(reaction.subs_atoms),
            'prods_atoms': safe_json_loads(reaction.prods_atoms),
            'subs_charge': safe_json_loads(reaction.subs_charge),
            'prods_charge': safe_json_loads(reaction.prods_charge),
            'symb_to_name': safe_json_loads(reaction.symb_to_name),
            'subs_found': safe_json_loads(reaction.subs_found),
            'subs_miriams': safe_json_loads(reaction.subs_miriams),
            'prod_found': safe_json_loads(reaction.prod_found),
            'prod_miriams': safe_json_loads(reaction.prod_miriams),
            'vmh_found': reaction.vmh_found,
            'vmh_found_similar': reaction.vmh_found_similar,
            'vmh_url': reaction.vmh_url,
            'vmh_formula': reaction.vmh_formula,
            'metabolite_names': safe_json_loads(reaction.metabolite_names),
            'metabolite_formulas': safe_json_loads(reaction.metabolite_formulas),
            'metabolite_charges': safe_json_loads(reaction.metabolite_charges),
            'metabolite_mol_file_strings': safe_json_loads(reaction.metabolite_mol_file_strings),
            'stereo_counts': safe_json_loads(reaction.stereo_counts),
            'stereo_locations_list': safe_json_loads(reaction.stereo_locations_list),
        }
        print(reaction_data['stereo_counts'])
        print(reaction_data['stereo_locations_list'])

        return JsonResponse(reaction_data)

    except Reaction.DoesNotExist:
        return JsonResponse({'error': 'Reaction not found'}, status=404)

    except Exception as e:
        return JsonResponse(
            {'error': 'An unexpected error occurred'}, status=500)


@csrf_exempt
def add_info_to_reaction(request):
    if request.method == 'POST':
        try:
            data = json.loads(request.body)
            user_id = data.get('userID')
            info_type = data.get('infoType')
            info_text = data.get('infoText')
            ext_link_type = data.get('extLinkType', '')
            ref_type = data.get('refType', '')
            reaction_id = data.get('reactionId')

            if not all([user_id, info_type, info_text]) or (
                    info_type != 'Gene Info' and not reaction_id):
                return JsonResponse({'status': 'error',
                                     'message': 'All fields are required.',
                                     'info_type': info_type},
                                    status=400)

            user = User.objects.get(pk=user_id)
            username = user.name

            if info_type == 'Gene Info' and not reaction_id:
                # Store gene info in session if no reaction_id is provided
                if 'gene_info' not in request.session:
                    request.session['gene_info'] = []

                info_data = {'info': info_text, 'user_name': username}
                request.session['gene_info'].append(info_data)
                request.session.modified = True

                return JsonResponse({'status': 'success',
                                     'message': 'Gene information added to session.',
                                     'info_type': info_type})

            # Fetch the reaction if reaction_id is provided
            reaction = Reaction.objects.get(id=reaction_id)

            if info_type == 'Reference':
                if reaction.references is None:
                    reaction.references = []
                info_data_template = {
                    'user_name': username, 'ref_type': ref_type}

                # Determine the delimiter based on the reference type
                if 'PMID' in ref_type:
                    ref_list = [ref.strip() for ref in info_text.split(';')]
                    existing_refs = {
                        ref['info'] for ref in reaction.references if 'PMID' in ref['info']}
                elif 'DOI' in ref_type:
                    ref_list = [ref.strip() for ref in info_text.split(',')]
                    existing_refs = {
                        ref['info'] for ref in reaction.references if 'DOI' in ref['info']}
                else:
                    ref_list = [info_text.strip()]
                    existing_refs = {ref['info']
                                     for ref in reaction.references}

                for ref in ref_list:
                    if ref not in existing_refs:
                        info_data = info_data_template.copy()
                        info_data['info'] = ref
                        reaction.references.append(info_data)
                        existing_refs.add(ref)

            elif info_type == 'External Link':
                if reaction.ext_links is None:
                    reaction.ext_links = []
                info_data = {
                    'info': info_text,
                    'user_name': username,
                    'ext_link_type': ext_link_type}
                reaction.ext_links.append(info_data)

            elif info_type == 'Gene Info':
                if reaction.gene_info is None:
                    reaction.gene_info = []
                info_data = {'info': info_text, 'user_name': username}
                reaction.gene_info.append(info_data)

            elif info_type == 'Comment':
                if reaction.comments is None:
                    reaction.comments = []
                info_data = {'info': info_text, 'user_name': username}
                reaction.comments.append(info_data)

            reaction.save()
            return JsonResponse({'status': 'success',
                                 'message': 'Information added successfully.',
                                 'reaction_id': reaction_id,
                                 'info_type': info_type})

        except User.DoesNotExist:
            return JsonResponse({'status': 'error',
                                 'message': 'Invalid user key.',
                                 'info_type': info_type},
                                status=404)
        except Reaction.DoesNotExist:
            return JsonResponse({'status': 'error',
                                 'message': 'Invalid reaction ID.',
                                 'info_type': info_type},
                                status=404)
        except json.JSONDecodeError:
            return JsonResponse({'status': 'error',
                                 'message': 'Invalid JSON.',
                                 'info_type': info_type},
                                status=400)
        except Exception as e:
            return JsonResponse({'status': 'error', 'message': str(
                e), 'info_type': info_type}, status=500)
    else:
        return JsonResponse({'status': 'error',
                             'message': 'Invalid request method.',
                             'info_type': info_type},
                            status=405)


def update_gene_info(request):
    try:
        data = json.loads(request.body)
        user_id = data.get('userID')
        reaction_id = data.get('reactionID')
        gene = data.get('gene')
        field_type = data.get('fieldType')
        updated_value = data.get('updatedValue')

        # Debugging statements to trace values

        # Check if any required field is missing or updated_value is empty
        if not all([user_id, reaction_id, gene, field_type]
                   ) or updated_value is None or updated_value.strip() == "":
            return JsonResponse(
                {
                    'status': 'error',
                    'message': 'All fields are required and updated value cannot be empty.'},
                status=400)

        user = User.objects.get(pk=user_id)
        reaction = Reaction.objects.get(id=reaction_id)

        if reaction.gene_info is None:
            return JsonResponse(
                {'status': 'error', 'message': 'No gene information found.'}, status=404)

        gene_info_updated = False
        for gene_info in reaction.gene_info:
            if gene in gene_info['info']:
                if field_type == 'Organs':
                    organ_match = re.search(
                        r'ORGAN\(([^)]+)\)', gene_info['info'])
                    if organ_match:
                        old_organs = organ_match.group(1)
                        new_info = gene_info['info'].replace(
                            f'ORGAN({old_organs})', f'ORGAN({updated_value})')
                        gene_info['info'] = new_info
                        gene_info_updated = True
                elif field_type == 'SubcellularLocations':
                    subcellular_match = re.search(
                        r'SUBCELLULAR\(([^)]+)\)', gene_info['info'])
                    if subcellular_match:
                        old_subcellular = subcellular_match.group(1)
                        new_info = gene_info['info'].replace(
                            f'SUBCELLULAR({old_subcellular})', f'SUBCELLULAR({updated_value})')
                        gene_info['info'] = new_info
                        gene_info_updated = True

        if gene_info_updated:
            reaction.save()
            return JsonResponse(
                {'status': 'success', 'message': 'Gene information updated successfully.'})
        else:
            return JsonResponse(
                {'status': 'error', 'message': 'Gene information not found or not updated.'}, status=404)

    except User.DoesNotExist:
        return JsonResponse(
            {'status': 'error', 'message': 'Invalid user key.'}, status=404)
    except Reaction.DoesNotExist:
        return JsonResponse(
            {'status': 'error', 'message': 'Invalid reaction ID.'}, status=404)
    except json.JSONDecodeError:
        return JsonResponse(
            {'status': 'error', 'message': 'Invalid JSON.'}, status=400)
    except Exception as e:
        return JsonResponse({'status': 'error', 'message': str(e)}, status=500)


@csrf_exempt
def delete_reaction_info(request):
    if request.method == 'POST':
        try:
            req_body = json.loads(request.body)
            reaction_id = req_body.get('reaction_id')
            tab_id = req_body.get('tab_id')
            item_to_delete = req_body.get('item_to_delete')
            # Fetch the reaction object
            react_obj = Reaction.objects.get(pk=reaction_id)

            # Handle deletion based on the tab_id
            if tab_id == 'refs-content':
                for ref in react_obj.references:
                    if ref['info'] == item_to_delete['info'] and ref['ref_type'] == item_to_delete['ref_type']:
                        react_obj.references.remove(ref)
                        react_obj.save()
                        return JsonResponse(
                            {'status': 'success', 'message': 'Reference deleted successfully'})

            elif tab_id == 'ext-links-content':
                for ext_link in react_obj.ext_links:
                    if ext_link['info'] == item_to_delete['info'] and ext_link['ext_link_type'] == item_to_delete['ext_link_type']:
                        react_obj.ext_links.remove(ext_link)
                        react_obj.save()
                        return JsonResponse(
                            {'status': 'success', 'message': 'External link deleted successfully'})

            elif tab_id == 'comments-content':
                for comment in react_obj.comments:
                    if comment['info'] == item_to_delete['info']:
                        react_obj.comments.remove(comment)
                        react_obj.save()
                        return JsonResponse(
                            {'status': 'success', 'message': 'Comment deleted successfully'})

            elif tab_id == 'gene-info-content':
                for gene_info in react_obj.gene_info:
                    if gene_info['info'] == item_to_delete['info']:
                        react_obj.gene_info.remove(gene_info)
                        react_obj.save()
                        return JsonResponse(
                            {'status': 'success', 'message': 'Gene info deleted successfully'})

            # If no valid tab_id matched
            return JsonResponse(
                {'error': True, 'message': 'Invalid tab_id or item not found'}, status=400)

        except Reaction.DoesNotExist:
            return JsonResponse(
                {'error': True, 'message': 'Reaction not found'}, status=404)
        except Exception as e:
            return JsonResponse({'error': True, 'message': str(e)}, status=500)
    else:
        return JsonResponse(
            {'error': True, 'message': 'Invalid request method'}, status=400)


def get_reaction_details(request):
    try:
        reaction_id = json.loads(request.body)
        reaction = Reaction.objects.get(pk=reaction_id)

        if reaction.references is None:
            reaction.references = []
        if reaction.ext_links is None:
            reaction.ext_links = []
        if reaction.gene_info is None:
            reaction.gene_info = []
        if reaction.comments is None:
            reaction.comments = []

        return JsonResponse({
            'references': reaction.references,
            'external_links': reaction.ext_links,
            'gene_info': reaction.gene_info,
            'comments': reaction.comments,
        })
    except Reaction.DoesNotExist:
        return JsonResponse({'error': 'Reaction not found'}, status=404)


def delete_reaction(request):
    if request.method == 'POST':
        userID = request.POST.get('userID')
        reaction_id = request.POST.get('reaction_id')
        user = validate_user_ID(userID)
        if user:
            reaction = Reaction.objects.get(pk=reaction_id)
            user.saved_reactions.remove(reaction)
            saved_reactions_url = reverse('saved_reactions')
            return redirect(saved_reactions_url)
    return render(request, '404.html')


def saved_reactions(request, modal=False):
    userID = request.session.get('userID')
    user = validate_user_ID(userID)
    if user:
        reactions = user.saved_reactions.all().order_by('id')
        reactions_json = serializers.serialize('json', reactions)
        combined_reactions_details = []

        for reaction in reactions:
            # Parse JSON fields
            subs_sch = json.loads(reaction.subs_sch)
            prods_sch = json.loads(reaction.prods_sch)
            subs_comps = json.loads(reaction.subs_comps)
            prods_comps = json.loads(reaction.prods_comps)
            short_name = reaction.short_name

            # Construct details strings
            subs_details = " + ".join(
                f"{float(sch)} {json.loads(reaction.substrates_names)[idx]} [{comp}]" for idx,
                (sch,
                 comp) in enumerate(
                    zip(
                        subs_sch,
                        subs_comps)))
            prods_details = " + ".join(
                f"{float(sch)} {json.loads(reaction.products_names)[idx]} [{comp}]" for idx,
                (sch,
                 comp) in enumerate(
                    zip(
                        prods_sch,
                        prods_comps)))

            # Check if gene_info is None
            if reaction.gene_info:
                if isinstance(reaction.gene_info, str):
                    gene_info_data = json.loads(reaction.gene_info)
                else:
                    gene_info_data = reaction.gene_info

                gene_info_list = []
                for gene in gene_info_data:
                    if 'info' in gene:
                        info = gene['info']
                        gpr_start = info.find('GPR: ')
                        if gpr_start != -1:
                            gpr_end = info.find(';', gpr_start)
                            if gpr_end != -1:
                                gpr_info = info[gpr_start + 5:gpr_end]
                                gene_info_list.append(gpr_info.strip())
                        else:
                            gene_info_list.append(info.strip())
            else:
                gene_info_list = []

            # Get associated flags and their colors
            flags = reaction.flags.all()
            flag_details = [{"name": flag.name_flag,
                             "color": flag.color} for flag in flags]

            combined_reactions_details.append({
                'reaction': reaction,
                'details': {
                    'short_name': short_name,
                    'subs_details': subs_details,
                    'prods_details': prods_details,
                    'molc_formula': reaction.molc_formula,
                    'balanced_count': json.loads(reaction.balanced_count)[0] if reaction.balanced_count else None,
                    'balanced_charge': json.loads(reaction.balanced_charge)[0] if reaction.balanced_charge else None,
                    'subsystem': reaction.subsystem,
                    'direction': reaction.direction,
                    'gene_info': gene_info_list,
                    'flags': flag_details,  # Include flag details with name and color
                }
            })

        context = {
            'reactions': reactions,
            'reactions_json': reactions_json,
            'userID': userID,
            'user_name': user.name,
            'combined_reactions_details': combined_reactions_details
        }

        if modal:
            return render(
                request,
                'reactions/saved_reactions_modal.html',
                context)
        else:
            return render(request, 'reactions/saved_reactions.html', context)
    else:
        return render(request, 'reactions/error.html',
                      {'error_message': 'Invalid key'})


def search_reactions(request):
    user = request.user
    query = request.GET.get('q', '')
    if query:
        reactions = user.saved_reactions.filter(
            substrates__icontains=query) | user.saved_reactions.filter(
            products__icontains=query) | user.saved_reactions.filter(
            short_name__icontains=query) | user.saved_reactions.filter(
            direction__icontains=query)
    else:
        reactions = user.saved_reactions.all()
    reactions_data = [
        {
            'substrates': reaction.substrates,
            'products': reaction.products,
            'short_name': reaction.short_name,
            'direction': reaction.direction,
            # Add more fields if needed
        }
        for reaction in reactions
    ]
    return JsonResponse({'reactions': reactions_data})


def identical_reaction(request):
    """
    Checks if the user already has a saved reaction that matches the
    posted substrates/products (and their stoichiometries).
    Returns a JSON response with 'exists', 'reaction_id', and 'reaction_name' if found.
    """
    user_id = request.POST.get('userID')
    if not user_id:
        # If no userID is provided, we consider it "not found."
        return JsonResponse({'exists': False, 'status': 'success'})

    # Attempt to load the user
    try:
        user = User.objects.get(pk=user_id)
    except User.DoesNotExist:
        # User not found => no saved reactions
        return JsonResponse({'exists': False, 'status': 'success'})

    # Gather substrates/products and their stoichiometries from form data
    substrates_list = request.POST.getlist('substrates')  # e.g. ["H2O", "ATP"]
    subs_sch_list = request.POST.getlist('subs_sch')    # e.g. ["1", "1"]
    products_list = request.POST.getlist('products')    # e.g. ["ADP", "Pi"]
    prods_sch_list = request.POST.getlist('prod_sch')    # e.g. ["1", "1"]

    # Convert to JSON-strings as stored in Reaction model
    # (Make sure subs_sch and prod_sch are converted to int before JSON-serializing
    #  if that's how they are stored in the Reaction model.)
    substrates_str = json.dumps(substrates_list)
    subs_sch_str = json.dumps([int(s) for s in subs_sch_list])
    products_str = json.dumps(products_list)
    prods_sch_str = json.dumps([int(s) for s in prods_sch_list])

    matching_reactions = user.saved_reactions.filter(
        substrates=substrates_str,
        subs_sch=subs_sch_str,
        products=products_str,
        prods_sch=prods_sch_str
    )

    if matching_reactions.exists():
        matches_data = []
        for r in matching_reactions:
            r_name = r.short_name or f"Reaction {r.id}"
            matches_data.append({
                'reaction_id': r.id,
                'reaction_name': r_name
            })
        return JsonResponse({
            'exists': True,
            'matches': matches_data,
            'status': 'success'
        })
    else:
        return JsonResponse({'exists': False, 'status': 'success'})


@require_POST
def clone_reaction_view(request):
    try:
        reaction_id = request.POST.get('reaction_id')
        user_id = request.POST.get('userID')
        reaction_name = request.POST.get('name')
        cloned_reaction = Reaction.objects.get(pk=reaction_id)

        cloned_reaction.pk = None  # Set the primary key to None to create a new instance
        cloned_reaction.short_name = reaction_name
        # Save the cloned reaction object to generate a new ID
        cloned_reaction.save()
        user = User.objects.get(pk=user_id)
        user.saved_reactions.add(cloned_reaction)
        user.save()
        return JsonResponse({'status': 'success',
                             'message': 'Reaction cloned successfully'})
    except Exception as e:
        return JsonResponse({'status': 'error', 'message': str(e)})


@require_GET
def get_user_reactions_and_vmh(request):
    """Return all existing users' full names, the number of reactions saved, reactions added to VMH, and reactions created."""
    try:
        users_data = []
        users = User.objects.all()

        for user in users:
            saved_reactions_count = user.saved_reactions.count()
            reactions_added_vmh_count = ReactionsAddedVMH.objects.filter(
                user=user).count()
            created_reactions_count = CreatedReaction.objects.filter(
                user=user).count()

            user_data = {
                'full_name': user.full_name,
                'saved': saved_reactions_count,
                'added': reactions_added_vmh_count,
                'created': created_reactions_count
            }
            users_data.append(user_data)

        return JsonResponse(users_data, safe=False)

    except Exception as e:
        return JsonResponse({'error': str(e)}, status=500)


@csrf_exempt
def create_reaction(request):
    if request.method == 'POST':
        data = json.loads(request.body)
        user_id = data.get('user_id')
        reaction_id = data.get('reaction_id')
        # Validate the user ID using the provided function
        user = validate_user_ID(user_id)

        if not user:
            return JsonResponse(
                {'success': False, 'error': 'User does not exist'}, status=400)

        # Fetch the reaction based on reaction_id
        reaction = get_object_or_404(Reaction, id=reaction_id)

        created_reaction = CreatedReaction.objects.create(
            user=user, reaction=reaction)
        return JsonResponse({'success': True,
                             'created_reaction_id': created_reaction.id})
    else:
        return JsonResponse(
            {'success': False, 'error': 'Invalid request method'}, status=400)


@csrf_exempt
def reaction_name_exists(request):
    # Check if the given reaction name already exists in the users saved
    # reactions
    if request.method == 'POST':
        try:
            data = json.loads(request.body)
            user_id = data.get('user_id')
            reaction_name = data.get('short_name')
            user = User.objects.get(pk=user_id)
            saved_reaction_names = list(
                user.saved_reactions.all().values_list(
                    'short_name', flat=True))
            is_name_saved = reaction_name in saved_reaction_names
            return JsonResponse({'is_name_saved': is_name_saved})
        except User.DoesNotExist:
            return JsonResponse({'error': 'User does not exist'}, status=404)
        except (KeyError, ValueError):
            return JsonResponse({'error': 'Invalid data'}, status=400)


def get_available_reactions(request):
    """Return the IDs of available reactions for the given user."""
    if request.method == 'POST':
        try:
            data = json.loads(request.body)
            user_id = data.get('user_id')

            user = User.objects.get(pk=user_id)
            saved_reaction_ids = list(
                user.saved_reactions.all().values_list(
                    'id', flat=True))
            available_reactions = Reaction.objects.exclude(
                id__in=saved_reaction_ids)
            available_reaction_ids = list(
                available_reactions.values_list(
                    'id', flat=True))
            last_index = saved_reaction_ids[-1] if saved_reaction_ids else None

            return JsonResponse(
                {'available_reaction_ids': available_reaction_ids, 'last_index': last_index})
        except User.DoesNotExist:
            return JsonResponse({'error': 'User does not exist'}, status=404)
        except (KeyError, ValueError):
            return JsonResponse({'error': 'Invalid data'}, status=400)
    return JsonResponse({'error': 'Invalid request method'}, status=405)


@csrf_exempt
def already_saved(request):
    """Check if the given reaction_id is in the user's saved reactions."""
    if request.method == 'POST':
        try:
            data = json.loads(request.body)
            user_id = data.get('user_id')
            # Convert reaction_id to integer
            reaction_id = int(data.get('reaction_id'))

            user = User.objects.get(pk=user_id)
            saved_reaction_ids = list(
                user.saved_reactions.all().values_list(
                    'id', flat=True))

            is_reaction_saved = reaction_id in saved_reaction_ids

            return JsonResponse({'is_reaction_saved': is_reaction_saved})
        except User.DoesNotExist:
            return JsonResponse({'error': 'User does not exist'}, status=404)
        except (KeyError, ValueError):
            return JsonResponse({'error': 'Invalid data'}, status=400)

    return JsonResponse({'error': 'Invalid request method'}, status=405)
