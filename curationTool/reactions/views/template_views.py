from django.http import JsonResponse
import json
from django.views.decorators.csrf import csrf_exempt
from reactions.models import ReactionTemplate, User


@csrf_exempt
def create_template(request):
    if request.method == 'POST':
        try:
            # Parse the form data
            form_data = request.POST
            user_id = form_data.get('userID')
            template_name = form_data.get('template_name')
            # Ensure user is authenticated
            if not user_id:
                return JsonResponse(
                    {'status': 'error', 'message': 'User not authenticated.'}, status=403)

            user = User.objects.get(id=user_id)

            # Check if template name already exists
            if (ReactionTemplate.objects.filter(name=template_name, user=user).exists() or
                    ReactionTemplate.objects.filter(is_default=True, name=template_name).exists()):
                return JsonResponse(
                    {'status': 'error',
                    'message': f'You already have a template named {template_name}.'},
                    status=400)

            # Collect the fields for the template
            substrates = form_data.getlist('substrates')
            substrates_types = form_data.getlist('substrates_type')
            subs_comps = form_data.getlist('subs_comps')
            subs_sch = form_data.getlist('subs_sch')

            products = form_data.getlist('products')
            products_types = form_data.getlist('products_type')
            prods_comps = form_data.getlist('prod_comps')
            prods_sch = form_data.getlist('prod_sch')

            direction = form_data.get('direction', 'forward')
            subsystem = form_data.get('subsystem', 'undefined')
            organs = json.loads(form_data.get('organs', '[]'))

            description = form_data.get('description', '')
            # Create and save the new template
            template = ReactionTemplate.objects.create(
                name=template_name,
                user=user,
                is_default=False,
                description=description,
                substrates=','.join(substrates),
                substrates_types=','.join(substrates_types),
                subs_comps=','.join(subs_comps),
                subs_sch=','.join(subs_sch),
                products=','.join(products),
                products_types=','.join(products_types),
                prods_comps=','.join(prods_comps),
                prods_sch=','.join(prods_sch),
                direction=direction,
                subsystem=subsystem,
                Organs=','.join(organs)
            )

            user.templates.add(template)

            return JsonResponse({'status': 'success',
                                 'message': 'Template created successfully.'})
        except Exception as e:
            return JsonResponse(
                {'status': 'error', 'message': str(e)}, status=500)
    else:
        return JsonResponse(
            {'status': 'error', 'message': 'Invalid request method.'}, status=400)


@csrf_exempt
def get_rxn_template(request):
    if request.method == 'POST':
        try:
            # Parse the request body
            data = json.loads(request.body)
            reaction_type = data.get('reaction_type', '')
            user_id = data.get('userID')

            # First, try to get the template created by the user
            template = ReactionTemplate.objects.filter(name=reaction_type, user=user_id).first()

            # If no user-specific template exists, get the first available one
            if not template:
                template = ReactionTemplate.objects.filter(name=reaction_type).first()

            # If no template is found, return an error
            if not template:
                return JsonResponse({'error': 'Template not found'}, status=404)

            # Prepare the response data
            result = {
                'substrates': template.substrates.split(',') if template.substrates else [],
                'subs_sch': template.subs_sch.split(',') if template.subs_sch else ['1'] * len(template.substrates.split(',')),
                'subs_comps': template.subs_comps.split(',') if template.subs_comps else ['-'] * len(template.substrates.split(',')),
                'subs_types': template.substrates_types.split(',') if template.substrates_types else ['vmh'] * len(template.substrates.split(',')),
                'products': template.products.split(',') if template.products else [],
                'prod_sch': template.prods_sch.split(',') if template.prods_sch else ['1'] * len(template.products.split(',')),
                'prods_comps': template.prods_comps.split(',') if template.prods_comps else ['-'] * len(template.products.split(',')),
                'prods_types': template.products_types.split(',') if template.products_types else ['vmh'] * len(template.products.split(',')),
                'direction': template.direction,
                'subsystem': template.subsystem,
                'description': template.description,
                'name': template.name,
                'Organs': template.Organs.split(',') if template.Organs else []
            }
            return JsonResponse(result)

        except ReactionTemplate.DoesNotExist:
            return JsonResponse({'error': 'Template not found'}, status=404)
        except KeyError:
            return JsonResponse(
                {'error': 'Invalid template structure'}, status=400)
    else:
        return JsonResponse({'error': 'Invalid request method'}, status=400)


@csrf_exempt
def list_templates(request):
    if request.method == 'POST':
        default_templates = ReactionTemplate.objects.filter(
            is_default=True).values_list('name', flat=True)

        user_id = json.loads(request.body).get('userID')

        if user_id:
            user = User.objects.get(id=user_id)
            user_templates = user.templates.values_list('name', flat=True)
            templates = list(user_templates) + list(default_templates)
        else:
            templates = list(default_templates)
        templates.sort()  # Optionally sort the combined list

        return JsonResponse({'templates': templates})
    else:
        return JsonResponse(
            {'status': 'error', 'message': 'Invalid request method.'}, status=400)


@csrf_exempt
def share_template(request):
    if request.method == 'POST':
        data = json.loads(request.body)
        user_id = data.get('userID')
        share_with_username = data.get('share_with_user')
        template_names = data.get('template_names', [])

        if not user_id or not share_with_username or not template_names:
            return JsonResponse({'status': 'error', 'message': 'Missing required fields.'}, status=400)
        try:
            sharing_user = User.objects.get(id=user_id)
        except User.DoesNotExist:
            return JsonResponse({'status': 'error', 'message': 'Sharing user not found.'}, status=404)

        try:
            recipient_user = User.objects.get(name=share_with_username)
        except User.DoesNotExist:
            return JsonResponse({'status': 'error', 'message': 'Recipient user not found.'}, status=404)

        # Get the templates to share from the sharing user's saved templates.
        templates_to_share = sharing_user.templates.filter(name__in=template_names)

        for template in templates_to_share:
            # Determine a unique name for the recipient’s copy
            if recipient_user.templates.filter(name=template.name).exists():
                base_name = template.name
                counter = 1
                new_name = f"{base_name}{counter}"
                while recipient_user.templates.filter(name=new_name).exists():
                    counter += 1
                    new_name = f"{base_name}{counter}"
            else:
                new_name = template.name

            # Always create a copy for the recipient
            new_template = ReactionTemplate.objects.create(
                name=new_name,
                user=recipient_user,
                is_default=False,
                substrates=template.substrates,
                products=template.products,
                direction=template.direction,
                substrates_types=template.substrates_types,
                products_types=template.products_types,
                subsystem=template.subsystem,
                subs_comps=template.subs_comps,
                prods_comps=template.prods_comps,
                subs_sch=template.subs_sch,
                prods_sch=template.prods_sch,
                Organs=template.Organs,
                description=template.description,
            )
            new_template.save()
            recipient_user.templates.add(new_template)
            recipient_user.save()
        return JsonResponse({'status': 'success', 'message': 'Templates shared successfully.'}, status=200)
    else:
        return JsonResponse({'status': 'error', 'message': 'Invalid request method.'}, status=400)


@csrf_exempt
def update_template(request):
    if request.method == 'POST':
        try:
            data = json.loads(request.body)
            old_name = data.get('old_name')
            new_name = data.get('new_name')
            user_id = data.get('userID')
            description = data.get('description', '')

            if not user_id or not old_name or not new_name:
                return JsonResponse({'status': 'error', 'message': 'Missing required fields.'}, status=400)

            user = User.objects.get(id=user_id)

            # First, try to find the template in the user's saved templates.
            template = user.templates.filter(name=old_name).first()

            # If found and the user owns it, update it directly.
            if template and template.user and template.user.id == user.id:
                # Check that new_name is not already used in the user's list (personal + default)
                if (user.templates.filter(name=new_name).exclude(id=template.id).exists() or
                        ReactionTemplate.objects.filter(is_default=True, name=new_name).exists()):
                    return JsonResponse({'status': 'error',
                                         'message': f'A template named {new_name} already exists.'},
                                        status=400)
                template.name = new_name
                template.description = description
                template.save()
                return JsonResponse({'status': 'success', 'message': 'Template updated successfully.'})

            # Otherwise, the template was not owned by the user (could be default or from another user).
            # Look for the default template.
            template = ReactionTemplate.objects.filter(name=old_name).first()
            if template.is_default:
                return JsonResponse({'status': 'error', 'message': 'Cannot update default templates.'}, status=403)
            if not template:
                return JsonResponse({'status': 'error', 'message': 'Template not found in your list.'}, status=404)

            # Ensure the new name does not conflict with any template in the user's combined list.
            if (user.templates.filter(name=new_name).exists() or
                    ReactionTemplate.objects.filter(is_default=True, name=new_name).exists()):
                return JsonResponse({'status': 'error',
                                     'message': f'A template named {new_name} already exists.'},
                                    status=400)

            # Create a copy for the user and update it.
            template.pk = None  # This makes a new record when saved.
            template.user = user
            template.name = new_name
            template.description = description
            template.is_default = False
            template.save()
            user.templates.add(template)
            return JsonResponse({'status': 'success', 'message': 'Template copied and updated successfully.'})

        except User.DoesNotExist:
            return JsonResponse({'status': 'error', 'message': 'User not found.'}, status=404)
        except Exception as e:
            return JsonResponse({'status': 'error', 'message': str(e)}, status=500)

    return JsonResponse({'status': 'error', 'message': 'Invalid request method.'}, status=400)


@csrf_exempt
def delete_template(request):
    if request.method == 'POST':
        try:
            data = json.loads(request.body)
            template_name = data.get('template_name')
            user_id = data.get('userID')
            if not user_id:
                return JsonResponse({'status': 'error', 'message': 'Login required.'}, status=403)
            user = User.objects.get(id=user_id)

            # Find the template within the user's saved templates.

            template = user.templates.filter(name=template_name).first()
            if not template:
                template = ReactionTemplate.objects.filter(name=template_name, is_default=True).first()
                if template:
                    return JsonResponse({'status': 'error', 'message': 'Cannot delete default templates.'}, status=403)
                return JsonResponse({'status': 'error', 'message': 'Template not found in your list.'}, status=404)
            # Remove the template from the user's list.
            user.templates.remove(template)
            # Optionally, if no other user is associated with this template, delete it from the database.
            if template.created_by_users.count() == 0:
                template.delete()
            return JsonResponse({'status': 'success', 'message': 'Template deleted successfully.'})

        except User.DoesNotExist:
            return JsonResponse({'status': 'error', 'message': 'User not found.'}, status=404)
        except Exception as e:
            return JsonResponse({'status': 'error', 'message': str(e)}, status=500)

    return JsonResponse({'status': 'error', 'message': 'Invalid request method.'}, status=400)
