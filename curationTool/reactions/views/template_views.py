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
            if ReactionTemplate.objects.filter(
                    name=template_name, user=user).exists():
                return JsonResponse(
                    {
                        'status': 'error',
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

            # Create and save the new template
            template = ReactionTemplate.objects.create(
                name=template_name,
                user=user,
                is_default=False,
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

            # Fetch the reaction template by name
            template = ReactionTemplate.objects.get(name=reaction_type)

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
    print(request.method)
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
        # Parse the request body
        data = json.loads(request.body)
        # The ID of the user sharing the templates
        user_id = data.get('userID')
        # The username of the recipient
        share_with_username = data.get('share_with_user')
        # List of template names to share
        template_names = data.get('template_names', [])

        # Ensure all required fields are present
        if not user_id or not share_with_username or not template_names:
            return JsonResponse(
                {'status': 'error', 'message': 'Missing required fields.'}, status=400)

        # Fetch the sharing user
        try:
            sharing_user = User.objects.get(id=user_id)
        except User.DoesNotExist:
            return JsonResponse(
                {'status': 'error', 'message': 'Sharing user not found.'}, status=404)

        # Fetch the recipient user by username
        try:
            recipient_user = User.objects.get(name=share_with_username)
        except User.DoesNotExist:
            return JsonResponse(
                {'status': 'error', 'message': 'Recipient user not found.'}, status=404)

        # Get templates from user object instead
        templates_to_share = sharing_user.templates.filter(
            name__in=template_names)

        # Handle conflicts and add templates
        for template in templates_to_share:
            if template in recipient_user.templates.all():
                return JsonResponse(
                    {
                        'status': 'error',
                        'message': f'Template {template.name} already exists for the recipient.'},
                    status=400)
            # Check if the recipient already has a template with the same name
            if ReactionTemplate.objects.filter(
                    user=recipient_user, name=template.name).exists():
                # Generate a unique name for the new template
                base_name = template.name
                counter = 1
                new_name = f"{base_name}{counter}"
                while ReactionTemplate.objects.filter(
                        user=recipient_user, name=new_name).exists():
                    counter += 1
                    new_name = f"{base_name}{counter}"

                # Create a new template with the modified name for the
                # recipient
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
                )
                recipient_user.templates.add(new_template)
            else:
                # Add the original template to the recipient if there's no
                # conflict
                recipient_user.templates.add(template)

        return JsonResponse(
            {'status': 'success', 'message': 'Templates shared successfully.'}, status=200)

    return JsonResponse(
        {'status': 'error', 'message': 'Invalid request method.'}, status=400)


@csrf_exempt
def rename_template(request):
    if request.method == 'POST':
        try:
            data = json.loads(request.body)
            old_name = data.get('old_name')
            new_name = data.get('new_name')
            userID = data.get('userID')
            if not userID:
                return JsonResponse(
                    {'status': 'error', 'message': 'Login required.'}, status=403)
            user = User.objects.get(id=userID)
            user_templates = User.objects.get(id=userID).templates

            if not old_name or not new_name:
                return JsonResponse(
                    {'status': 'error', 'message': 'Missing required fields.'}, status=400)

            template = ReactionTemplate.objects.get(name=old_name)
            user_templates.remove(template)
            template.pk = None

            # Check if a template with the new name already exists
            if ReactionTemplate.objects.filter(
                    name=new_name, user=template.user).exists():
                return JsonResponse(
                    {'status': 'error', 'message': f'A template named {new_name} already exists.'}, status=400)
            # Rename the template
            template.name = new_name
            template.save()
            # Add the renamed template
            user_templates.add(template)
            user.save()
            return JsonResponse({'status': 'success',
                                 'message': 'Template renamed successfully.'})
        except ReactionTemplate.DoesNotExist:
            return JsonResponse(
                {'status': 'error', 'message': 'Template not found.'}, status=404)
        except Exception as e:
            return JsonResponse(
                {'status': 'error', 'message': str(e)}, status=500)

    return JsonResponse(
        {'status': 'error', 'message': 'Invalid request method.'}, status=400)


@csrf_exempt
def delete_template(request):
    if request.method == 'POST':
        try:
            data = json.loads(request.body)
            template_name = data.get('template_name')
            userID = data.get('userID')
            user = User.objects.get(id=userID)
            if not userID:
                return JsonResponse(
                    {'status': 'error', 'message': 'Login required.'}, status=403)
            user_templates = User.objects.get(id=userID).templates
            # remove the template from the user's templates
            user_templates.remove(
                ReactionTemplate.objects.get(
                    name=template_name))
            user.save()
            return JsonResponse({'status': 'success',
                                 'message': 'Template deleted successfully.'})
        except ReactionTemplate.DoesNotExist:
            return JsonResponse(
                {'status': 'error', 'message': 'Template not found.'}, status=404)
        except Exception as e:
            return JsonResponse(
                {'status': 'error', 'message': str(e)}, status=500)

    return JsonResponse(
        {'status': 'error', 'message': 'Invalid request method.'}, status=400)
