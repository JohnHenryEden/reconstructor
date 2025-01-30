from reactions.models import Flag, User, Reaction
import json
from django.http import JsonResponse
from django.views.decorators.csrf import csrf_exempt
from django.shortcuts import get_object_or_404


def get_user_flags(request, user_id):
    try:
        user = User.objects.get(pk=user_id)
        flags = user.flags.all()
        flags_data = [{'id': flag.id,
                       'name_flag': flag.name_flag,
                       'color': flag.color} for flag in flags]
        return JsonResponse({'status': 'success', 'flags': flags_data})
    except Exception as e:
        return JsonResponse({'status': 'error', 'message': str(e)}, status=500)


@csrf_exempt
def add_flag(request):
    if request.method == 'POST':
        try:
            data = json.loads(request.body)
            user_id = data.get('user_id')
            flag_name = data.get('name_flag')
            flag_color = data.get('color')
            user = User.objects.get(pk=user_id)
        except User.DoesNotExist:
            return JsonResponse({'status': 'error', 'message': 'Invalid user'})

        if flag_name and flag_color:
            flag_name = flag_name.strip()
            flag, created = Flag.objects.get_or_create(
                name_flag=flag_name, color=flag_color, user=user)
            return JsonResponse({
                'status': 'success',
                'message': 'Flag added successfully',
                'flag': {'id': flag.id, 'name_flag': flag.name_flag, 'color': flag.color}
            })
        else:
            return JsonResponse({'status': 'error',
                                 'message': 'Flag name and color are required'})
    return JsonResponse(
        {'status': 'error', 'message': 'Invalid request method'})


@csrf_exempt
def save_flags_in_saved_reactions(request):
    if request.method == 'POST':
        try:
            data = json.loads(request.body)
            user_id = data.get('userID')
            reaction_ids = data.get('reaction_ids', [])
            flag_name = data.get('flag_name')
            flag_color = data.get('flag_color')

            user = User.objects.get(pk=user_id)
            flag = Flag.objects.get(
                name_flag=flag_name, color=flag_color, user=user)
            for reaction_id in reaction_ids:
                reaction = get_object_or_404(Reaction, pk=reaction_id)
                reaction.flags.add(flag)

            return JsonResponse({'status': 'success'})
        except Exception as e:
            return JsonResponse({'status': 'error', 'message': str(e)})
    return JsonResponse(
        {'status': 'error', 'message': 'Invalid request method'})
