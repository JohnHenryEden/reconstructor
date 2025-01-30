from reactions.models import User, Reaction, Flag
from django.http import JsonResponse
import json


def get_user(request):
    username = request.POST.get('username', '')
    password = request.POST.get('password', '')
    try:
        user = User.objects.get(name=username)
        if user.check_password(password):
            # Store the user ID in the session
            request.session['userID'] = user.pk
            return JsonResponse(
                {'status': 'success', 'userName': user.name, 'userID': user.pk})
        else:
            return JsonResponse(
                {'status': 'error', 'message': 'Invalid password'})
    except User.DoesNotExist:
        return JsonResponse(
            {'status': 'error', 'message': 'User does not exist'})


def register_user(request):
    username = request.POST.get('username', '')
    password = request.POST.get('password', '')
    email = request.POST.get('email', '')
    orchid_id = request.POST.get('orchid_id', '')
    full_name = request.POST.get('full_name', '')   # New field
    users = User.objects.all()
    usernames = [user.name for user in users]
    if username in usernames:
        return JsonResponse(
            {'status': 'error', 'message': 'Username already exists'})
    if not username or not password or not email or not full_name:
        return JsonResponse(
            {'status': 'error', 'message': 'Ysername, password, email, and full name are required'})

    try:
        if User.objects.filter(name=username).exists():
            return JsonResponse(
                {'status': 'error', 'message': 'Username already exists'})

        user = User(
            name=username,
            email=email,
            password=password,  # Hash the password before saving
            orchid_id=orchid_id,
            full_name=full_name
        )
        user.save()

        request.session['userID'] = user.pk  # Store the user ID in the session
        return JsonResponse(
            {'status': 'success', 'userName': user.name, 'userID': user.pk})
    except Exception as e:
        return JsonResponse({'status': 'error', 'message': str(e)})


def validate_user_ID(user_id):
    """Validate the user ID."""
    try:
        user = User.objects.get(pk=user_id)
        return user
    except User.DoesNotExist:
        return None


def set_session_user(request):
    req_body = json.loads(request.body)
    userID = int(req_body.get('userID'))
    request.session['userID'] = userID
    if request.session.get('userID') == userID and validate_user_ID(userID):
        return JsonResponse(
            {'status': 'success', 'message': 'User set in session'})
    else:
        return JsonResponse(
            {'status': 'error', 'message': 'User not set in session'})


def save_user_reaction(request):
    if request.method == 'POST':
        reaction_id = request.POST.get('reaction_id')
        userID = request.POST.get('userID')
        short_name = request.POST.get('short_name')
        flag_name = request.POST.get('flag_name')
        flag_color = request.POST.get('flag_color')

        try:
            reaction = Reaction.objects.get(pk=reaction_id)
        except Reaction.DoesNotExist:
            return JsonResponse(
                {'status': 'error', 'message': 'Invalid reaction'})

        user = validate_user_ID(userID)

        if user and reaction:
            reaction.short_name = short_name
            if flag_name.strip() not in [
                    'None', 'Choose a flag'] and flag_color != 'null':
                # Get or create the flag
                flag = Flag.objects.get(name_flag=flag_name, color=flag_color)

                # Associate the new flag with the reaction
                reaction.flags.add(flag)

            reaction.save()
            user.saved_reactions.add(reaction)
            return JsonResponse({'status': 'success'})

        return JsonResponse(
            {'status': 'error', 'message': 'Invalid user or reaction'})

    return JsonResponse({'status': 'error', 'message': 'Invalid request'})
