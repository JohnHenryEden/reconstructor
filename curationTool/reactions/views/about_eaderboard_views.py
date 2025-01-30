from django.shortcuts import render


def leader_board(request):
    return render(request, 'reactions/leader_board.html')


def about_view(request):
    return render(request, 'reactions/about.html')
