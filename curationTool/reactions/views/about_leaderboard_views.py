"""
This module provides views for rendering the leaderboard and about pages.
"""

from django.shortcuts import render


def leader_board(request):
    """
    Render the leaderboard page.

    Returns:
        HttpResponse: The rendered 'leader_board.html' template.
    """
    return render(request, 'reactions/leader_board.html')


def about_view(request):
    """
    Render the about page.

    Returns:
        HttpResponse: The rendered 'about.html' template.
    """
    return render(request, 'reactions/about.html')
    