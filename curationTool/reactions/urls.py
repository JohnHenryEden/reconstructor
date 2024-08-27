from django.contrib import admin
from django.urls import path
from . import views
from django.conf import settings
from reactions.views import saved_reactions  # Correctly import from the app where the views.py is located

from django.conf.urls.static import static
urlpatterns = [
    path('flags/<int:user_id>/', views.get_user_flags, name='get_user_flags'),
    path('add_flag/', views.add_flag, name='add_user_flag'),
    path('', views.input_reaction, name='input_reaction'),
    path('saved_reactions/', views.saved_reactions, name='saved_reactions'),
    path('get_ai_response/', views.get_ai_response, name='get_ai_response'),
    path('check_reaction_vmh/', views.check_reaction_vmh, name='check_reaction_vmh'),
    path('get_from_vmh/', views.get_from_vmh, name='get_from_vmh'),
    path('save_reaction/', views.save_user_reaction, name='save_reaction'),
    path('delete_reaction/', views.delete_reaction, name='delete_reaction'),
    path('chemdoodle_sketcher/', views.chemdoodle_sketcher, name='chemdoodle_sketcher'),
    path('add_info_to_reaction/', views.add_info_to_reaction, name='add_info_to_reaction'),
    path('get_reaction_details/', views.get_reaction_details, name='get_reaction_details'),
    path('add_to_vmh/', views.add_to_vmh, name='add_to_vmh'),
    path('get_user/', views.get_user, name='get_user'),
    path('get_subsystems/', views.get_subsystems, name='get_subsystems'),
    path('set_session_user/', views.set_session_user, name='set_session_user'),
    path('verify_metabolite/', views.verify_metabolite, name='verify_metabolite'),
    path('prepare_add_to_vmh/', views.prepare_add_to_vmh, name='prepare_add_to_vmh'),
    path('get_pubmed_info/<str:pmid>/', views.get_pubmed_info, name='get_pubmed_info'),
    path('get_doi_info/<path:doi>/', views.get_doi_info, name='get_doi_info'),
    path('get_gene_info/', views.get_gene_info, name='get_gene_info'),
    path('update_subsystems/', views.update_subsystems, name='update_subsystems'),
    path('delete_reaction_info/', views.delete_reaction_info, name='delete_reaction_info'),
    path('gene_parsing/', views.gene_parsing, name='gene_parsing'),
    path('register_user',views.register_user,name='register_user'),
    path('check-reaction', views.is_reaction_in_user_saved_reactions, name='check_reaction'),
    path('available_reactions', views.get_available_reactions, name='available_reactions'),
    path('reaction_view', views.reaction_view, name='reaction_view'),
    path('check-session/', views.check_session_data, name='check_session'),
    path('delete-gene-info/', views.delete_gene_info_from_session, name='delete_gene_info_from_session'),
    path('clear-session/', views.clear_session, name='clear_session'),
    path('create-formula-abbr/', views.create_formula_abbr, name='create_formula_abbr'),
    path('save-formula/', views.save_formula, name='save_formula'),

    path('saved_reactions/', views.saved_reactions_page, name='saved_reactions_page'),
    path('search_reactions/', views.search_reactions, name='search_reactions'),
    path('user-reactions-vmh/', views.get_user_reactions_and_vmh, name='get_user_reactions_and_vmh'),
    path('get_reaction/<int:reaction_id>/', views.get_reaction, name='get_reaction'),
    path('stats/', views.leader_board, name='leader_board'),
    path('create-reaction/', views.create_reaction, name='create_reaction'),
    path('parse_formula_with_compartments/', views.parse_formula_with_compartments, name='parse_formula_with_compartments'),
    path("gene_details_view/", views.gene_details_view, name="gene_details_view"),
    path("parse_gene_info/", views.parse_gene_info, name="parse_gene_info"),
    path('update_gene_info/', views.update_gene_info, name='update_gene_info'),
    path('get_rxn_template/', views.get_rxn_template, name='get_rxn_template'),
    path('update_reaction_data', views.update_reaction_data, name='update_reaction_data'),
    path('check_subsystem_similarity/', views.check_subsystem_similarity, name='check_subsystem_similarity'),
    path('temp_gene_details/',views.temp_gene_details,name='temp_gene_details'),
    path('convert_to_smiles/', views.convert_to_smiles, name='convert_to_smiles'),
    path('saved-reactions-modal/', lambda request: saved_reactions(request, modal=True), name='saved_reactions_modal_content'),
    path('saved_reactions_view/', views.saved_reactions_view, name='saved_reactions'),


    ]+ static(settings.MEDIA_URL, document_root=settings.MEDIA_ROOT)
