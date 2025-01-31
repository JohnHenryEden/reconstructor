import json
import os
import uuid
from django.conf import settings
from django.core.management.base import BaseCommand
from reactions.models import Reaction
from reactions.utils.to_mol import any_to_mol
from rdkit.Chem import MolFromMolFile, MolToInchiKey
from tqdm import tqdm  

class Command(BaseCommand):
    help = "Generate and store reaction_signature for all existing Reaction objects"

    def handle(self, *args, **kwargs):
        reactions_updated = 0
        failed_reactions = 0
        total_reactions = Reaction.objects.count()

        with tqdm(total=total_reactions, desc="Processing Reactions") as pbar:
            for reaction in Reaction.objects.all():
                pbar.update(1)  

                if reaction.reaction_signature:
                    continue  

                try:
                    # Debug: Print the reaction ID before processing
                    self.stdout.write(f"🔍 Processing Reaction ID: {reaction.id}")

                    # Load stored data with extra validation
                    substrates_list = self.safe_json_load(reaction.substrates, "substrates", reaction.id)
                    products_list = self.safe_json_load(reaction.products, "products", reaction.id)
                    subs_sch = self.safe_json_load(reaction.subs_sch, "subs_sch", reaction.id)
                    prods_sch = self.safe_json_load(reaction.prods_sch, "prods_sch", reaction.id)
                    substrates_types = self.safe_json_load(reaction.substrates_types, "substrates_types", reaction.id)
                    products_types = self.safe_json_load(reaction.products_types, "products_types", reaction.id)

                    if any(v is None for v in [substrates_list, products_list, subs_sch, prods_sch, substrates_types, products_types]):
                        failed_reactions += 1
                        continue

                    # Process molecules
                    subs_mols = self.get_mol_objects(substrates_list, substrates_types)
                    prod_mols = self.get_mol_objects(products_list, products_types)

                    if any(mol is None for mol in subs_mols + prod_mols):
                        failed_reactions += 1
                        pbar.set_postfix(missing_mols=failed_reactions)

                    subs_inchi_keys = [(MolToInchiKey(mol), int(stoich)) for mol, stoich in zip(subs_mols, subs_sch) if mol]
                    prods_inchi_keys = [(MolToInchiKey(mol), int(stoich)) for mol, stoich in zip(prod_mols, prods_sch) if mol]

                    subs_inchi_keys.sort()
                    prods_inchi_keys.sort()

                    reaction_signature = json.dumps({
                        "substrates": subs_inchi_keys,
                        "products": prods_inchi_keys
                    }, sort_keys=True)

                    reaction.reaction_signature = reaction_signature
                    reaction.save()
                    reactions_updated += 1

                except Exception as e:
                    self.stderr.write(f"⚠️ Error processing reaction {reaction.id}: {str(e)}")

        self.stdout.write(f"\n✅ Updated {reactions_updated} reactions with reaction_signature.")
        self.stdout.write(f"⚠️ {failed_reactions} reactions had missing or invalid Mol objects.")

    def safe_json_load(self, json_str, field_name, reaction_id):
        """ Safely loads JSON while handling empty or malformed fields. """
        try:
            if not json_str or json_str.strip() == "":
                self.stderr.write(f"⚠️ Warning: Empty {field_name} for reaction {reaction_id}. Skipping.")
                return None
            return json.loads(json_str)
        except json.JSONDecodeError as e:
            self.stderr.write(f"❌ JSON Error in field {field_name} for reaction {reaction_id}: {str(e)}")
            return None

    def get_mol_objects(self, mols_list, mols_types):
        mol_objects = []
        for i in range(len(mols_list)):
            if mols_types[i] == "MDL Mol file":
                mol_file_path = os.path.join(settings.MEDIA_ROOT, mols_list[i].replace(settings.MEDIA_URL, ""))
                if os.path.exists(mol_file_path):
                    mol = MolFromMolFile(mol_file_path)
                else:
                    mol = None  
            else:
                mol, _, _ = any_to_mol([mols_list[i]], [mols_types[i]], None, side='unknown')
                mol = mol[0] if mol else None
            
            mol_objects.append(mol)
        return mol_objects
