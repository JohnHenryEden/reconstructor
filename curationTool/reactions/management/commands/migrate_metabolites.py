import json
from django.core.management.base import BaseCommand
from django.db.models import Q
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem.Descriptors import MolWt
from rdkit.Chem.rdMolDescriptors import CalcMolFormula
from reactions.models import User, Reaction, SavedMetabolite
from tqdm import tqdm
from reactions.utils.to_mol import any_to_mol
class Command(BaseCommand):
    help = 'Migrate non-VMH metabolites in user reactions to SavedMetabolite objects'

    def handle(self, *args, **options):
        users = User.objects.all()
        total_users = users.count()
        
        # List to track issues
        issues = []

        # Outer progress bar for users
        with tqdm(total=total_users, desc="Processing Users", unit="user") as pbar_users:
            for user in users:
                self.stdout.write(f"\nProcessing user: {user.name}")
                reactions = user.saved_reactions.all()
                total_reactions = reactions.count()
                
                # Inner progress bar for reactions
                with tqdm(total=total_reactions, desc="Processing Reactions", unit="reaction", leave=False) as pbar_reactions:
                    for reaction in reactions:
                        try:
                            self.process_reaction(user, reaction, issues)
                        except Exception as e:
                            issues.append({
                                "reaction_id": reaction.id,
                                "user": user.name,
                                "metabolites": None,
                                "reason": str(e)
                            })
                        pbar_reactions.update(1)  # Update reaction progress
                        pbar_reactions.set_postfix({"User": user.id})  # Show current user ID
                    
                pbar_users.update(1)  # Update user progress
                pbar_users.set_postfix({"Reactions Processed": total_reactions})  # Show reactions processed for the user
        
        # After all users and reactions are processed, print the issues
        if issues:
            self.stdout.write("\n\nReactions with Issues:")
            for issue in issues:
                self.stdout.write(f"Reaction ID: {issue['reaction_id']}, User: {issue['user']}")
                if issue['metabolites']:
                    for metabolite in issue['metabolites']:
                        self.stdout.write(f"\tMetabolite: {metabolite['name']}, Reason: {metabolite['reason']}")
                else:
                    self.stdout.write(f"\tReason: {issue['reason']}")

    def process_reaction(self, user, reaction, issues):
        # Process substrates
        substrates = json.loads(reaction.substrates)
        substrates_types = json.loads(reaction.substrates_types)
        substrates_names = json.loads(reaction.substrates_names)
        updated = False

        for i in range(len(substrates)):
            mtype = substrates_types[i]
            if mtype in ['VMH', 'Saved']:
                continue
            original_id = substrates[i]
            mols, _, _ = any_to_mol(mols=[original_id], types=[mtype], request=None, side=None)
            mol = mols[0]
            if not mol:
                issues.append({
                    "reaction_id": reaction.id,
                    "user": user.name,
                    "metabolites": [{"name": substrates_names[i], "reason": f"Invalid {mtype}: {original_id}"}],
                    "reason": "Invalid metabolite"
                })
                continue
            inchi_key = Chem.MolToInchiKey(mol)
            existing = SavedMetabolite.objects.filter(
                Q(owner=user) | Q(shared_with=user),
                inchi_key=inchi_key
            ).first()
            if not existing:
                # Check name conflict
                name = substrates_names[i]
                if SavedMetabolite.objects.filter(Q(owner=user) | Q(shared_with=user), name=name.strip()).exists():
                    issues.append({
                        "reaction_id": reaction.id,
                        "user": user.name,
                        "metabolites": [{"name": name.strip(), "reason": "Name conflict"}],
                        "reason": "Name conflict"
                    })
                    continue
                # Create new metabolite (same as original code, but errors added for failed creations)
                try:
                    saved_metabolite = self.create_metabolite(user, name.strip(), mol, original_id, mtype)
                    existing = saved_metabolite
                except Exception as e:
                    issues.append({
                        "reaction_id": reaction.id,
                        "user": user.name,
                        "metabolites": [{"name": substrates_names[i], "reason": f"Failed to create metabolite: {str(e)}"}],
                        "reason": "Metabolite creation failed"
                    })
                    continue

            substrates[i] = str(existing.id)
            substrates_types[i] = 'Saved'
            updated = True

        # Process products (similar to substrates, repeating the same logic)
        products = json.loads(reaction.products)
        products_types = json.loads(reaction.products_types)
        products_names = json.loads(reaction.products_names)

        for i in range(len(products)):
            mtype = products_types[i]
            if mtype in ['VMH', 'Saved']:
                continue
            original_id = products[i]
            mols, _, _ = any_to_mol(mols=[original_id], types=[mtype], request=None, side=None)
            mol = mols[0]
            if not mol:
                issues.append({
                    "reaction_id": reaction.id,
                    "user": user.name,
                    "metabolites": [{"name": products_names[i], "reason": f"Invalid {mtype}: {original_id}"}],
                    "reason": "Invalid metabolite"
                })
                continue
            inchi_key = Chem.MolToInchiKey(mol)
            existing = SavedMetabolite.objects.filter(
                Q(owner=user) | Q(shared_with=user),
                inchi_key=inchi_key
            ).first()
            if not existing:
                name = products_names[i]
                if SavedMetabolite.objects.filter(Q(owner=user) | Q(shared_with=user), name=name.strip()).exists():
                    issues.append({
                        "reaction_id": reaction.id,
                        "user": user.name,
                        "metabolites": [{"name": name.strip(), "reason": "Name conflict"}],
                        "reason": "Name conflict"
                    })
                    continue
                try:
                    saved_metabolite = self.create_metabolite(user, name.strip(), mol, original_id, mtype)
                    existing = saved_metabolite
                except Exception as e:
                    issues.append({
                        "reaction_id": reaction.id,
                        "user": user.name,
                        "metabolites": [{"name": products_names[i], "reason": f"Failed to create metabolite: {str(e)}"}],
                        "reason": "Metabolite creation failed"
                    })
                    continue

            products[i] = str(existing.id)
            products_types[i] = 'Saved'
            updated = True

        if updated:
            # Save the updated reaction
            reaction.substrates = json.dumps(substrates)
            reaction.substrates_types = json.dumps(substrates_types)
            reaction.products = json.dumps(products)
            reaction.products_types = json.dumps(products_types)
            reaction.save()
            self.stdout.write(f"Updated reaction {reaction.id} for user {user.name}")

    def create_metabolite(self, user, name, mol, original_id, mtype):
        # Same metabolite creation logic as before, but moved to its own method to avoid duplication
        inchi_key = Chem.MolToInchiKey(mol)
        smiles = Chem.MolToSmiles(mol)
        inchi = Chem.MolToInchi(mol)
        mol_w = MolWt(mol)
        mol_formula = CalcMolFormula(mol)
        try:
            if AllChem.EmbedMolecule(mol, AllChem.ETKDG()) != 0:
                raise ValueError("Embedding failed.")
            AllChem.UFFOptimizeMolecule(mol)
            mol_file = Chem.MolToMolBlock(mol)
        except Exception as e:
            mol_file = Chem.MolToMolBlock(mol)
        saved_metabolite = SavedMetabolite.objects.create(
            owner=user,
            name=name,
            inchi_key=inchi_key,
            inchi=inchi,
            smiles=smiles,
            mol_w=mol_w,
            mol_formula=mol_formula,
            mol_file=mol_file,
            source_type=mtype,
            original_identifier=original_id
        )
        return saved_metabolite
