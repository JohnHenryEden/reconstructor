from rdkit.Chem.rdMolDescriptors import CalcMolFormula
import rdkit.Chem as Chem
from reactions.reaction_info import calculate_total_charge
from reactions.utils.to_mol import smiles_with_explicit_hydrogens

def get_mol_info(mols):
    formulas, charges, mol_file_strings = [], [], []
    stereo_counts, stereo_locations_list = [], []

    for mol in mols:
        mol.UpdatePropertyCache(strict=False)
        mol_formula = CalcMolFormula(mol)
        mol_charge = calculate_total_charge([mol])

        try:
            mol_smiles = Chem.MolToSmiles(mol, allHsExplicit=True)
            mol_smiles = smiles_with_explicit_hydrogens(mol_smiles)
            mol = Chem.MolFromSmiles(mol_smiles)
            mol = Chem.AddHs(mol, explicitOnly=True)
        except:
            pass

        # Convert to MolBlock
        mol_block = Chem.MolToMolBlock(mol)
        mol_file_strings.append(mol_block)
        formulas.append(mol_formula)
        charges.append(mol_charge)

        # Detect chiral centers (including unassigned centers)
        chiral_centers = Chem.FindMolChiralCenters(mol, includeUnassigned=True)
        # Filter only the centers labeled with '?', i.e., unspecified
        unspecified_centers = [center for center in chiral_centers if center[1] == '?']

        # Count how many centers are unspecified
        stereo_count = len(unspecified_centers)
        # Collect the atom indices where they occur
        stereo_locations = [center[0] for center in unspecified_centers]

        stereo_counts.append(stereo_count)
        stereo_locations_list.append(stereo_locations)

    # Return stereo info alongside existing data
    return formulas,charges,mol_file_strings,stereo_counts,stereo_locations_list
