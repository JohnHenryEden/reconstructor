from django.core.management.base import BaseCommand
from reactions.models import ReactionTemplate
import pandas as pd
import re

class Command(BaseCommand):
    """
    Management command to load all default reaction templates into the database.
    """
    help = "Load all default reaction templates into the database"

    def handle(self, *args, **kwargs):
        DEFAULT_TEMPLATES = {
            'hydrolysis': {
                'substrates': 'empty,h2o,h',
                'products': 'empty,empty,h',
                'direction': 'forward',
                'substrates_types': 'vmh,vmh,vmh',
                'products_types': 'vmh,vmh,vmh',
                'subsystem': 'undefined',
                'subs_comps': 'c,c,c',
                'prods_comps': 'c,c,c',
                'subs_sch': '1,1,1',
                'prods_sch': '1,1,1',
            },
            'O2NADPHOX': {
                'substrates': 'empty,o2,nadph,h',
                'products': 'empty,h2o,nadp',
                'direction': 'forward',
                'substrates_types': 'vmh,vmh,vmh,vmh',
                'products_types': 'vmh,vmh,vmh',
                'subsystem': 'undefined',
                'subs_comps': 'c,c,c,c',
                'prods_comps': 'c,c,c',
                'subs_sch': '1,1,1,1',
                'prods_sch': '1,1,1',
            },
            'SULT': {
                'substrates': 'empty,paps',
                'products': 'empty,pap,h',
                'direction': 'forward',
                'substrates_types': 'vmh,vmh',
                'products_types': 'vmh,vmh,vmh',
                'subsystem': 'undefined',
                'subs_comps': 'c,c',
                'prods_comps': 'c,c,c',
                'subs_sch': '1,1',
                'prods_sch': '1,1,1',
            },
            'UGT': {
                'substrates': 'empty,udpglcur',
                'products': 'empty,udp,h',
                'direction': 'forward',
                'substrates_types': 'vmh,vmh',
                'products_types': 'vmh,vmh,vmh',
                'subsystem': 'undefined',
                'subs_comps': 'c,c',
                'prods_comps': 'c,c,c',
                'subs_sch': '1,1',
                'prods_sch': '1,1,1',
            },
            'UGT glucose': {
                'substrates': 'empty,udpg',
                'products': 'empty,udp,h',
                'direction': 'forward',
                'substrates_types': 'vmh,vmh',
                'products_types': 'vmh,vmh,vmh',
                'subsystem': 'undefined',
                'subs_comps': 'c,c',
                'prods_comps': 'c,c,c',
                'subs_sch': '1,1',
                'prods_sch': '1,1,1',
            },
            'UGT carb glucur': {
                'substrates': 'empty,udpglcur,co2',
                'products': 'empty,udp,h',
                'direction': 'forward',
                'substrates_types': 'vmh,vmh,vmh',
                'products_types': 'vmh,vmh,vmh',
                'subsystem': 'undefined',
                'subs_comps': 'c,c,c',
                'prods_comps': 'c,c,c',
                'subs_sch': '1,1,1',
                'prods_sch': '1,1,1',
            },
            'CoA': {
                'substrates': 'empty,coa,atp',
                'products': 'empty,amp,ppi',
                'direction': 'forward',
                'substrates_types': 'vmh,vmh,vmh',
                'products_types': 'vmh,vmh,vmh',
                'subsystem': 'undefined',
                'subs_comps': 'c,c,c',
                'prods_comps': 'c,c,c',
                'subs_sch': '1,1,1',
                'prods_sch': '1,1,1',
            },
            'FAOXhd': {
                'substrates': 'empty,h2o',
                'products': 'empty',
                'direction': 'forward',
                'substrates_types': 'vmh,vmh',
                'products_types': 'vmh',
                'subsystem': 'undefined',
                'subs_comps': 'c,c',
                'prods_comps': 'c',
                'subs_sch': '1,1',
                'prods_sch': '1',
            },
            'FAOXnad': {
                'substrates': 'empty,nad',
                'products': 'empty,nadh,h',
                'direction': 'forward',
                'substrates_types': 'vmh,vmh',
                'products_types': 'vmh,vmh,vmh',
                'subsystem': 'undefined',
                'subs_comps': 'c,c',
                'prods_comps': 'c,c,c',
                'subs_sch': '1,1',
                'prods_sch': '1,1,1',
            },
            'FAOXcoa': {
                'substrates': 'empty,coa',
                'products': 'empty,accoa',
                'direction': 'forward',
                'substrates_types': 'vmh,vmh',
                'products_types': 'vmh,vmh',
                'subsystem': 'undefined',
                'subs_comps': 'c,c',
                'prods_comps': 'c,c',
                'subs_sch': '1,1',
                'prods_sch': '1,1',
            },
            'Nad ox': {
                'substrates': 'empty,nad',
                'products': 'empty,nadh,h',
                'direction': 'forward',
                'substrates_types': 'vmh,vmh',
                'products_types': 'vmh,vmh,vmh',
                'subsystem': 'undefined',
                'subs_comps': 'c,c',
                'prods_comps': 'c,c,c',
                'subs_sch': '1,1',
                'prods_sch': '1,1,1',
            },
            'NADH red': {
                'substrates': 'empty,nadh,h',
                'products': 'empty,nad',
                'direction': 'forward',
                'substrates_types': 'vmh,vmh,vmh',
                'products_types': 'vmh,vmh',
                'subsystem': 'undefined',
                'subs_comps': 'c,c,c',
                'prods_comps': 'c,c',
                'subs_sch': '1,1,1',
                'prods_sch': '1,1',
            },
            'cycl': {
                'substrates': 'empty,h',
                'products': 'empty,h2o',
                'direction': 'forward',
                'substrates_types': 'vmh,vmh',
                'products_types': 'vmh,vmh',
                'subsystem': 'undefined',
                'subs_comps': 'c,c',
                'prods_comps': 'c,c',
                'subs_sch': '1,1',
                'prods_sch': '1,1',
            },
            'NADPH red': {
                'substrates': 'empty,nadph,h',
                'products': 'empty,nadp',
                'direction': 'forward',
                'substrates_types': 'vmh,vmh,vmh',
                'products_types': 'vmh,vmh',
                'subsystem': 'undefined',
                'subs_comps': 'c,c,c',
                'prods_comps': 'c,c',
                'subs_sch': '1,1,1',
                'prods_sch': '1,1',
            },
            'AT': {
                'substrates': 'taur,empty',
                'products': 'empty,coa,h',
                'direction': 'forward',
                'substrates_types': 'vmh,vmh',
                'products_types': 'vmh,vmh,vmh',
                'subsystem': 'undefined',
                'subs_comps': 'c,c',
                'prods_comps': 'c,c,c',
                'subs_sch': '1,1',
                'prods_sch': '1,1,1',
            },
        }

        try:
            df = pd.read_excel('/home/saleh/Downloads/Rxn Recap.xlsx')
        except Exception as e:
            self.stdout.write(self.style.WARNING(
                f"Spreadsheet not found or error loading spreadsheet: {e}"
            ))
            df = None

        # Create a dictionary of templates from the spreadsheet.
        # We key by the template name.
        spreadsheet_templates = {}
        if df is not None:
            for idx, row in df.iterrows():
                template_name = str(row.get('Template', '')).strip()
                if not template_name:
                    continue

                formula = str(row.get('Formula VMH IDs', '')).strip()
                if not formula:
                    continue

                # Parse the formula to extract fields.
                parsed_template = self.parse_reaction_formula(formula)
                # Add the description from the spreadsheet ("Rxn Name")
                parsed_template['description'] = str(row.get('Rxn Name', '')).strip()

                # If duplicate exists, update the default template with the description.
                if template_name in spreadsheet_templates:
                    if parsed_template['description']:
                        spreadsheet_templates[template_name]['description'] = parsed_template['description']
                else:
                    spreadsheet_templates[template_name] = parsed_template

        # Merge: For duplicates, keep the DEFAULT_TEMPLATES entry (which may have extra fields)
        # and update its description with the one from the spreadsheet.
        for name, sp_template in spreadsheet_templates.items():
            if name in DEFAULT_TEMPLATES:
                DEFAULT_TEMPLATES[name]['description'] = sp_template.get('description', '')
            else:
                DEFAULT_TEMPLATES[name] = sp_template

        # Now update or create ReactionTemplate objects.
        for name, template_data in DEFAULT_TEMPLATES.items():
            ReactionTemplate.objects.update_or_create(
                name=name,
                is_default=True,
                user=None,
                defaults=template_data
            )

        self.stdout.write(self.style.SUCCESS('All default templates loaded successfully!'))

    def parse_reaction_formula(self, formula):
        """
        Parse a reaction formula (e.g., "[c] + paps[c] -> [c] + pap[c] + h[c]")
        to extract:
          - substrates and products names (with 'empty' for missing names or 'metID')
          - stoichiometry (defaults to '1' for each species)
          - compartment info (e.g., 'c' from "[c]")
          - reaction direction
        Returns a dictionary with keys:
          substrates, products, direction,
          substrates_types, products_types,
          subs_comps, prods_comps, subs_sch, prods_sch,
          and subsystem.
        """
        # Determine the reaction arrow and split sides.
        if '<=>' in formula:
            left_side, right_side = formula.split('<=>')
            direction = 'bidirectional'
        elif '->' in formula:
            left_side, right_side = formula.split('->')
            direction = 'forward'
        elif '<=' in formula:
            left_side, right_side = formula.split('<=')
            direction = 'reverse'
        else:
            return {}

        def parse_side(side):
            parts = [part.strip() for part in side.split('+') if part.strip()]
            names = []
            comps = []
            sch = []
            for part in parts:
                # Assume part is like "[c]" or "paps[c]" or "metID[c]".
                m = re.match(r'^(.*?)(\[[^\]]+\])$', part)
                if m:
                    name = m.group(1).strip().lower()
                    comp = m.group(2).strip().strip('[]')
                else:
                    name = part.strip().lower()
                    comp = 'c'
                if name == '' or name == 'metid':
                    name = 'empty'
                names.append(name)
                comps.append(comp)
                sch.append('1')
            return names, comps, sch

        subs, subs_comps, subs_sch = parse_side(left_side)
        prods, prods_comps, prods_sch = parse_side(right_side)

        return {
            'substrates': ','.join(subs),
            'products': ','.join(prods),
            'direction': direction,
            'substrates_types': ','.join(['vmh'] * len(subs)),
            'products_types': ','.join(['vmh'] * len(prods)),
            'subsystem': 'undefined',
            'subs_comps': ','.join(subs_comps),
            'prods_comps': ','.join(prods_comps),
            'subs_sch': ','.join(subs_sch),
            'prods_sch': ','.join(prods_sch),
        }