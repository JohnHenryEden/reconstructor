# Generated by Django 5.0.1 on 2024-08-25 15:07

from django.db import migrations, models


class Migration(migrations.Migration):

    dependencies = [
        ('reactions', '0007_rename_formulas_reaction_rxn_formula_abbr'),
    ]

    operations = [
        migrations.RemoveField(
            model_name='reaction',
            name='rxn_formula_abbr',
        ),
        migrations.AddField(
            model_name='reaction',
            name='rxn_formula',
            field=models.CharField(blank=True, max_length=10, null=True),
        ),
    ]
