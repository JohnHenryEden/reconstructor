# Generated by Django 5.0.1 on 2025-01-31 18:45

from django.db import migrations, models


class Migration(migrations.Migration):

    dependencies = [
        ("reactions", "0015_reaction_reaction_signature"),
    ]

    operations = [
        migrations.AddField(
            model_name="reaction",
            name="description",
            field=models.TextField(blank=True, null=True),
        ),
    ]
