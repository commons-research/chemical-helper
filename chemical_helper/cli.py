import click

from .translator import to_cansmiles, to_inchi


@click.group()
def cli():
    pass


@cli.command()
@click.option("--smiles", required=True, type=str, help="The SMILES string.")
def to_cansmiles_cli(smiles: str):
    """
    Command-line function that generates a canonical SMILES from any input SMILES.
    Args:
    smiles (str): a string representing a SMILES molecule

    Outputs the Canonical SMILES or Error message.
    """
    try:
        canonical_smiles_str = to_cansmiles(smiles)
        click.echo(click.style(f"{canonical_smiles_str}", fg="green"))
    except Exception as e:
        click.echo(click.style(f"Error: {e!s}", fg="red"))


@cli.command()
@click.option("--cansmiles", required=True, type=str, help="The Canonical SMILES string.")
def to_inchi_cli(cansmiles: str):
    """
    Command-line function that generates an InChI from an input Canonical SMILES.
    Args:
    cansmiles (str): a string representing a Canonical SMILES molecule

    Outputs the InChI or Error message.
    """
    try:
        inchi_str = to_inchi(cansmiles)
        click.echo(click.style(f"{inchi_str}", fg="green"))
    except Exception as e:
        click.echo(click.style(f"Error: {e!s}", fg="red"))


if __name__ == "__main__":
    cli()
