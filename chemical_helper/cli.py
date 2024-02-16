import click

from .translator import is_canonical, to_canonical_smiles, to_inchi


@click.group()
def cli():
    pass


@cli.command(name="is-canonical")
@click.option("--smiles", required=True, type=str, help="The SMILES string to check.")
def is_canonical_cli(smiles: str):
    """
    Command-line function to check if a given SMILES string is canonical.
    Args:
    smiles (str): a string representing a SMILES molecule

    Outputs whether the SMILES is canonical or not.
    """
    try:
        if is_canonical(smiles):
            click.echo(click.style("The SMILES is canonical.", fg="green"))
        else:
            click.echo(click.style("The SMILES is not canonical.", fg="yellow"))
    except Exception as e:
        click.echo(click.style(f"Error: {e}", fg="red"))


if __name__ == "__main__":
    cli()


@cli.command(name="to-canonical-smiles")
@click.option("--smiles", required=True, type=str, help="The SMILES string.")
@click.option("--isomeric", is_flag=True, help="Include isomeric information in the canonical SMILES.")
def to_canonical_smiles_cli(smiles: str, isomeric: bool):
    """
    Command-line function that generates a canonical SMILES from any input SMILES, with an option to include isomeric information.
    Args:
    smiles (str): a string representing a SMILES molecule
    isomeric (bool): flag to include isomeric information in the canonical SMILES

    Outputs the Canonical SMILES or Error message.
    """
    try:
        canonical_smiles_str = to_canonical_smiles(smiles, isomericSmiles=isomeric)
        click.echo(click.style(f"{canonical_smiles_str}", fg="green"))
    except Exception as e:
        click.echo(click.style(f"Error: {e}", fg="red"))


@cli.command(name="to-inchi")
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
        click.echo(click.style(f"Error: {e}", fg="red"))
