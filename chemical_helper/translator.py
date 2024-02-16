from rdkit import Chem


class CustomError(Exception):
    pass


def to_cansmiles(input_smiles: str) -> str:
    """
    This function generates a canonical SMILES from any SMILES.

    Args:
    input_smiles (str): a string representing a SMILES molecule

    Raises:
    ValueError: If the SMILES is not valid

    Returns:
    str: a string that represents the canonical SMILES
    """
    mol = Chem.MolFromSmiles(input_smiles)
    if not mol:
        error_message = "Invalid SMILES."
        raise CustomError(error_message)
    canonical_smiles = Chem.MolToSmiles(mol, isomericSmiles=True)
    return canonical_smiles


def to_inchi(canonical_smiles: str) -> str:
    """
    This function generates an InChI from a Canonical SMILES.

    Args:
    canonical_smiles (str): a string representing a Canonical SMILES molecule

    Raises:
    ValueError: If the Canonical SMILES is not valid

    Returns:
    str: a string that represents the InChI of the molecule
    """
    mol = Chem.MolFromSmiles(canonical_smiles)
    if not mol:
        error_message = "Invalid canonical SMILES."
        raise CustomError(error_message)
    inchi = Chem.MolToInchi(mol)
    return inchi
