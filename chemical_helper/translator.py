from rdkit import Chem


class CustomError(Exception):
    pass


def is_canonical_smiles(input_smiles: str) -> bool:
    """
    Checks whether a SMILES string is canonical.

    Args:
    input_smiles (str): A string representing a SMILES molecule.

    Returns:
    bool: True if the SMILES string is canonical, False otherwise.
    """
    mol = Chem.MolFromSmiles(input_smiles)
    if not mol:
        # If RDKit cannot parse the SMILES, it's considered not canonical.
        # Alternatively, you could raise an exception or handle the error differently.
        return False

    canonical_smiles = Chem.MolToSmiles(mol, isomericSmiles=False, canonical=True)
    return input_smiles == canonical_smiles


def to_canonical_smiles(input_smiles: str, isomericSmiles: bool = False) -> str:
    """
    This function generates a canonical SMILES from any SMILES, with an option to include isomeric information.

    Args:
    input_smiles (str): A string representing a SMILES molecule.
    isomericSmiles (bool, optional): Whether to include isomeric information in the canonical SMILES. Defaults to False.

    Raises:
    CustomError: If the SMILES is not valid.

    Returns:
    str: A string that represents the canonical SMILES.
    """
    mol = Chem.MolFromSmiles(input_smiles)
    if not mol:
        error_message = "Invalid SMILES."
        raise CustomError(error_message)
    canonical_smiles = Chem.MolToSmiles(mol, isomericSmiles=isomericSmiles, canonical=True)
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
