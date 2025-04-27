"""
 Descriptors derived from a molecule's 3D structure

"""
from __future__ import annotations
from rdkit.Chem.Descriptors import _isCallable
from rdkit.Chem import rdMolDescriptors
__all__ = ['CalcMolDescriptors3D', 'descList', 'rdMolDescriptors']
def CalcMolDescriptors3D(mol, confId = None):
    """
    
        Compute all 3D descriptors of a molecule
        
        Arguments:
        - mol: the molecule to work with
        - confId: conformer ID to work with. If not specified the default (-1) is used
        
        Return:
        
        dict
            A dictionary with decriptor names as keys and the descriptor values as values
    
        raises a ValueError 
            If the molecule does not have conformers
        
    """
def _setupDescriptors(namespace):
    ...
descList: list  # value = [('PMI1', <function <lambda> at 0x100a52cb0>), ('PMI2', <function <lambda> at 0x10387ad40>), ('PMI3', <function <lambda> at 0x10387add0>), ('NPR1', <function <lambda> at 0x10387ae60>), ('NPR2', <function <lambda> at 0x10387aef0>), ('RadiusOfGyration', <function <lambda> at 0x10387af80>), ('InertialShapeFactor', <function <lambda> at 0x10387b010>), ('Eccentricity', <function <lambda> at 0x10387b0a0>), ('Asphericity', <function <lambda> at 0x10387b130>), ('SpherocityIndex', <function <lambda> at 0x10387b1c0>), ('PBF', <function <lambda> at 0x10387b250>)]
