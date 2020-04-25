from numbers import Integral
import numpy as np
from xml.etree import ElementTree as ET

import openmc.checkvalue as cv
from .mixin import EqualityMixin, IDManagerMixin


class Sensitivity(EqualityMixin, IDManagerMixin):
    """A material perturbation derivative to apply to a tally.

    Parameters
    ----------
    derivative_id : int, optional
        Unique identifier for the tally derivative. If none is specified, an
        identifier will automatically be assigned
    variable : str, optional
        Accepted values are 'density', 'nuclide_density', and 'temperature'
    material : int, optional
        The perturbed material ID
    nuclide : str, optional
        The perturbed nuclide. Only needed for 'nuclide_density' derivatives.
        Ex: 'Xe135'

    Attributes
    ----------
    id : int
        Unique identifier for the tally derivative
    variable : str
        Accepted values are 'density', 'nuclide_density', and 'temperature'
    material : int
        The perturubed material ID
    nuclide : str
        The perturbed nuclide. Only needed for 'nuclide_density' derivatives.
        Ex: 'Xe135'

    """

    next_id = 1
    used_ids = set()

    def __init__(self, sensitivity_id=None, variable=None, material=None,
                 nuclide=None, reaction=None,energy=None):
        # Initialize Tally class attributes
        self.id = sensitivity_id
        self.variable = variable
        self.nuclide = nuclide
        self.reaction = reaction
        self.energy = np.asarray(energy)

    def __repr__(self):
        string = 'Sensitivity\n'
        string += '{: <16}=\t{}\n'.format('\tID', self.id)
        string += '{: <16}=\t{}\n'.format('\tVariable', self.variable)
        string += '{: <16}=\t{}\n'.format('\tNuclide', self.nuclide)
        if self.variable == 'cross_section':
            string += '{: <16}=\t{}\n'.format('\tReaction', self.reaction)

        return string

    @property
    def variable(self):
        return self._variable

    @property
    def nuclide(self):
        return self._nuclide

    @variable.setter
    def variable(self, var):
        if var is not None:
            cv.check_type('sensitivity variable', var, str)
            cv.check_value('sensitivity variable', var,
                           ('cross_section', 'multipole'))
        self._variable = var

    @nuclide.setter
    def nuclide(self, nuc):
        if nuc is not None:
            cv.check_type('sensitivity nuclide', nuc, str)
        self._nuclide = nuc

    def to_xml_element(self):
        """Return XML representation of the tally derivative

        Returns
        -------
        element : xml.etree.ElementTree.Element
            XML element containing derivative data

        """

        element = ET.Element("sensitivity")
        element.set("id", str(self.id))

        subelement = ET.SubElement(element, 'nuclide')
        subelement.text = self.nuclide

        subelement = ET.SubElement(element, 'variable')
        subelement.text = self.variable
        
        if self.variable == 'cross_section':
            subelement = ET.SubElement(element, 'reaction')
            subelement.text = self.reaction

            subelement = ET.SubElement(element, 'energy')
            subelement.text = ' '.join(str(e) for e in self.energy)


        return element