# Licensed under a 3-clause BSD style license - see LICENSE.rst
from __future__ import (absolute_import, division, print_function,
                        unicode_literals)
import json
import os
from astropy.utils.data import download_file
from astropy.io import ascii
from astropy.table import QTable
from astropy.coordinates import SkyCoord
import astropy.units as u

__all__ = ['NasaExoplanetArchive']

EXOPLANETS_CSV_URL = ('http://exoplanetarchive.ipac.caltech.edu/cgi-bin/'
                      'nstedAPI/nph-nstedAPI?table=exoplanets')
KOIS_CSV_URL = ('http://exoplanetarchive.ipac.caltech.edu/cgi-bin/'
                      'nstedAPI/nph-nstedAPI?table=q1_q17_dr25_koi')
STELLAR_CSV_URL = ('http://exoplanetarchive.ipac.caltech.edu/cgi-bin/'
                      'nstedAPI/nph-nstedAPI?table=keplerstellar')
TIME_ATTRS = {'rowupdate': 'iso'}
BOOL_ATTRS = ('pl_kepflag', 'pl_ttvflag', 'pl_k2flag', 'st_massblend',
              'st_optmagblend', 'st_radblend', 'st_teffblend')


class NasaExoplanetArchiveClass(object):
    """
    Exoplanet Archive querying object. Use the ``get_confirmed_planets_table``
    or ``query_planet`` methods to get information about exoplanets via the NASA
    Exoplanet Archive.
    """
    def __init__(self):
        self._param_units = None
        self._table = None

    @property
    def param_units(self):
        if self._param_units is None:
            module_dir = os.path.dirname(os.path.abspath(__file__))
            units_file = open(os.path.join(module_dir, 'data',
                                           'exoplanet_nexsci_units.json'))
            self._param_units = json.load(units_file)

        return self._param_units

    def get_table(self, url=None, cache=True, show_progress=True,
                                    table_path=None):
        """
        Download (and optionally cache) a table from the `NExScI Exoplanet Archive 
                                    <http://exoplanetarchive.ipac.caltech.edu/index.html>`_.
                                    
        Parameters
        ----------
        url : str (optional)
            Web location of the table to be downloaded.
        cache : bool (optional)
            Cache table to local astropy cache? Default is `True`.
        show_progress : bool (optional)
            Show progress of table download (if no cached copy is
            available). Default is `True`.
        table_path : str (optional)
            Path to a local table file. Default `None` will trigger a
            download of the table from the internet.
        Returns
        -------
        table : `~astropy.table.QTable`
            Astropy table of requested data.
        """
        if table_path is None:
            table_path = download_file(url, cache=cache,
                                       show_progress=show_progress,
                                       timeout=120)
        table = ascii.read(table_path)

        # Create sky coordinate mixin column
        table['sky_coord'] = SkyCoord(ra=table['ra_str'], dec=table['dec_str'],
                                        unit=(u.hourangle, u.deg))

        # Assign units to columns where possible
        for col in table.colnames:
            if col in self.param_units:
                # Check that unit is implemented in this version of astropy
                if hasattr(u, self.param_units[col]):
                    table[col].unit = u.Unit(self.param_units[col])

        return QTable(table)

        
    def get_confirmed_planets_table(self, cache=True, show_progress=True,
                                    table_path=None):
        """
        Download (and optionally cache) the `NExScI Exoplanet Archive Confirmed
        Planets table <http://exoplanetarchive.ipac.caltech.edu/index.html>`_.

        The Exoplanet Archive table returns lots of columns of data. A full
        description of the columns can be found `here
        <https://exoplanetarchive.ipac.caltech.edu/docs/API_exoplanet_columns.html>`_

        Parameters
        ----------
        cache : bool (optional)
            Cache exoplanet table to local astropy cache? Default is `True`.
        show_progress : bool (optional)
            Show progress of exoplanet table download (if no cached copy is
            available). Default is `True`.
        table_path : str (optional)
            Path to a local table file. Default `None` will trigger a
            download of the table from the internet.
        Returns
        -------
        table : `~astropy.table.QTable`
            Table of exoplanet properties.
        """
        exoplanet_table = self.get_table(url=EXOPLANETS_CSV_URL, cache=cache, 
                                        show_progress=show_progress, 
                                        table_path=table_path)
        
        # Store column of lowercase names for indexing:
        lowercase_names = [host_name.lower().replace(' ', '') + letter
                           for host_name, letter in
                           zip(exoplanet_table['pl_hostname'].data,
                               exoplanet_table['pl_letter'].data)]
        exoplanet_table['NAME_LOWERCASE'] = lowercase_names
        exoplanet_table.add_index('NAME_LOWERCASE')
        
        self._exoplanet_table = exoplanet_table
        return self._exoplanet_table
        
    def get_kois_table(self, cache=True, show_progress=True,
                                    table_path=None):
        """
        Download (and optionally cache) the `NExScI Exoplanet Archive Kepler
        Objects of Interest table <http://exoplanetarchive.ipac.caltech.edu/index.html>`_.

        The Exoplanet Archive table returns lots of columns of data. A full
        description of the columns can be found `here
        <https://exoplanetarchive.ipac.caltech.edu/docs/API_kepcandidate_columns.html>`_

        Parameters
        ----------
        cache : bool (optional)
            Cache exoplanet table to local astropy cache? Default is `True`.
        show_progress : bool (optional)
            Show progress of exoplanet table download (if no cached copy is
            available). Default is `True`.
        table_path : str (optional)
            Path to a local table file. Default `None` will trigger a
            download of the table from the internet.
        Returns
        -------
        table : `~astropy.table.QTable`
            Table of KOI properties.
        """
        koi_table = self.get_table(url=KOIS_CSV_URL, cache=cache, 
                                        show_progress=show_progress, 
                                        table_path=table_path)
        
        self._koi_table = koi_table
        return self._koi_table
        
    def get_keplerstellar_table(self, cache=True, show_progress=True,
                                    table_path=None):
        """
        Download (and optionally cache) the `NExScI Exoplanet Archive Kepler
        Stellar table <http://exoplanetarchive.ipac.caltech.edu/index.html>`_.

        The Exoplanet Archive table returns lots of columns of data. A full
        description of the columns can be found `here
        <https://exoplanetarchive.ipac.caltech.edu/docs/API_keplerstellar_columns.html>`_

        Parameters
        ----------
        cache : bool (optional)
            Cache exoplanet table to local astropy cache? Default is `True`.
        show_progress : bool (optional)
            Show progress of exoplanet table download (if no cached copy is
            available). Default is `True`.
        table_path : str (optional)
            Path to a local table file. Default `None` will trigger a
            download of the table from the internet.
        Returns
        -------
        table : `~astropy.table.QTable`
            Table of Kepler stellar properties.
        """
        keplerstellar_table = self.get_table(url=STELLAR_CSV_URL, cache=cache, 
                                        show_progress=show_progress, 
                                        table_path=table_path)
        
        self._keplerstellar_table = keplerstellar_table
        return self._keplerstellar_table

                                    

    def query_planet(self, planet_name, table_path=None):
        """
        Get table of exoplanet properties.

        Parameters
        ----------
        planet_name : str
            Name of planet
        table_path : str (optional)
            Path to a local table file. Default `None` will trigger a
            download of the table from the internet.
        Return
        ------
        table : `~astropy.table.QTable`
            Table of one exoplanet's properties.
        """

        exoplanet_table = self.get_confirmed_planets_table(table_path=table_path)
        return exoplanet_table.loc[planet_name.strip().lower().replace(' ', '')]
        
    


NasaExoplanetArchive = NasaExoplanetArchiveClass()
