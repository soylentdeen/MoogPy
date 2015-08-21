import os
import glob

def parse_config(filename, defaults=None):
    """ Read options for a configuration file.
        
        filename    the configuration filename
        defaults    a dictionary of allowed options with default values.
        
        The options are returned as a dictionary than can also be indexed
        by attribute.
        
        Sample format:
        
            lines = lines.dat
            x = 20
            save = T

            
        Ignores blank lines or lines starting with '#'
        
        Attempts to convert the values to int, float, boolean then string.
        

        Borrowed from Neil Crighton
    """
    if defaults is None:
        cfg = dict()
    else:
        cfg = dict(**defaults)
    fh = open(filename)
    for row in fh:
        if not row.strip() or row.lstrip().startswith('#'):
            continue
        option, value = [r.strip() for r in row.split('=')]
        if defaults is not None:
            assert option in defaults, 'Unknown option %r' % option
        try:
            value = int(value)
        except ValueError:
            try:
                value = float(value)
            except ValueError:
                if value == 'True':
                    value = True
                elif value == 'False':
                    value = False

        cfg[option] = value
    fh.close()
    return cfg

