import platform

if platform.system() == 'Windows':
    import os
    for path in os.environ['Path'].split(';')[:-1]:
        try:
            os.add_dll_directory(path)
        except (FileNotFoundError, OSError) as error:
            pass


from . import ion_crystals
from . import ms_pulse_shaping
